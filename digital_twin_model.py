# digital_twin_model.py

import numpy as np
from scipy.integrate import solve_ivp
from collections import deque

from Alarms.alarmModule import alarmModule
import time, json
from datetime import datetime, timedelta

## pathologies
from Events.Pathologies.pathology import Pathology
## therapy
from Events.Therapy.therapies import Therapy

class DigitalTwinModel:
    def __init__(self, patient_id, param_file="sepsis.json", 
                 data_callback=None, 
                 sleep=True,
                 time_step=0.01
                 ): ## modify for brute-force jobs
        self.patient_id = patient_id
        self.running = False
        self.t = 0
        self.dt = time_step  # Time step
        window_sec      = 200                         # how many seconds of history
        self.window_n   = int(window_sec / self.dt)  # number of samples
        self.P_buffer   = deque(maxlen=self.window_n)
        self.print_interval = 5  # Interval for printing heart rate
        self.data_callback = data_callback  # Callback function to emit data
        self.output_frequency = 1  # Output frequency for data callback -> 1 Hz

        self.pathologies = Pathology()  # Initialize pathologie-events 
        self.therapeutic = Therapy()
        self.alarmModule = alarmModule()  # Initialize alarm module
        self.sleep = sleep  # Sleep between iterations, boolean

        self.events = [] ## Actionable event introduction

        ## set datapoint amount in local memory
        self.data_points = 120 ## 2 minutes 

        self.data_epoch = {}
        self.start_timestamp = datetime.now()

        # Load parameters from JSON
        self._load_parameters(param_file)

        # Initialize model parameters before state
        self.initialize_model_parameters()

        # Initialize state variables
        self.current_state = self.initialize_state()

        self.current_heart_rate = 0  # Initialize monitored value
        #self.master_parameters = {}  # Initialize master parameters

        # Initialize Inference class
        self._compute_all_derived_params()


        # ── A) helper to flatten master_parameters access ─────────────────────────
    def get(self, key: str) -> float:
        """Fetch .master_parameters[key]['value'] or raise KeyError."""
        try:
            return self.master_parameters[key]['value']
        except KeyError:
            raise KeyError(f"Parameter '{key}' not found in master_parameters")

    def set(self, key: str, val: float):
        """Write .master_parameters[key]['value']=val, creating entry if needed."""
        if key in self.master_parameters:
            self.master_parameters[key]['value'] = val
        else:
            self.master_parameters[key] = {'value': val}

    def _compute_all_derived_params(self):
        """Populate all the M_*, V_*, D_* and bloodflow fractions in master_parameters."""
        mp = self.master_parameters
        M_O2 = self.get('params.M_O2')

        # 1) Base production rates
        derived = {
            'params.M_B_CO2': 0.2 * M_O2,
            'params.M_CO2':   0.85 * M_O2,
            'params.M_B_O2': -0.2 * M_O2,
        }
        derived['params.M_S_CO2'] = derived['params.M_CO2'] - derived['params.M_B_CO2']
        derived['params.M_S_O2']  = -M_O2 - derived['params.M_B_O2']

        # 2) Volumes
        Vtot = {'CO2': self.get('params.V_CO2'),
                'O2':  self.get('params.V_O2')}
        Vbt  = {'CO2': self.get('params.V_Btis_CO2'),
                'O2':  self.get('params.V_Btis_O2')}
        fcap = self.get('params.f_V_cap')

        derived['params.V_Stis_CO2'] = Vtot['CO2'] - Vbt['CO2']
        derived['params.V_Stis_O2']  = Vtot['O2']  - Vbt['O2']
        derived['params.V_Bcap_CO2'] = fcap * Vbt['CO2']
        derived['params.V_Bcap_O2']  = fcap * Vbt['O2']
        derived['params.V_Scap_CO2'] = fcap * derived['params.V_Stis_CO2']
        derived['params.V_Scap_O2']  = fcap * derived['params.V_Stis_O2']

        # 3) Diffusion constants
        w    = self.get('params.w')
        K_CO2 = self.get('params.K_CO2')
        K_O2τ = self.get('params.K_O2_tau')
        derived['params.D_T_CO2'] = 9/60 * w / K_CO2
        derived['params.D_T_O2']  = 9/60 * w / K_O2τ

        # 4) Blood‐flow fractions
        q_p = self.get('bloodflows.q_p')
        derived['bloodflows.q_Bv'] = 0.2 * q_p
        derived['bloodflows.q_S']  = 0.8 * q_p

        # 5) initial capillary concentrations (if not present)
        def maybe_init(src_key, tgt_key, op):
            if tgt_key not in mp:
                c_src = self.get(src_key)
                val = op(c_src)
                derived[tgt_key] = val

        # c_Scap_CO2 = c_Stis_CO2 - M_S_CO2/D_T_CO2
        maybe_init('initial_conditions.c_Stis_CO2',
                   'initial_conditions.c_Scap_CO2',
                   lambda c: c - derived['params.M_S_CO2']/derived['params.D_T_CO2'])
        # c_Scap_O2
        maybe_init('initial_conditions.c_Stis_O2',
                   'initial_conditions.c_Scap_O2',
                   lambda c: c + derived['params.M_S_O2']/derived['params.D_T_O2'])
        # c_Bcap_CO2
        maybe_init('initial_conditions.c_Btis_CO2',
                   'initial_conditions.c_Bcap_CO2',
                   lambda c: c - derived['params.M_B_CO2']/derived['params.D_T_CO2'])
        # c_Bcap_O2
        maybe_init('initial_conditions.c_Btis_O2',
                   'initial_conditions.c_Bcap_O2',
                   lambda c: c + derived['params.M_B_O2']/derived['params.D_T_O2'])

        # 6) gas_exchange_params ⇒ just mirror D_T_*
        derived.update({
            'gas_exchange_params.D_S_CO2': derived['params.D_T_CO2'],
            'gas_exchange_params.D_B_CO2': derived['params.D_T_CO2'],
            'gas_exchange_params.D_S_O2':  derived['params.D_T_O2'],
            'gas_exchange_params.D_B_O2':  derived['params.D_T_O2'],
        })

        # Finally, write all derived back into master_parameters
        for key, val in derived.items():
            self.set(key, val)        

    def add_disease(self, disease, severity):
        package = {
            "disease": disease,
            "severity": severity}
        self.events.append(package)

    def _load_parameters(self, param_file):
        """
        Load parameters from a JSON configuration file.
        Raise exceptions with detailed messages for missing or invalid files.
        """
        try:
            print(f"Loading parameter file: {param_file}")
            with open(param_file, "r") as file:
                initialHealthyParams = json.load(file)
            self.master_parameters = initialHealthyParams
            
            print(f"Successfully loaded parameters for patient {self.patient_id}.")
        except FileNotFoundError:
            raise FileNotFoundError(f"Parameter file '{param_file}' not found. Please verify the file path.")
        except json.JSONDecodeError as e:
            raise ValueError(f"Parameter file '{param_file}' is not a valid JSON file: {e}")

        # Process parameters that involve expressions
        self.process_parameter_expressions()

    def process_parameter_expressions(self):
        """Evaluate expressions in parameters and update the dictionaries."""

        # Ensure 't_eval' is defined
        if 't_eval' not in self.master_parameters:
            self.master_parameters['misc_constants.t_eval'] = {"value": None, "min": None, "max": None}
            self.master_parameters['misc_constants.t_eval']['value'] = np.arange(
                self.master_parameters['misc_constants.tmin']['value'],
                self.master_parameters['misc_constants.tmax']['value'] + self.master_parameters['misc_constants.T']['value'],
                self.master_parameters['misc_constants.T']['value']
            )

    def _compute_all_derived_params(self):


        # --- 1. Base O2 / CO2 production rates ---
        M_O2   = self.master_parameters['params.M_O2']['value']                   # 5.2
        # M_B_CO2 = 0.2 * M_O2
        self.master_parameters['params.M_B_CO2'] = {'value': 0.2 * M_O2}
        # M_CO2 = 0.85 * M_O2
        self.master_parameters['params.M_CO2']   = {'value': 0.85 * M_O2}
        # M_S_CO2 = M_CO2 - M_B_CO2
        self.master_parameters['params.M_S_CO2'] = {
            'value': self.master_parameters['params.M_CO2']['value']
                   - self.master_parameters['params.M_B_CO2']['value']
        }
        # M_B_O2 = -0.2 * M_O2
        self.master_parameters['params.M_B_O2']  = {'value': -0.2 * M_O2}
        # M_S_O2 = -M_O2 - M_B_O2  (matches your old “-5.2 - M_B_O2”)
        self.master_parameters['params.M_S_O2']  = {
            'value': -M_O2 - self.master_parameters['params.M_B_O2']['value']
        }

        # --- 2. Volumes ---
        V_CO2      = self.master_parameters['params.V_CO2']['value']
        V_O2       = self.master_parameters['params.V_O2']['value']
        V_Btis_CO2 = self.master_parameters['params.V_Btis_CO2']['value']
        V_Btis_O2  = self.master_parameters['params.V_Btis_O2']['value']
        fVcap      = self.master_parameters['params.f_V_cap']['value']

        # V_Stis_* = total - blood-tissue
        self.master_parameters['params.V_Stis_CO2'] = {'value': V_CO2 - V_Btis_CO2}
        self.master_parameters['params.V_Stis_O2']  = {'value': V_O2  - V_Btis_O2 }

        # capillary volumes = f_V_cap * tissue volumes
        self.master_parameters['params.V_Bcap_CO2'] = {'value': fVcap * V_Btis_CO2}
        self.master_parameters['params.V_Bcap_O2']  = {'value': fVcap * V_Btis_O2 }
        self.master_parameters['params.V_Scap_CO2'] = {'value': fVcap * self.master_parameters['params.V_Stis_CO2']['value']}
        self.master_parameters['params.V_Scap_O2']  = {'value': fVcap * self.master_parameters['params.V_Stis_O2']['value']}

        # --- 3. Diffusion constants ---
        w       = self.master_parameters['params.w']['value']
        K_CO2   = self.master_parameters['params.K_CO2']['value']
        K_O2τ   = self.master_parameters['params.K_O2_tau']['value']
        # D_T_CO2 = 9/60 * w / K_CO2
        self.master_parameters['params.D_T_CO2'] = {'value': 9/60 * w / K_CO2}
        # D_T_O2  = 9/60 * w / K_O2_tau
        self.master_parameters['params.D_T_O2']  = {'value': 9/60 * w / K_O2τ}

        # --- 4. Blood-flow fractions ---
        q_p = self.master_parameters['bloodflows.q_p']['value']
        self.master_parameters['bloodflows.q_Bv'] = {'value': 0.2 * q_p}
        self.master_parameters['bloodflows.q_S']  = {'value': 0.8 * q_p}
    # --- 5. Derived initial concentrations ---
    # pull directly from self.master_parameters (mp)
        if 'initial_conditions.c_Scap_CO2' not in self.master_parameters:
            c_Stis_CO2 = self.master_parameters['initial_conditions.c_Stis_CO2']['value']
    
    def _compute_all_derived_params(self):
        """Populate all the M_*, V_*, D_* and bloodflow fractions in master_parameters."""
        mp = self.master_parameters
        M_O2 = self.get('params.M_O2')

        # 1) Base production rates
        derived = {
            'params.M_B_CO2': 0.2 * M_O2,
            'params.M_CO2':   0.85 * M_O2,
            'params.M_B_O2': -0.2 * M_O2,
        }
        derived['params.M_S_CO2'] = derived['params.M_CO2'] - derived['params.M_B_CO2']
        derived['params.M_S_O2']  = -M_O2 - derived['params.M_B_O2']

        # 2) Volumes
        Vtot = {'CO2': self.get('params.V_CO2'),
                'O2':  self.get('params.V_O2')}
        Vbt  = {'CO2': self.get('params.V_Btis_CO2'),
                'O2':  self.get('params.V_Btis_O2')}
        fcap = self.get('params.f_V_cap')

        derived['params.V_Stis_CO2'] = Vtot['CO2'] - Vbt['CO2']
        derived['params.V_Stis_O2']  = Vtot['O2']  - Vbt['O2']
        derived['params.V_Bcap_CO2'] = fcap * Vbt['CO2']
        derived['params.V_Bcap_O2']  = fcap * Vbt['O2']
        derived['params.V_Scap_CO2'] = fcap * derived['params.V_Stis_CO2']
        derived['params.V_Scap_O2']  = fcap * derived['params.V_Stis_O2']

        # 3) Diffusion constants
        w    = self.get('params.w')
        K_CO2 = self.get('params.K_CO2')
        K_O2τ = self.get('params.K_O2_tau')
        derived['params.D_T_CO2'] = 9/60 * w / K_CO2
        derived['params.D_T_O2']  = 9/60 * w / K_O2τ

        # 4) Blood‐flow fractions
        q_p = self.get('bloodflows.q_p')
        derived['bloodflows.q_Bv'] = 0.2 * q_p
        derived['bloodflows.q_S']  = 0.8 * q_p

        # 5) initial capillary concentrations (if not present)
        def maybe_init(src_key, tgt_key, op):
            if tgt_key not in mp:
                c_src = self.get(src_key)
                val = op(c_src)
                derived[tgt_key] = val

        # c_Scap_CO2 = c_Stis_CO2 - M_S_CO2/D_T_CO2
        maybe_init('initial_conditions.c_Stis_CO2',
                   'initial_conditions.c_Scap_CO2',
                   lambda c: c - derived['params.M_S_CO2']/derived['params.D_T_CO2'])
        # c_Scap_O2
        maybe_init('initial_conditions.c_Stis_O2',
                   'initial_conditions.c_Scap_O2',
                   lambda c: c + derived['params.M_S_O2']/derived['params.D_T_O2'])
        # c_Bcap_CO2
        maybe_init('initial_conditions.c_Btis_CO2',
                   'initial_conditions.c_Bcap_CO2',
                   lambda c: c - derived['params.M_B_CO2']/derived['params.D_T_CO2'])
        # c_Bcap_O2
        maybe_init('initial_conditions.c_Btis_O2',
                   'initial_conditions.c_Bcap_O2',
                   lambda c: c + derived['params.M_B_O2']/derived['params.D_T_O2'])

        # 6) gas_exchange_params ⇒ just mirror D_T_*
        derived.update({
            'gas_exchange_params.D_S_CO2': derived['params.D_T_CO2'],
            'gas_exchange_params.D_B_CO2': derived['params.D_T_CO2'],
            'gas_exchange_params.D_S_O2':  derived['params.D_T_O2'],
            'gas_exchange_params.D_B_O2':  derived['params.D_T_O2'],
        })

        # Finally, write all derived back into master_parameters
        for key, val in derived.items():
            self.set(key, val)

    def compute_cardiac_parameters(self):
        """Extract elastance‐min/max, resistance, uvolume arrays from master_parameters."""
        # we assume keys are sorted: elastance_min_0…N, elastance_max_0…N, resistance_0…N, uvolume_0…N
        elast_min = sorted(k for k in self.master_parameters if 'cardio' in k and 'elastance' in k and 'min' in k)
        elast_max = sorted(k for k in self.master_parameters if 'cardio' in k and 'elastance' in k and 'max' in k)
        resist    = sorted(k for k in self.master_parameters if 'cardio' in k and 'resistance' in k)
        uv        = sorted(k for k in self.master_parameters if 'cardio' in k and 'uvolume' in k)

        self.elastance = np.array([
            [self.get(k) for k in elast_min],
            [self.get(k) for k in elast_max],
        ])
        self.resistance = np.array([self.get(k) for k in resist])
        self.uvolume    = np.array([self.get(k) for k in uv])


    def initialize_model_parameters(self):
        """Initialize model parameters like elastance, resistance, uvolume, etc."""
        # load cardiac computed parameters (elastance, resistance, uVolme)
        self.compute_cardiac_parameters()

        # Lung model equations
        # Mechanical system parameters
        self.C_cw = 0.2445  # l/cmH2O

        self.A_mechanical = np.array([
            [-1 / (self.master_parameters['respi_constants.C_l']['value'] * ( 1.021 * 1.5 )) ## R_ml constant
                     - 1 / (self.master_parameters['respi_constants.R_lt']['value'] * self.master_parameters['respi_constants.C_l']['value']), 
                     1 / (self.master_parameters['respi_constants.R_lt']['value'] * self.master_parameters['respi_constants.C_l']['value']), 0, 0, 0],
            [1 / (self.master_parameters['respi_constants.R_lt']['value'] * self.master_parameters['respi_constants.C_tr']['value']),
                 -1 / (self.master_parameters['respi_constants.C_tr']['value'] * self.master_parameters['respi_constants.R_lt']['value']) - 1 / (self.master_parameters['respi_constants.R_tb']['value'] * self.master_parameters['respi_constants.C_tr']['value']), 
                 1 / (self.master_parameters['respi_constants.R_tb']['value'] * self.master_parameters['respi_constants.C_tr']['value']), 0, 0],
            [0, 1 / (self.master_parameters['respi_constants.R_tb']['value'] * self.master_parameters['respi_constants.C_b']['value']),
             -1 / (self.master_parameters['respi_constants.C_b']['value'] * self.master_parameters['respi_constants.R_tb']['value']) - 1 / (self.master_parameters['respi_constants.R_bA']['value'] * self.master_parameters['respi_constants.C_b']['value']), 
             1 / (self.master_parameters['respi_constants.R_bA']['value'] * self.master_parameters['respi_constants.C_b']['value']), 0],
            [0, 0, 1 / (self.master_parameters['respi_constants.R_bA']['value'] * self.master_parameters['respi_constants.C_A']['value']), -1 / (self.master_parameters['respi_constants.C_A']['value'] * self.master_parameters['respi_constants.R_bA']['value']), 0],
            [1 / (self.master_parameters['respi_constants.R_lt']['value'] * self.C_cw), -1 / (self.C_cw * self.master_parameters['respi_constants.R_lt']['value']), 0, 0, 0]
        ])

        self.B_mechanical = np.array([
            [1 / ( ( 1.021 * 1.5 ) ## R_ml constant
                  * self.master_parameters['respi_constants.C_l']['value']), 0, 0],
            [0, 1, 0],
            [0, 1, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])

        # Initialize heart period parameters
        self.HR = self.master_parameters['misc_constants.HR']['value']
        self.update_heart_period(self.HR)

        # Initialize buffer for sliding window
        self.dt = self.master_parameters['misc_constants.T']['value']
        window_duration = 20  # seconds
        self.window_size = int(window_duration / self.dt)
        self.P_store = np.zeros((10, self.window_size))
        self.F_store = np.zeros((10, self.window_size))
        self.HR_store = np.zeros(self.window_size)
        self.buffer_index = 0

    def initialize_state(self):
        """Set initial state variables based on loaded parameters."""
        # Initialize state variables for the combined model
        state = np.zeros(29)  # Adjust the size if needed
        # Initialize blood volumes based on unstressed volumes
        TBV = self.master_parameters['misc_constants.TBV']['value']
        state[:10] = TBV * (self.uvolume / np.sum(self.uvolume))
        # try to adjust the initial volumes
        state[0]=state[0]+200
        state[1]=state[1]+100
        state[2]=state[2]+100   
        # Initialize mechanical states (indices 10 to 14)
        state[10:15] = np.zeros(5)  # Adjust initial values if necessary
        # Initialize other states as before
        state[15] = 157 / 731  # FD_O2
        state[16] = 7 / 731    # FD_CO2
        state[17] = self.master_parameters['initial_conditions.p_a_CO2']['value']
        state[18] = self.master_parameters['initial_conditions.p_a_O2']['value']
        state[19] = self.master_parameters['initial_conditions.c_Stis_CO2']['value']

        self.M_B_CO2 = 0.2 * self.master_parameters['params.M_O2']['value']
        self.M_S_CO2 = (self.master_parameters['params.M_O2']['value'] *0.85 ) - self.M_B_CO2
        self.D_T_CO2 = 9 / 60 * self.master_parameters['params.w']['value'] / self.master_parameters['params.K_CO2']['value']
        state[20] = 0.543 - self.M_S_CO2 / self.D_T_CO2 

        state[21] = self.master_parameters['initial_conditions.c_Stis_O2']['value']
        self.D_T_O2 = 9 / 60 * self.master_parameters['params.w']['value'] / self.master_parameters['params.K_O2_tau']['value']
        self.M_S_O2 = -5.2 - (-0.2 * self.master_parameters['params.M_O2']['value']) 

        state[22] = 0.128 - self.M_S_O2 / self.D_T_O2
        state[23] = self.master_parameters['respiratory_control_params.Delta_RR_c']['value']
        state[24] = self.master_parameters['respiratory_control_params.Delta_Pmus_c']['value']
        state[25] = -2  # Pmus
        state[26] = 0   # Delta_HR_c
        state[27] = 0   # Delta_R_c
        state[28] = 0   # Delta_UV_c
        return state

    def update_heart_period(self, HR):
        """Calculate heart period parameters."""
        self.HR = HR
        self.HP = 60 / HR
        self.Tas = 0.03 + 0.09 * self.HP
        self.Tav = 0.01
        self.Tvs = 0.16 + 0.2 * self.HP

    def calculate_elastances(self, t):
        """
        Compute time-varying elastances for atria (era), ventricles (elv),
        and “active” elastance for arteries (ela) from a single normalized phase φ.
        """
        # 1) normalized cardiac phase φ ∈ [0,1)
        φ = (t % self.HP) / self.HP

        # 2) precompute fractional durations
        a_frac = self.Tas / self.HP              # active (systolic) phase end
        c_start = (self.Tas + self.Tav) / self.HP # diastolic start
        v_frac = self.Tvs / self.HP              # duration of ventricle contraction

        # 3) atrial/arterial “active” factor aaf
        if φ <= a_frac and a_frac > 0:
            aaf = np.sin(np.pi * (φ / a_frac))
        else:
            aaf = 0.0

        # 4) ventricular factor vaf
        #    φ ∈ (c_start, c_start+v_frac]
        if c_start < φ <= c_start + v_frac and v_frac > 0:
            vaf = np.sin(np.pi * ((φ - c_start) / v_frac))
        else:
            vaf = 0.0

        # 5) pull the four compartments by index
        #    0: min-elastance, 1: max-elastance
        #    [8] arterial, [4] atrial, [9] ventricular, [5] ventricular
        emin, emax = self.elastance
        ela = emin[8] + (emax[8] - emin[8]) * aaf
        era = emin[4] + (emax[4] - emin[4]) * aaf
        elv = emin[9] + (emax[9] - emin[9]) * vaf
        erv = emin[5] + (emax[5] - emin[5]) * vaf

        return ela, elv, era, erv


    def get_inputs(self, t):
        """Get inputs for the cardiovascular system."""
        ela, elv, era, erv = self.calculate_elastances(t)

        return ela, elv, era, erv

    def compute_variables(self, t, y):
        """Compute pressure, flow, and heart rate based on state variables."""
        V = y[:10]  # Volumes for cardiovascular model

        # Inputs for cardiovascular model
        ela, elv, era, erv = self.get_inputs(t)

        # Need to get x[25], which is Pmus
        Pmus = y[25]

        # Calculate the pressures for cardiovascular model
        P = np.zeros(10)
        P[0] = self.elastance[0, 0] * (V[0] - self.uvolume[0]) + Pmus
        P[1] = self.elastance[0, 1] * (V[1] - self.uvolume[1])
        UV_c = self.master_parameters['cardio_control_params.UV_n']['value'] + y[28]
        P[2] = self.elastance[0, 2] * (V[2] - self.uvolume[2] * UV_c)
        P[3] = self.elastance[0, 3] * (V[3] - self.uvolume[3] * UV_c) + Pmus
        P[4] = era * (V[4] - self.uvolume[4]) + Pmus
        P[5] = erv * (V[5] - self.uvolume[5]) + Pmus
        P[6] = self.elastance[0, 6] * (V[6] - self.uvolume[6]) + Pmus
        P[7] = self.elastance[0, 7] * (V[7] - self.uvolume[7]) + Pmus
        P[8] = ela * (V[8] - self.uvolume[8]) + Pmus
        P[9] = elv * (V[9] - self.uvolume[9]) + Pmus

        # Calculate the flows for cardiovascular model
        R_c = self.master_parameters['cardio_control_params.R_n']['value'] - y[27]  # Delta_R_c is y[27]
        F = np.zeros(10)
        F[0] = (P[0] - P[1]) / (self.resistance[0] * R_c)
        F[1] = (P[1] - P[2]) / (self.resistance[1] * R_c)
        F[2] = (P[2] - P[3]) / (self.resistance[2] * R_c)
        F[3] = (P[3] - P[4]) / self.resistance[3] if P[3] - P[4] > 0 else (P[3] - P[4]) / (10 * self.resistance[3])
        F[4] = (P[4] - P[5]) / self.resistance[4] if P[4] - P[5] > 0 else 0
        F[5] = (P[5] - P[6]) / self.resistance[5] if P[5] - P[6] > 0 else 0
        F[6] = (P[6] - P[7]) / self.resistance[6]
        F[7] = (P[7] - P[8]) / self.resistance[7] if P[7] - P[8] > 0 else (P[7] - P[8]) / (10 * self.resistance[7])
        F[8] = (P[8] - P[9]) / self.resistance[8] if P[8] - P[9] > 0 else 0
        F[9] = (P[9] - P[0]) / self.resistance[9] if P[9] - P[0] > 0 else 0

        # Compute heart rate
        HR = self.master_parameters['cardio_control_params.HR_n']['value'] - y[26]  # Delta_HR_c is y[26]

        # Compute SaO2
        p_a_O2 = y[18]
        CaO2 = (self.master_parameters['params.K_O2']['value'] * np.power((1 - np.exp(-self.master_parameters['params.k_O2']['value'] * min(p_a_O2, 700))), 2)) * 100
        Sa_O2 = np.round(((CaO2 - p_a_O2 * 0.003 / 100) / (self.master_parameters['misc_constants.Hgb']['value'] * 1.34)) * 100)
        return P, F, HR, Sa_O2

    def ventilator_pressure(self, t):
        """Ventilator pressure as a function of time (simple square wave for demonstration)."""
        RR = self.master_parameters['respiratory_control_params.RR_0']['value']  # Respiratory Rate (breaths per minute)
        PEEP = 5  # Positive End-Expiratory Pressure (cm H2O)
        T = 60 / RR  # period of one respiratory cycle in seconds
        IEratio = 1
        TI = T * IEratio / (1 + IEratio)
        cycle_time = t % T
        if 0 <= cycle_time <= TI:
            return (PEEP + 15) * 0.735  # Inhalation , 1 cmH2O = 0.735 mmHg
        else:
            return PEEP * 0.735  # Exhalation

    def input_function(self, t, RR, Pmus_min, IEratio=1.0):
        """Input function for mechanical states."""
        T = 60 / RR
        TI = T * IEratio / (1 + IEratio)
        TE = T - TI
        exp_time = TE / 5
        cycle_time = t % T

        if 0 <= cycle_time <= TI:
            dPmus_dt = 2 * (-Pmus_min / (TI * TE)) * cycle_time + (Pmus_min * T) / (TI * TE)
        else:
            dPmus_dt = -Pmus_min / (exp_time * (1 - np.exp(-TE / exp_time))) * np.exp(-(cycle_time - TI) / exp_time)

        return np.array([0, dPmus_dt])

    def extended_state_space_equations(self, t, x):
        """
        Combined ODE for cardiovascular + respiratory system,
        with reflex gains computed from a sliding MAP buffer.
        """
        # ── 1) Unpack state vector ────────────────────────────────────────────────
        V = x[0:10]                          # blood volumes
        mech = x[10:15]                      # lung mechanical states
        FD_O2, FD_CO2 = x[15], x[16]
        p_a_CO2, p_a_O2 = x[17], x[18]
        c_Stis_CO2, c_Scap_CO2, c_Stis_O2, c_Scap_O2 = x[19:23]
        Δ_RR_c, Δ_Pmus_c = x[23], x[24]
        Pmus = x[25]
        Δ_HR_c, Δ_R_c, Δ_UV_c = x[26], x[27], x[28]

        # ── 2) Pull “setpoint” parameters ─────────────────────────────────────────
        HR_n    = self.master_parameters['cardio_control_params.HR_n']['value']
        R_n     = self.master_parameters['cardio_control_params.R_n']['value']
        UV_n    = self.master_parameters['cardio_control_params.UV_n']['value']
        RR0     = self.master_parameters['respiratory_control_params.RR_0']['value']
        FI_O2   = self.master_parameters['gas_exchange_params.FI_O2']['value']
        FI_CO2  = self.master_parameters['gas_exchange_params.FI_CO2']['value']
        CO_nom  = self.master_parameters['bloodflows.CO']['value']
        shunt   = self.master_parameters['bloodflows.sh']['value']
        MV_mode = self.master_parameters['misc_constants.MV']['value']

        # Apply baroreflex deltas
        HR = HR_n   - Δ_HR_c
        R_c = R_n   - Δ_R_c
        UV_c = UV_n + Δ_UV_c
        RR  = RR0  + Δ_RR_c

        # Update heart/respiratory rates
        self.HR, self.RR = HR, RR
        self.update_heart_period(HR)

        # ── 3) Ventilator vs. Spontaneous Breathing ───────────────────────────────
        if MV_mode == 0:
            # spontaneous
            P_ao = 0
            Pmus_min = (
                self.master_parameters['respiratory_control_params.Pmus_0']['value']
                + Δ_Pmus_c
            )
            _, Pmus_dt = self.input_function(t, RR, Pmus_min)
        else:
            # ventilator
            P_ao = self.ventilator_pressure(t)
            Pmus_dt = 0

        # ── 4) Cardiovascular Pressures P_i ──────────────────────────────────────
        ela, elv, era, erv = self.get_inputs(t)
        P = np.zeros(10)
        P[0] = self.elastance[0,0] * (V[0] - self.uvolume[0]) + mech[4]
        P[1] = self.elastance[0,1] * (V[1] - self.uvolume[1])
        P[2] = self.elastance[0,2] * (V[2] - self.uvolume[2] * UV_c)
        P[3] = self.elastance[0,3] * (V[3] - self.uvolume[3] * UV_c) + mech[4]
        P[4] = era * (V[4] - self.uvolume[4]) + mech[4]
        P[5] = erv * (V[5] - self.uvolume[5]) + mech[4]
        P[6] = self.elastance[0,6] * (V[6] - self.uvolume[6]) + mech[4]
        P[7] = self.elastance[0,7] * (V[7] - self.uvolume[7]) + mech[4]
        P[8] = ela * (V[8] - self.uvolume[8]) + mech[4]
        P[9] = elv * (V[9] - self.uvolume[9]) + mech[4]

        # ── 5) Flows F_i with simple reversal handling ────────────────────────────
        F = np.zeros(10)
        for i in range(9):
            R_eff = self.resistance[i] * (R_c if i < 3 else 1)
            ΔP = P[i] - P[i+1]
            if ΔP > 0:
                F[i] = ΔP / R_eff
            else:
                # for valves (i==3 or 7) allow tiny reverse flow
                F[i] = ΔP / (10 * R_eff) if i in (3,7) else 0
        ΔP = P[9] - P[0]
        F[9] = ΔP / self.resistance[9] if ΔP > 0 else 0

        # ── 6) Update MAP buffer ──────────────────────────────────────────────────
        self.P_store[0, self.buffer_index] = P[0]
        self.P_buffer.append(P[0])
        self.buffer_index = (self.buffer_index + 1) % self.window_size

        # ── 7) Chemo‐ & Baroreflex updates via one call ──────────────────────────
        dΔ_RR_c, dΔ_Pmus_c, dΔ_HR_c, dΔ_R_c, dΔ_UV_c = self.apply_reflexes(
            p_a_CO2, Δ_RR_c, Δ_Pmus_c, Δ_HR_c, Δ_R_c, Δ_UV_c
        )

        # ── 8) Volume derivatives dV/dt ──────────────────────────────────────────
        dVdt = np.zeros(10)
        dVdt[0] = F[9] - F[0]
        for i in range(1,10):
            dVdt[i] = F[i-1] - F[i]

        # ── 9) Lung mechanics ─────────────────────────────────────────────────────
        R_lt = self.master_parameters['respi_constants.R_lt']['value']
        C_cw = self.C_cw
        Ppl_dt = mech[0]/(R_lt*C_cw) - mech[1]/(C_cw*R_lt) + Pmus_dt
        dxdt_mech = (
            self.A_mechanical.dot(mech)
            + self.B_mechanical.dot([P_ao, Ppl_dt, Pmus_dt])
        )

        # ── 10) Gas exchange & alveolar ODEs ──────────────────────────────────────
        R_ml = 1.021 * 1.5
        Vdot_l = (P_ao - mech[0]) / R_ml
        R_bA   = self.master_parameters['respi_constants.R_bA']['value']
        Vdot_A = (mech[2] - mech[3]) / R_bA

        p_D_CO2 = FD_CO2 * 713
        p_D_O2  = FD_O2  * 713
        FA_CO2  = p_a_CO2 / 713
        FA_O2   = p_a_O2  / 713

        K_CO2, k_CO2 = (
            self.master_parameters['params.K_CO2']['value'],
            self.master_parameters['params.k_CO2']['value']
        )
        c_a_CO2 = K_CO2 * p_a_CO2 + k_CO2

        if p_a_O2 > 700:
            p_a_O2 = 700
        K_O2, k_O2 = (
            self.master_parameters['params.K_O2']['value'],
            self.master_parameters['params.k_O2']['value']
        )
        c_a_O2 = K_O2 * (1 - np.exp(-k_O2 * p_a_O2))**2

        c_v_CO2, c_v_O2 = c_Scap_CO2, c_Scap_O2

        V_D = self.master_parameters['gas_exchange_params.V_D']['value']
        V_A = self.master_parameters['gas_exchange_params.V_A']['value']

        if (MV_mode == 1 and P_ao > 6 * 0.735) or (MV_mode == 0 and mech[0] < 0):
            dFD_O2_dt = Vdot_l * 1000 * (FI_O2  - FD_O2 ) / (V_D * 1000)
            dFD_CO2_dt= Vdot_l * 1000 * (FI_CO2 - FD_CO2) / (V_D * 1000)
            dp_a_CO2  = (
                863 * CO_nom * (1-shunt) * (c_v_CO2 - c_a_CO2)
                + Vdot_A * 1000 * (p_D_CO2 - p_a_CO2)
            ) / (V_A * 1000)
            dp_a_O2   = (
                863 * CO_nom * (1-shunt) * (c_v_O2 - c_a_O2)
                + Vdot_A * 1000 * (p_D_O2 - p_a_O2)
            ) / (V_A * 1000)
        else:
            dFD_O2_dt = Vdot_A * 1000 * (FD_O2  - FA_O2 ) / (V_D * 1000)
            dFD_CO2_dt= Vdot_A * 1000 * (FD_CO2 - FA_CO2) / (V_D * 1000)
            dp_a_CO2  = 863 * CO_nom * (1-shunt) * (c_v_CO2 - c_a_CO2) / (V_A * 1000)
            dp_a_O2   = 863 * CO_nom * (1-shunt) * (c_v_O2 - c_a_O2) / (V_A * 1000)

        # ── 11) Systemic tissue ODEs ───────────────────────────────────────────────
        M_S_CO2 = self.master_parameters['params.M_S_CO2']['value']
        D_S_CO2 = self.master_parameters['gas_exchange_params.D_S_CO2']['value']
        V_Stis_CO2 = self.master_parameters['params.V_Stis_CO2']['value']
        V_Scap_CO2 = self.master_parameters['params.V_Scap_CO2']['value']

        dc_Stis_CO2 = (M_S_CO2 - D_S_CO2 * (c_Stis_CO2 - c_Scap_CO2)) / V_Stis_CO2
        dc_Scap_CO2 = (
            CO_nom * 0.8 * (c_a_CO2 - c_Scap_CO2)
            + D_S_CO2 * (c_Stis_CO2 - c_Scap_CO2)
        ) / V_Scap_CO2

        M_S_O2 = self.master_parameters['params.M_S_O2']['value']
        D_S_O2 = self.master_parameters['gas_exchange_params.D_S_O2']['value']
        V_Stis_O2 = self.master_parameters['params.V_Stis_O2']['value']
        V_Scap_O2 = self.master_parameters['params.V_Scap_O2']['value']

        dc_Stis_O2 = (M_S_O2 - D_S_O2 * (c_Stis_O2 - c_Scap_O2)) / V_Stis_O2
        dc_Scap_O2 = (
            CO_nom * 0.8 * (c_a_O2 - c_Scap_O2)
            + D_S_O2 * (c_Stis_O2 - c_Scap_O2)
        ) / V_Scap_O2

        # ── 12) Pack dx/dt ─────────────────────────────────────────────────────────
        dxdt = np.concatenate([
            dVdt,
            dxdt_mech,
            [
                dFD_O2_dt, dFD_CO2_dt,
                dp_a_CO2, dp_a_O2,
                dc_Stis_CO2, dc_Scap_CO2,
                dc_Stis_O2, dc_Scap_O2,
                dΔ_RR_c, dΔ_Pmus_c,
                Pmus_dt,
                dΔ_HR_c, dΔ_R_c, dΔ_UV_c
            ]
        ])

        return dxdt



    def apply_reflexes(self, p_a_CO2, Delta_RR_c, Delta_Pmus_c, 
                             Delta_HR_c, Delta_R_c, Delta_UV_c):
        """uses the last self.window_n samples of arterial pressure (P[0])
           to compute mean MAP for baroreflex, and the current p_a_CO2
           for chemoreflex."""
        mp = self.master_parameters

        # 1) baroreflex on HR, R & UV
        if len(self.P_buffer)>0:
            recent_MAP = float(np.mean(self.P_buffer))
        else:
            recent_MAP = mp['cardio_control_params.ABP_n']['value']

        #print(recent_MAP)
        ABP_n   = mp['cardio_control_params.ABP_n']['value']
        Gc_hr   = mp['cardio_control_params.Gc_hr']['value']
        tau_hr  = mp['cardio_control_params.tau_hr']['value']
        Gc_r    = mp['cardio_control_params.Gc_r']['value']
        tau_r   = mp['cardio_control_params.tau_r']['value']
        Gc_uv   = mp['cardio_control_params.Gc_uv']['value']
        tau_uv  = mp['cardio_control_params.tau_uv']['value']

        hr_err     = recent_MAP - ABP_n
        dDelta_HR_c = (-Delta_HR_c + Gc_hr * hr_err) / tau_hr
        dDelta_R_c  = (-Delta_R_c  + Gc_r  * hr_err) / tau_r
        dDelta_UV_c = (-Delta_UV_c + Gc_uv * hr_err) / tau_uv

        # 2) chemoreflex on RR & Pmus
        PaCO2_n = mp['respiratory_control_params.PaCO2_n']['value']
        Gc_A    = mp['respiratory_control_params.Gc_A']['value']
        tau_c_A = mp['respiratory_control_params.tau_c_A']['value']
        Gc_f    = mp['respiratory_control_params.Gc_f']['value']
        tau_p_f = mp['respiratory_control_params.tau_p_f']['value']

        u_c            = p_a_CO2 - PaCO2_n
        dDelta_Pmus_c = (-Delta_Pmus_c + Gc_A * u_c) / tau_c_A
        dDelta_RR_c   = (-Delta_RR_c   + Gc_f * u_c) / tau_p_f

        #print(recent_MAP, self.HR)

        return dDelta_RR_c, dDelta_Pmus_c, dDelta_HR_c, dDelta_R_c, dDelta_UV_c






    ## Start-stop calls

    def start_simulation(self):
        """Solve the ODEs and compute variables after integration."""
        self.running = True
        last_print_time = self.t  # Track last print time
        last_emit_time = self.t  # Track last emit time

        while self.running:

            
            # Define time span for the current step
            t_span = [self.t, self.t + self.dt]
            t_eval = [self.t + self.dt]  # Evaluate at the end of dt

            ## Starting computations, variables emitted - Check for disease progression
            for processedEvent in self.events:
                outcome = self.processEvent(processedEvent)
                if outcome is False:
                    self.events.remove(processedEvent)  # Remove the event after processing, else keep it


            # 1) integrate one step
            sol = solve_ivp(self.extended_state_space_equations,
                            [self.t, self.t + self.dt],
                            self.current_state,
                            t_eval=[self.t + self.dt],
                            method='RK45')

            # 2) update your state
            self.current_state = sol.y[:, -1]

            # 3) extract the reflex‐related state variables
            p_a_CO2      = self.current_state[17]
            Delta_RR_c   = self.current_state[23]
            Delta_Pmus_c = self.current_state[24]
            Delta_HR_c   = self.current_state[26]
            Delta_R_c    = self.current_state[27]
            Delta_UV_c   = self.current_state[28]

            # 4) now invoke your reflex‐logic
            self.apply_reflexes(p_a_CO2,
                                Delta_RR_c,
                                Delta_Pmus_c,
                                Delta_HR_c,
                                Delta_R_c,
                                Delta_UV_c)
            # Compute variables at the latest time point
            P, F, HR, Sa_O2 = self.compute_variables(sol.t[-1], self.current_state)
            #self.current_heart_rate = self.HR  # Update monitored value
            #self.current_SaO2 = Sa_O2  # Update monitored value
            #self.current_MAP = np.mean(self.P_store[0])  # Update monitored value

            # Update buffer index
            self.P_store[:, self.buffer_index] = P[0]
            self.F_store[:, self.buffer_index] = F
            self.HR_store[self.buffer_index] = HR
            self.buffer_index = (self.buffer_index + 1) % 2000

            self.current_heart_rate = self.HR  # Update monitored value
            self.current_SaO2 = Sa_O2  # Update monitored value
            self.current_MAP = float(np.mean(self.P_buffer))  # Update monitored value
            self.current_RR = self.RR  # Update monitored value

            #print(self.buffer_index)
            #print(self.P_store[0])
            # Print heart rate at intervals
            if self.t - last_print_time >= self.print_interval:
                # print MAP, SAP, DAP
                MAP = np.mean(self.P_store[0])
                SAP = np.max(self.P_store[0])
                DAP = np.min(self.P_store[0])
                RR = self.RR
                Sa_O2 = self.current_SaO2
                print(f"MAP, SAP, DAP, HR, SaO2, RR at time {self.t:.2f}s for patient {self.patient_id}: {MAP} {SAP} {DAP} {HR} {Sa_O2} {RR}")
                #print(f"HR, p_a_O2 at time {self.t:.2f}s for patient {self.patient_id}: {HR} {self.current_state[18]} ")
                last_print_time = self.t
 

            # Emit data to the client every 1 seconds = output_frequency
            if self.t - last_emit_time >= self.output_frequency:
                data = {'time': self.start_timestamp + timedelta(seconds=np.round(self.t)),  ## Round timestep to the second
                        'values':{
                            "heart_rate": np.round(self.current_heart_rate, 2),
                            "SaO2": np.round(self.current_SaO2, 2),
                            "MAP": np.round(np.mean(self.P_store), 2),
                            "SAP": np.round(np.max(self.P_store), 2), ## CHECK LEX
                            "DAP": np.round(np.min(self.P_store), 2), ## CHECK LEX
                            "RR": np.round(self.RR, 2),
                            "etCO2": np.round(self.current_state[17], 2)
                    } }
                
                ## Send data to the client
                #self.redisClient.add_vital_sign(self.patient_id, data)

                ## Save latest data as epoch: Keep n data points for callback purposes
                if len(self.data_epoch) == 0:
                    self.data_epoch = [data]
                else: 
                    self.data_epoch.append(data)
                    self.data_epoch = self.data_epoch[ -self.data_points:] # Keep last N data points

                ## Send the data to the Alarm-module 
                self.alarmModule.evaluate_data(curr_data = data, historic_data = self.data_epoch)
                last_emit_time = self.t



            # Update simulation time and pause
            self.t += self.dt
            if self.sleep:
                time.sleep(self.dt)  # Control simulation speed

    def addProcessedEvent(self, event):
        """Add a processed event to the event queue."""
        if event:
            self.events.append(event)
            print(f"Event added for patient {self.patient_id}: {event['eventType']}")

    def processEvent(self, eventContent):
        outcome = None
        error = None
        """Process an event and update the model parameters."""
        if eventContent["eventType"] == "common": ## routine, no specific things -> Special events re-route to custom functions
            if eventContent["timeCategorical"] == "continuous" or eventContent["timeCategorical"] == "limited":
                ## ongoing event, check delta T 
                outcome = True ## keep the event 
                lastEmissionTime = eventContent["lastEmission"] if eventContent["lastEmission"] else 0
                if self.t - lastEmissionTime <= int(timedelta(**{eventContent["timeUnit"]: eventContent["timeInterval"]}).total_seconds()):
                    ## Not enough time has passed: Break early
                    return
                else:
                    ## Enough time has passed: Update last emission time
                    eventContent["lastEmission"] = self.t
                    eventContent["eventCount"] -= 1
                    if eventContent["eventCount"] == 0:
                        outcome = False ## remove the event after this processing round

            ## continue with the event processing
            for paramChange in eventContent['parameters']:
                paramName = paramChange['name']
                paramValChange = paramChange['value']
                ## paramValChange is in percentages compared to [-100 | 100]
                if paramName in self.master_parameters:
                    ## A] Type = 'relative' -> convert current value to a percentage - add percentage and recalculate
                    if paramChange['type'] == 'relative':
                        ## get the current value of the parameter
                        currValue = self.master_parameters[paramName]['value']
                        ## get the current percentage of the parameter
                        currPercent = self.pathologies.CalcParamPercent(currValue, paramName)

                        if paramChange['action'] == 'decay':
                            ## decay the parameter
                            currPercent += paramValChange ## decay goes both ways: in- and decrease facilitated
                        elif paramChange['action'] == 'set':
                            currPercent = paramValChange
                        else:
                            error = f"Error: Unknown action {paramChange} for parameter {paramName} in patient {self.patient_id}"

                        ## get the actual value for the parameter in non-percentualness
                        splineFunction = self.pathologies.splineFunctions[paramName]
                        newValue = splineFunction(currPercent)
                        ## Update the parameter in the model
                        self.master_parameters[paramName]['value'] = newValue
                        print(f" === EVENT ACTION: param {paramName} updated to {self.master_parameters[paramName]['value']} ===")

                    ## B] Type = 'absolute' -> set directly, no inference of spline
                    elif paramChange['type'] == 'absolute':
                        if paramChange['action'] == 'decay':
                            ## decay the parameter
                            self.master_parameters[paramName]['value'] += paramChange['value']
                        elif paramChange['action'] == 'set':
                            ## set parameter direcly
                            self.master_parameters[paramName]['value'] = paramChange['value']
                        else:
                            error = f"Error: Unknown action {paramChange} for parameter {paramName} in patient {self.patient_id}"
                    else:
                        error = f"Error: Unknown action {paramChange} for parameter {paramName} in patient {self.patient_id}"
                        return False
            
            ## recalculate cardiac elastances, resistance and uVolme
            self.compute_cardiac_parameters()
        
        elif eventContent["eventType"] == "special": 
            ## Custom-made function for special events
            outcome = False
            print(f"Special events not yet incorporated, removing from list")
            ## to be filled later!

        if error:
             return error 
        else: 
            return outcome ## True = event can stay for next iteration, False = event to be deleted
            
    def stop_simulation(self):
        """Stop the simulation loop."""
        if self.running:
            self.running = False
            print(f"Simulation stopping for patient {self.patient_id}...")
        else:
            print(f"Simulation for patient {self.patient_id} is not running.")
        return {"status": f"Simulation stopped for patient {self.patient_id}"}
    

    def update_param(self, patient_id, param, value):
        """ Update non-calculated parameters due to outbound communications"""
        param = str(param)
        loc = None ## Initialize

        ## Mapping of parameter names to their indices in the state vector
        if param.split('|')[0] == "resistance" or param.split('|')[0] == "uvolume" or param.split('|')[0] == "elastance":
            loc = int(param.split('|')[1]) -1 ## correct for 0'th position
            param = str(param.split('|')[0])
        possibilities = ["params", "initial_conditions","cardio_parameters", "bloodflows", "respi_constants", "gas_exchange_params", "respiratory_control_params", "cardio_control_params", "misc_constants"]
        for key in possibilities:
            attr = getattr(self, key)
            if param in attr.keys():

                if loc:
                    for i in range(len(attr[param])):
                        if i == loc:
                            attr[param][loc] = float(value)
                else:
                    attr[param] = float(value)
                return {"status": f"Parameter {param} updated to {value} for patient {patient_id}"}
                #break
        else: 
            return(f"Error: Parameter {param} not found in the parameter dictionaries for patient {patient_id}")