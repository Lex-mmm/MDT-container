# digital_twin_model.py

import numpy as np
from scipy.integrate import solve_ivp
from Inference.inference_calc import Inference
import time
import json
from datetime import datetime
import threading


class DigitalTwinModel:
    def __init__(self, patient_id, param_file="sepsis.json", 
                 data_callback=None, 
                 alarm_callback=None, 
                 sleep=True,
                 time_step=0.01
                 ): ## modify for brute-force jobs
        self.patient_id = patient_id
        self.running = False
        self.t = 0
        self.dt = time_step  # Time step
        self.print_interval = 5  # Interval for printing heart rate
        self.data_callback = data_callback  # Callback function to emit data
        self.output_frequency = 1  # Output frequency for data callback -> 1 Hz

        ## Integrate starting levels of measurement uncertainty
        self.inference = Inference()

        self.sleep = sleep  # Sleep between iterations, boolean

        ## set datapoint amount in local memory
        self.data_points = 120 ## 2 minutes 

        self.data_epoch = {}

        # Load parameters from JSON
        self._load_parameters(param_file)

        # Initialize model parameters before state
        self.initialize_model_parameters()

        # Initialize state variables
        self.current_state = self.initialize_state()

        self.current_heart_rate = 0  # Initialize monitored value

    ## Add definition for API setting new parameter-settings dependent on .json format
    def update_parameters(self, param_file_new):
        self._load_parameters(param_file_new)


    def _load_parameters(self, param_file):
        """
        Load parameters from a JSON configuration file.
        Raise exceptions with detailed messages for missing or invalid files.
        """
        try:
            print(f"Loading parameter file: {param_file}")
            with open(param_file, "r") as file:
                config = json.load(file)
            self.params = config.get("params", {})
            self.initial_conditions = config.get("initial_conditions", {})
            self.bloodflows = config.get("bloodflows", {})
            self.respi_constants = config.get("respi_constants", {})
            self.cardio_parameters = config.get("cardio_parameters", {})
            self.gas_exchange_params = config.get("gas_exchange_params", {})
            self.derived_gas_exchange_params = config.get("derived_gas_exchange_params", {})
            self.respiratory_control_params = config.get("respiratory_control_params", {})
            self.cardio_control_params = config.get("cardio_control_params", {})
            self.misc_constants = config.get("misc_constants", {})
            print(f"Successfully loaded parameters for patient {self.patient_id}.")
        except FileNotFoundError:
            raise FileNotFoundError(f"Parameter file '{param_file}' not found. Please verify the file path.")
        except json.JSONDecodeError as e:
            raise ValueError(f"Parameter file '{param_file}' is not a valid JSON file: {e}")

        # Process parameters that involve expressions
        self.process_parameter_expressions()

    def process_parameter_expressions(self):
        """Evaluate expressions in parameters and update the dictionaries."""
        # Combine all parameters into a single context for evaluation
        context = {}
        context.update({"np": np, "exp": np.exp, "min": min, "max": max})

        # Collect all parameter dictionaries into one for evaluation
        all_parameters = {}
        parameter_dicts = [
            self.params,
            self.initial_conditions,
            self.bloodflows,
            self.respi_constants,
            self.cardio_parameters,
            self.gas_exchange_params,
            self.derived_gas_exchange_params,
            self.respiratory_control_params,
            self.cardio_control_params,
            self.misc_constants,
        ]

        for param_dict in parameter_dicts:
            all_parameters.update(param_dict)

        # Keep track of unresolved parameters
        unresolved = set(all_parameters.keys())
        resolved = set()
        max_iterations = 10  # Prevent infinite loops

        for iteration in range(max_iterations):
            progress_made = False
            for key in list(unresolved):
                value = all_parameters[key]
                if isinstance(value, str):
                    try:
                        # Attempt to evaluate the expression
                        evaluated_value = eval(value, context)
                        all_parameters[key] = evaluated_value
                        context[key] = evaluated_value
                        unresolved.remove(key)
                        resolved.add(key)
                        progress_made = True
                    except Exception:
                        # Cannot evaluate yet; dependencies might not be resolved
                        continue
                else:
                    # Value is already a number or list
                    context[key] = value
                    unresolved.remove(key)
                    resolved.add(key)
                    progress_made = True

            if not progress_made:
                break  # No further progress can be made

        if unresolved:
            print(f"Could not resolve the following parameters after {iteration + 1} iterations:")
            for key in unresolved:
                print(f" - {key}")
            raise ValueError("Parameter evaluation failed due to unresolved dependencies.")
        else:
            print("All parameters resolved successfully.")

        # Update the original dictionaries with evaluated parameters
        for param_dict in parameter_dicts:
            for key in param_dict.keys():
                param_dict[key] = all_parameters[key]

        # Ensure 't_eval' is defined
        if 't_eval' not in self.misc_constants:
            self.misc_constants['t_eval'] = np.arange(
                self.misc_constants['tmin'],
                self.misc_constants['tmax'] + self.misc_constants['T'],
                self.misc_constants['T']
            )

    def initialize_model_parameters(self):
        """Initialize model parameters like elastance, resistance, uvolume, etc."""
        # Initialize elastance, resistance, uvolume from parameters
        elastance_list = self.cardio_parameters['elastance']
        resistance_list = self.cardio_parameters['resistance']
        uvolume_list = self.cardio_parameters['uvolume']

        # Convert lists to numpy arrays
        self.elastance = np.array(elastance_list, dtype=np.float64)
        self.resistance = np.array(resistance_list, dtype=np.float64)
        self.uvolume = np.array(uvolume_list, dtype=np.float64)

        # Lung model equations
        # Mechanical system parameters
        self.C_cw = 0.2445  # l/cmH2O

        self.A_mechanical = np.array([
            [-1 / (self.respi_constants['C_l'] * self.respi_constants['R_ml']) - 1 / (self.respi_constants['R_lt'] * self.respi_constants['C_l']), 1 / (self.respi_constants['R_lt'] * self.respi_constants['C_l']), 0, 0, 0],
            [1 / (self.respi_constants['R_lt'] * self.respi_constants['C_tr']), -1 / (self.respi_constants['C_tr'] * self.respi_constants['R_lt']) - 1 / (self.respi_constants['R_tb'] * self.respi_constants['C_tr']), 1 / (self.respi_constants['R_tb'] * self.respi_constants['C_tr']), 0, 0],
            [0, 1 / (self.respi_constants['R_tb'] * self.respi_constants['C_b']), -1 / (self.respi_constants['C_b'] * self.respi_constants['R_tb']) - 1 / (self.respi_constants['R_bA'] * self.respi_constants['C_b']), 1 / (self.respi_constants['R_bA'] * self.respi_constants['C_b']), 0],
            [0, 0, 1 / (self.respi_constants['R_bA'] * self.respi_constants['C_A']), -1 / (self.respi_constants['C_A'] * self.respi_constants['R_bA']), 0],
            [1 / (self.respi_constants['R_lt'] * self.C_cw), -1 / (self.C_cw * self.respi_constants['R_lt']), 0, 0, 0]
        ])

        self.B_mechanical = np.array([
            [1 / (self.respi_constants['R_ml'] * self.respi_constants['C_l']), 0, 0],
            [0, 1, 0],
            [0, 1, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])

        # Initialize heart period parameters
        self.HR = self.misc_constants.get('HR', 75)
        self.update_heart_period(self.HR)

        # Initialize buffer for sliding window
        self.dt = self.misc_constants.get('T', 0.01)
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
        TBV = self.misc_constants.get('TBV', 5000)
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
        state[17] = self.initial_conditions.get('p_a_CO2', 40)
        state[18] = self.initial_conditions.get('p_a_O2', 95)
        state[19] = self.initial_conditions.get('c_Stis_CO2', 0.5)
        state[20] = self.initial_conditions.get('c_Scap_CO2', 0.5)
        state[21] = self.initial_conditions.get('c_Stis_O2', 0.2)
        state[22] = self.initial_conditions.get('c_Scap_O2', 0.2)
        state[23] = self.respiratory_control_params.get('Delta_RR_c', 0)
        state[24] = self.respiratory_control_params.get('Delta_Pmus_c', 0)
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
        """Calculate heart elastances."""
        HP = self.HP
        Tas = self.Tas
        Tav = self.Tav
        Tvs = self.Tvs
        T = self.misc_constants['T']
        ncc = (t % HP) / T

        if ncc <= round(Tas / T):
            aaf = np.sin(np.pi * ncc / (Tas / T))
        else:
            aaf = 0

        ela = self.elastance[0, 8] + (self.elastance[1, 8] - self.elastance[0, 8]) * aaf
        era = self.elastance[0, 4] + (self.elastance[1, 4] - self.elastance[0, 4]) * aaf

        if ncc <= round((Tas + Tav) / T):
            vaf = 0
        elif ncc <= round((Tas + Tav + Tvs) / T):
            vaf = np.sin(np.pi * (ncc - (Tas + Tav) / T) / (Tvs / T))
        else:
            vaf = 0

        elv = self.elastance[0, 9] + (self.elastance[1, 9] - self.elastance[0, 9]) * vaf
        erv = self.elastance[0, 5] + (self.elastance[1, 5] - self.elastance[0, 5]) * vaf

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
        UV_c = self.cardio_control_params['UV_n'] + y[28]
        P[2] = self.elastance[0, 2] * (V[2] - self.uvolume[2] * UV_c)
        P[3] = self.elastance[0, 3] * (V[3] - self.uvolume[3] * UV_c) + Pmus
        P[4] = era * (V[4] - self.uvolume[4]) + Pmus
        P[5] = erv * (V[5] - self.uvolume[5]) + Pmus
        P[6] = self.elastance[0, 6] * (V[6] - self.uvolume[6]) + Pmus
        P[7] = self.elastance[0, 7] * (V[7] - self.uvolume[7]) + Pmus
        P[8] = ela * (V[8] - self.uvolume[8]) + Pmus
        P[9] = elv * (V[9] - self.uvolume[9]) + Pmus

        # Calculate the flows for cardiovascular model
        R_c = self.cardio_control_params['R_n'] - y[27]  # Delta_R_c is y[27]
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
        HR = self.cardio_control_params['HR_n'] - y[26]  # Delta_HR_c is y[26]

        # Compute SaO2
        p_a_O2 = y[18]
        CaO2 = (self.params['K_O2'] * np.power((1 - np.exp(-self.params['k_O2'] * min(p_a_O2, 700))), 2)) * 100
        Sa_O2 = np.round(((CaO2 - p_a_O2 * 0.003 / 100) / (self.misc_constants['Hgb'] * 1.34)) * 100)

        return P, F, HR, Sa_O2

    def ventilator_pressure(self, t):
        """Ventilator pressure as a function of time (simple square wave for demonstration)."""
        RR = self.respiratory_control_params['RR_0']  # Respiratory Rate (breaths per minute)
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
        """Define the combined system of differential equations."""
        # Use self.HR_store, self.buffer_index, self.window_size, self.HR
        # Split x into cardiovascular and lung model components
        V = x[:10]  # Volumes for cardiovascular model
        mechanical_states = x[10:15]  # Mechanical states
        FD_O2, FD_CO2, p_a_CO2, p_a_O2 = x[15], x[16], x[17], x[18]
        c_Stis_CO2, c_Scap_CO2, c_Stis_O2, c_Scap_O2 = x[19:23]
        Delta_RR_c = x[23]
        Delta_Pmus_c = x[24]
        Pmus = x[25]
        Delta_HR_c = x[26]
        Delta_R_c = x[27]
        Delta_UV_c = x[28]

        # Inputs for cardiovascular model
        ela, elv, era, erv = self.get_inputs(t)

        HR = self.cardio_control_params['HR_n'] - Delta_HR_c
        R_c = self.cardio_control_params['R_n'] - Delta_R_c
        UV_c = self.cardio_control_params['UV_n'] + Delta_UV_c

        # Update heart period parameters
        self.HR = HR
        self.update_heart_period(HR)

        if self.misc_constants['MV'] == 0:
            P_ao = 0
            self.RR = self.respiratory_control_params['RR_0'] + Delta_RR_c
            Pmus_min = self.respiratory_control_params['Pmus_0'] + Delta_Pmus_c
            driver = self.input_function(t, self.RR, Pmus_min)
            Pmus_dt = driver[1]
            FI_O2 = self.gas_exchange_params['FI_O2']
            FI_CO2 = self.gas_exchange_params['FI_CO2']
        else:
            P_ao = self.ventilator_pressure(t)
            driver = np.array([P_ao, 0])
            RR = 12
            FI_O2 = self.gas_exchange_params['FI_O2']
            FI_CO2 = self.gas_exchange_params['FI_CO2']
            Pmus_dt = 0

        #print(RR, Pmus_min, FI_O2, FI_CO2)

        # Calculate the pressures for cardiovascular model
        P = np.zeros(10)
        P[0] = self.elastance[0, 0] * (V[0] - self.uvolume[0]) + mechanical_states[4]
        P[1] = self.elastance[0, 1] * (V[1] - self.uvolume[1])
        P[2] = self.elastance[0, 2] * (V[2] - self.uvolume[2] * UV_c)
        P[3] = self.elastance[0, 3] * (V[3] - self.uvolume[3] * UV_c) + mechanical_states[4]
        P[4] = era * (V[4] - self.uvolume[4]) + mechanical_states[4]
        P[5] = erv * (V[5] - self.uvolume[5]) + mechanical_states[4]
        P[6] = self.elastance[0, 6] * (V[6] - self.uvolume[6]) + mechanical_states[4]
        P[7] = self.elastance[0, 7] * (V[7] - self.uvolume[7]) + mechanical_states[4]
        P[8] = ela * (V[8] - self.uvolume[8]) + mechanical_states[4]
        P[9] = elv * (V[9] - self.uvolume[9]) + mechanical_states[4]

        # Calculate the flows for cardiovascular model
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

        # Store values in the circular buffer
        #self.P_store[:, self.buffer_index] = P
        #self.F_store[:, self.buffer_index] = F
        #self.HR_store[self.buffer_index] = HR

        # Update buffer index
        #self.buffer_index = (self.buffer_index + 1) % self.window_size

        # Calculate the derivatives of volumes for cardiovascular model
        dVdt = np.zeros(10)
        dVdt[0] = F[9] - F[0]
        dVdt[1] = F[0] - F[1]
        dVdt[2] = F[1] - F[2]
        dVdt[3] = F[2] - F[3]
        dVdt[4] = F[3] - F[4]
        dVdt[5] = F[4] - F[5]
        dVdt[6] = F[5] - F[6]
        dVdt[7] = F[6] - F[7]
        dVdt[8] = F[7] - F[8]
        dVdt[9] = F[8] - F[9]

        # Lung cardio to respi model by exchanging the Cardiac output with the lung model
        # Calculate CO as area under the curve of F[0]
        #cycle_length = 60 / HR
        #cycle_points = int(cycle_length / self.dt)
        idx = self.buffer_index  # Current buffer index
        #if idx >= cycle_points:
        #    F0_cycle = self.F_store[0, idx - cycle_points:idx]
        #    CO = np.trapz(F0_cycle, dx=self.dt)   # in ml/sec
        #else:
        #    CO = self.bloodflows['CO']  # Initial values

        CO = self.bloodflows['CO']
        q_p = CO
        sh = 0.02  # Shunt fraction
        q_Bv = 0.2 * CO
        q_S = 0.8 * CO

        # Compute mechanical derivatives
        Pmus_dt = driver[1]
        Ppl_dt = (mechanical_states[0] / (self.respi_constants['R_lt'] * self.C_cw)) - (mechanical_states[1] / (self.C_cw * self.respi_constants['R_lt'])) + Pmus_dt
        dxdt_mechanical = np.dot(self.A_mechanical, mechanical_states) + np.dot(self.B_mechanical, [P_ao, Ppl_dt, Pmus_dt])

        # Compute the ventilation rates based on pressures
        Vdot_l = (P_ao - mechanical_states[0]) / self.respi_constants['R_ml']
        Vdot_A = (mechanical_states[2] - mechanical_states[3]) / self.respi_constants['R_bA']

        p_D_CO2 = FD_CO2 * 713
        p_D_O2 = FD_O2 * 713
        FA_CO2 = p_a_CO2 / 713
        FA_O2 = p_a_O2 / 713

        c_a_CO2 = self.params['K_CO2'] * p_a_CO2 + self.params['k_CO2']  # Arterial CO2 concentration

        # Safeguard to prevent overflow
        if p_a_O2 > 700:
            p_a_O2 = 700

        c_a_O2 = self.params['K_O2'] * np.power((1 - np.exp(-self.params['k_O2'] * p_a_O2)), 2)

        c_v_CO2 = c_Scap_CO2
        c_v_O2 = c_Scap_O2

        # Determine inspiration or expiration
        if (self.misc_constants['MV'] == 1 and P_ao > 6 * 0.735) or (self.misc_constants['MV'] == 0 and mechanical_states[0] < 0):
            dFD_O2_dt = Vdot_l * 1000 * (FI_O2 - FD_O2) / (self.gas_exchange_params['V_D'] * 1000)
            dFD_CO2_dt = Vdot_l * 1000 * (FI_CO2 - FD_CO2) / (self.gas_exchange_params['V_D'] * 1000)

            dp_a_CO2 = (863 * q_p * (1 - sh) * (c_v_CO2 - c_a_CO2) + Vdot_A * 1000 * (p_D_CO2 - p_a_CO2)) / (self.gas_exchange_params['V_A'] * 1000)
            dp_a_O2 = (863 * q_p * (1 - sh) * (c_v_O2 - c_a_O2) + Vdot_A * 1000 * (p_D_O2 - p_a_O2)) / (self.gas_exchange_params['V_A'] * 1000)
        else:
            # Expiration
            dFD_O2_dt = Vdot_A * 1000 * (FD_O2 - FA_O2) / (self.gas_exchange_params['V_D'] * 1000)
            dFD_CO2_dt = Vdot_A * 1000 * (FD_CO2 - FA_CO2) / (self.gas_exchange_params['V_D'] * 1000)

            dp_a_CO2 = 863 * q_p * (1 - self.bloodflows['sh']) * (c_v_CO2 - c_a_CO2) / (self.gas_exchange_params['V_A'] * 1000)
            dp_a_O2 = 863 * q_p * (1 - self.bloodflows['sh']) * (c_v_O2 - c_a_O2) / (self.gas_exchange_params['V_A'] * 1000)

        # The systemic tissue compartment
        dc_Stis_CO2 = (self.params['M_S_CO2'] - self.gas_exchange_params['D_S_CO2'] * (c_Stis_CO2 - c_Scap_CO2)) / self.params['V_Stis_CO2']
        dc_Scap_CO2 = (q_S * (c_a_CO2 - c_Scap_CO2) + self.gas_exchange_params['D_S_CO2'] * (c_Stis_CO2 - c_Scap_CO2)) / self.params['V_Scap_CO2']
        dc_Stis_O2 = (self.params['M_S_O2'] - self.gas_exchange_params['D_S_O2'] * (c_Stis_O2 - c_Scap_O2)) / self.params['V_Stis_O2']
        dc_Scap_O2 = (q_S * (c_a_O2 - c_Scap_O2) + self.gas_exchange_params['D_S_O2'] * (c_Stis_O2 - c_Scap_O2)) / self.params['V_Scap_O2']

        # Central control
        u_c = p_a_CO2 - self.respiratory_control_params['PaCO2_n']
        hr_c = P[0] - self.cardio_control_params['ABP_n']
        dDelta_Pmus_c = (-Delta_Pmus_c + self.respiratory_control_params['Gc_A'] * u_c) / self.respiratory_control_params['tau_c_A']
        dDelta_RR_c = (-Delta_RR_c + self.respiratory_control_params['Gc_f'] * u_c) / self.respiratory_control_params['tau_p_f']

        dDelta_HR_c = (-Delta_HR_c + self.cardio_control_params['Gc_hr'] * hr_c) / self.cardio_control_params['tau_hr']  # Heart rate
        dDelta_R_c = (-Delta_R_c + self.cardio_control_params['Gc_r'] * hr_c) / self.cardio_control_params['tau_r']  # Resistance
        dDelta_UV_c = (-Delta_UV_c + self.cardio_control_params['Gc_uv'] * hr_c) / self.cardio_control_params['tau_uv']  # Unstressed volume

        # Combine all derivatives into a single array
        dxdt = np.concatenate([
            dVdt,
            dxdt_mechanical,
            [
                dFD_O2_dt,
                dFD_CO2_dt,
                dp_a_CO2,
                dp_a_O2,  
                dc_Stis_CO2,
                dc_Scap_CO2,
                dc_Stis_O2,
                dc_Scap_O2,
                dDelta_RR_c,
                dDelta_Pmus_c,
                Pmus_dt,
                dDelta_HR_c,
                dDelta_R_c,
                dDelta_UV_c
            ]
        ])
        return dxdt
    


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
            sol = solve_ivp(self.extended_state_space_equations, t_span, self.current_state, t_eval=t_eval, method='RK45')
            self.current_state = sol.y[:, -1]  # Update state to the latest solution

            # Compute variables at the latest time point
            P, F, HR, Sa_O2 = self.compute_variables(sol.t[-1], self.current_state)
            self.current_heart_rate = self.HR  # Update monitored value
            self.current_SaO2 = Sa_O2  # Update monitored value

            # Update buffer index
            self.P_store[:, self.buffer_index] = P[0]
            self.F_store[:, self.buffer_index] = F
            self.HR_store[self.buffer_index] = HR
            self.buffer_index = (self.buffer_index + 1) % 2000

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
                data = {'time': np.round(self.t),  ## Round timestep to the second
                        'values':{
                            "heart_rate": np.round(self.current_heart_rate, 2),
                            "SaO2": np.round(self.current_SaO2, 2),
                            "MAP": np.round(np.mean(self.P_store), 2),
                            "RR": np.round(self.RR, 2),
                            "etCO2": np.round(self.current_state[17], 2)
                    } }
                ## Do inference checks, add to existing values
                data['values'] = self.inference.check_inference(data['values'] ) 


                ## Save latest data as epoch: Keep n data points for callback purposes
                if len(self.data_epoch) == 0:
                    self.data_epoch = [data]
                else: 
                    self.data_epoch.append(data)
                    self.data_epoch = self.data_epoch[ -self.data_points:] # Keep last N data points

                last_emit_time = self.t

            # Update simulation time and pause
            self.t += self.dt
            if self.sleep:
                time.sleep(self.dt)  # Control simulation speed





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