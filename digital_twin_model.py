import numpy as np
from scipy.signal import butter, filtfilt
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
    def __init__(self, patient_id, param_file="healthyFlat", 
                 data_callback=None, 
                 sleep=True,
                 time_step=0.02):  # time_step argument will be effectively overridden by parameters
        self._initialize_basic_attributes(patient_id, data_callback, sleep)
        self._initialize_timing_attributes()
        self._initialize_modules()
        self._initialize_events_and_data()

        # Load parameters and initialize the model
        self._load_and_process_parameters(param_file)
        self._initialize_simulation_environment()

    def _initialize_basic_attributes(self, patient_id, data_callback, sleep):
        """Initialize basic attributes like patient ID, callback, and sleep mode."""
        self.patient_id = patient_id
        self.running = False
        self.t = 0
        self.data_callback = data_callback
        self.sleep = sleep

    def _initialize_timing_attributes(self):
        """Initialize timing-related attributes."""
        self.P_buffer_sum = 0  # Store the running sum
        self.print_interval = 5  # Interval for printing heart rate
        self.output_frequency = 10  # Output frequency for data callback -> 1 Hz

    def _initialize_modules(self):
        """Initialize external modules like pathologies, therapies, and alarms."""
        self.pathologies = Pathology()  # Initialize pathology events
        self.therapeutic = Therapy()
        self.alarmModule = alarmModule(self.patient_id)  # Initialize alarm module

    def _initialize_events_and_data(self):
        """Initialize event and data-related attributes."""
        self.events = []  # Actionable event introduction
        self.data_points = 120  # 2 minutes
        self.data_epoch = []  # Initialize as a list to avoid AttributeError
        self.start_timestamp = datetime.now()

    def _load_and_process_parameters(self, param_file):
        """Load parameters from a file and process them."""
        self._load_parameters(param_file)
        self.dt = self.master_parameters['misc_constants.T']['value']
        self.window_size = int(30 / self.dt)  # 30 seconds window
        self.P_buffer = deque([0.0] * self.window_size, maxlen=self.window_size)
        
    def _initialize_simulation_environment(self):
        """Initialize the simulation environment."""
        self._cache_baroreflex_parameters()
        self._cache_chemoreceptor_parameters()
        self._compute_all_derived_params()
        self.compute_cardiac_parameters()
        self._cache_ode_parameters()
        self._setup_simulation_environment()
        self.current_state = self.initialize_state()
        self.current_heart_rate = 0  # Initialize monitored value
        self.use_reflex = True  # Set to False if you want to skip reflex calculations
        self.AF = 0  # Atrial Fibrillation flag
        self.Fes_delayed = [2.66] * int(2 / self.dt)
        self.Fev_delayed = [4.66] * int(0.2 / self.dt)
        
        # Chemoreceptor delayed buffers
        self.p_CO2_delayed = [40.0] * int(self._chemo_delay_CO2 / self.dt)
        self.p_O2_delayed = [100.0] * int(self._chemo_delay_O2 / self.dt)
        
        # Simple cycle tracking
        self.last_cycle_start = 0.0  # Track when last cycle started

    def _cache_baroreflex_parameters(self):
        """Cache parameters used frequently in baroreceptor_control for speed."""
        self._baro_tz = self.master_parameters['baroreflex_params.tz']['value']
        self._baro_tp = self.master_parameters['baroreflex_params.tp']['value']
        self._baro_Fas_min = self.master_parameters['baroreflex_params.Fas_min']['value']
        self._baro_Fas_max = self.master_parameters['baroreflex_params.Fas_max']['value']
        self._baro_Ka = self.master_parameters['baroreflex_params.Ka']['value']
        # P_set is derived from ABP_n
        self._baro_P_set = self.master_parameters['cardio_control_params.ABP_n']['value']
        self._baro_Fes_inf = self.master_parameters['baroreflex_params.Fes_inf']['value']
        self._baro_Fes_0 = self.master_parameters['baroreflex_params.Fes_0']['value']
        self._baro_Kes = self.master_parameters['baroreflex_params.Kes']['value']
        self._baro_Fev_0 = self.master_parameters['baroreflex_params.Fev_0']['value']
        self._baro_Fev_inf = self.master_parameters['baroreflex_params.Fev_inf']['value']
        self._baro_Kev = self.master_parameters['baroreflex_params.Kev']['value']
        self._baro_Fas_0 = self.master_parameters['baroreflex_params.Fas_0']['value']

    def _cache_chemoreceptor_parameters(self):
        """Cache parameters used frequently in chemoreceptor_control for speed."""
        # CO2 chemoreceptor parameters - REDUCED GAINS
        self._chemo_tau_CO2 = self.master_parameters.get('chemoreceptor_params.tau_CO2', {'value': 8.0})['value']  # Increased time constant
        self._chemo_G_CO2 = self.master_parameters.get('chemoreceptor_params.G_CO2', {'value': 0.05})['value']  # Reduced from 0.20 to 0.05
        self._chemo_p_CO2_set = self.master_parameters.get('chemoreceptor_params.p_CO2_set', {'value': 40.0})['value']
        self._chemo_delay_CO2 = self.master_parameters.get('chemoreceptor_params.delay_CO2', {'value': 2.0})['value']  # Increased delay
        
        # O2 chemoreceptor parameters - REDUCED GAINS  
        self._chemo_tau_O2 = self.master_parameters.get('chemoreceptor_params.tau_O2', {'value': 6.0})['value']  # Increased time constant
        self._chemo_G_O2 = self.master_parameters.get('chemoreceptor_params.G_O2', {'value': -0.3})['value']  # Reduced from -1.5 to -0.3
        self._chemo_p_O2_threshold = self.master_parameters.get('chemoreceptor_params.p_O2_threshold', {'value': 105.0})['value']
        self._chemo_delay_O2 = self.master_parameters.get('chemoreceptor_params.delay_O2', {'value': 1.0})['value']  # Increased delay

    def _cache_ode_parameters(self):
        """Cache parameters used frequently in extended_state_space_equations for speed."""
        self._ode_HR_n = self.master_parameters['cardio_control_params.HR_n']['value']
        self._ode_HR_n_max = self.master_parameters['cardio_control_params.HR_n']['max']
        self._ode_HR_n_min = self.master_parameters['cardio_control_params.HR_n']['min']
        self._ode_R_n = self.master_parameters['cardio_control_params.R_n']['value']
        self._ode_R_n_max = self.master_parameters['cardio_control_params.R_n']['max']
        self._ode_R_n_min = self.master_parameters['cardio_control_params.R_n']['min']
        self._ode_UV_n = self.master_parameters['cardio_control_params.UV_n']['value']
        self._ode_RR0 = self.master_parameters['respiratory_control_params.RR_0']['value']
        self._ode_RR0_max = self.master_parameters['respiratory_control_params.RR_0']['max']
        self._ode_RR0_min = self.master_parameters['respiratory_control_params.RR_0']['min']
        self._ode_FI_O2 = self.master_parameters['gas_exchange_params.FI_O2']['value']
        self._ode_FI_CO2 = self.master_parameters['gas_exchange_params.FI_CO2']['value']
        self._ode_CO_nom = self.master_parameters['bloodflows.CO']['value']
        self._ode_shunt = self.master_parameters['bloodflows.sh']['value']
        self._ode_MV_mode = self.master_parameters['misc_constants.MV']['value']
        self._ode_Pmus_0 = self.master_parameters['respiratory_control_params.Pmus_0']['value']
        self._ode_R_lt = self.master_parameters['respi_constants.R_lt']['value']
        self._ode_R_bA = self.master_parameters['respi_constants.R_bA']['value']
        self._ode_K_CO2 = self.master_parameters['params.K_CO2']['value']
        self._ode_k_CO2 = self.master_parameters['params.k_CO2']['value']
        self._ode_K_O2 = self.master_parameters['params.K_O2']['value']
        self._ode_k_O2 = self.master_parameters['params.k_O2']['value']
        self._ode_V_D = self.master_parameters['gas_exchange_params.V_D']['value']
        self._ode_V_A = self.master_parameters['gas_exchange_params.V_A']['value']
        self._ode_M_S_CO2 = self.master_parameters['params.M_S_CO2']['value']
        self._ode_D_S_CO2 = self.master_parameters['gas_exchange_params.D_S_CO2']['value']
        self._ode_V_Stis_CO2 = self.master_parameters['params.V_Stis_CO2']['value']
        self._ode_V_Scap_CO2 = self.master_parameters['params.V_Scap_CO2']['value']
        self._ode_M_S_O2 = self.master_parameters['params.M_S_O2']['value']
        self._ode_D_S_O2 = self.master_parameters['gas_exchange_params.D_S_O2']['value']
        self._ode_V_Stis_O2 = self.master_parameters['params.V_Stis_O2']['value']
        self._ode_V_Scap_O2 = self.master_parameters['params.V_Scap_O2']['value']


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
            M_S_CO2    = self.master_parameters['params.M_S_CO2']['value']
            D_T_CO2    = self.master_parameters['params.D_T_CO2']['value']
            self.master_parameters['initial_conditions.c_Scap_CO2'] = {
                'value': c_Stis_CO2 - M_S_CO2 / D_T_CO2
            }

        if 'initial_conditions.c_Scap_O2' not in self.master_parameters:
            c_Stis_O2 = self.master_parameters['initial_conditions.c_Stis_O2']['value']
            M_S_O2    = self.master_parameters['params.M_S_O2']['value']
            D_T_O2    = self.master_parameters['params.D_T_O2']['value']
            self.master_parameters['initial_conditions.c_Scap_O2'] = {
                'value': c_Stis_O2 + M_S_O2 / D_T_O2
            }

        if 'initial_conditions.c_Bcap_CO2' not in self.master_parameters:
            c_Btis_CO2 = self.master_parameters['initial_conditions.c_Btis_CO2']['value']
            M_B_CO2    = self.master_parameters['params.M_B_CO2']['value']
            D_T_CO2    = self.master_parameters['params.D_T_CO2']['value']
            self.master_parameters['initial_conditions.c_Bcap_CO2'] = {
                'value': c_Btis_CO2 - M_B_CO2 / D_T_CO2
            }

        if 'initial_conditions.c_Bcap_O2' not in self.master_parameters:
            c_Btis_O2 = self.master_parameters['initial_conditions.c_Btis_O2']['value']
            M_B_O2    = self.master_parameters['params.M_B_O2']['value']
            D_T_O2    = self.master_parameters['params.D_T_O2']['value']
            self.master_parameters['initial_conditions.c_Bcap_O2'] = {
                'value': c_Btis_O2 + M_B_O2 / D_T_O2
            }

            # --- 3.b populate gas_exchange_params entries so extended_state_space_equations can find them ---
        self.master_parameters['gas_exchange_params.D_S_CO2'] = {'value': self.master_parameters['params.D_T_CO2']['value']}
        self.master_parameters['gas_exchange_params.D_B_CO2'] = {'value': self.master_parameters['params.D_T_CO2']['value']}
        self.master_parameters['gas_exchange_params.D_S_O2']  = {'value': self.master_parameters['params.D_T_O2']['value']}
        self.master_parameters['gas_exchange_params.D_B_O2'] = {'value': self.master_parameters['params.D_T_O2']['value']}

    def compute_cardiac_parameters(self):
        # Initialize elastance, resistance, uvolume from parameters

        cardio_elastance_min = [key for key in self.master_parameters if 'cardio' in key and 'min' in key and 'elastance' in key]
        cardio_elastance_max = [key for key in self.master_parameters if 'cardio' in key and 'max' in key and 'elastance' in key]

        cardio_resistance = [key for key in self.master_parameters if 'cardio' in key and 'resistance' in key]
        cardio_uvolume = [key for key in self.master_parameters if 'cardio' in key and 'uvolume' in key]

        self.elastance = np.array([ [self.master_parameters[key]['value'] for key in cardio_elastance_min], 
                             [self.master_parameters[key]['value'] for key in cardio_elastance_max]])
        
        self.resistance = np.array([self.master_parameters[key]['value'] for key in cardio_resistance])
        
        self.uvolume = np.array([self.master_parameters[key]['value'] for key in cardio_uvolume])


    def _setup_simulation_environment(self):
        """Initialize model structures like mechanical matrices, initial HR/HP, and data stores."""
        self._initialize_mechanical_system()
        self._initialize_heart_period()
        self._initialize_data_buffers()

    def _initialize_mechanical_system(self):
        """Set up the mechanical system matrices."""
        self.C_cw = 0.2445  # l/cmH2O

        # Define A_mechanical matrix
        self.A_mechanical = np.array([
            [-1 / (self.master_parameters['respi_constants.C_l']['value'] * (1.021 * 1.5))  # R_ml constant
             - 1 / (self.master_parameters['respi_constants.R_lt']['value'] * self.master_parameters['respi_constants.C_l']['value']),
             1 / (self.master_parameters['respi_constants.R_lt']['value'] * self.master_parameters['respi_constants.C_l']['value']), 0, 0, 0],
            [1 / (self.master_parameters['respi_constants.R_lt']['value'] * self.master_parameters['respi_constants.C_tr']['value']),
             -1 / (self.master_parameters['respi_constants.C_tr']['value'] * self.master_parameters['respi_constants.R_lt']['value']) - 1 / (self.master_parameters['respi_constants.R_tb']['value'] * self.master_parameters['respi_constants.C_tr']['value']),
             1 / (self.master_parameters['respi_constants.R_tb']['value'] * self.master_parameters['respi_constants.C_tr']['value']), 0, 0],
            [0, 1 / (self.master_parameters['respi_constants.R_tb']['value'] * self.master_parameters['respi_constants.C_b']['value']),
             -1 / (self.master_parameters['respi_constants.C_b']['value'] * self.master_parameters['respi_constants.R_tb']['value']) - 1 / (self.master_parameters['respi_constants.R_bA']['value'] * self.master_parameters['respi_constants.C_b']['value']),
             1 / (self.master_parameters['respi_constants.R_bA']['value'] * self.master_parameters['respi_constants.C_b']['value']), 0],
            [0, 0, 1 / (self.master_parameters['respi_constants.R_bA']['value'] * self.master_parameters['respi_constants.C_A']['value']),
             -1 / (self.master_parameters['respi_constants.C_A']['value'] * self.master_parameters['respi_constants.R_bA']['value']), 0],
            [1 / (self.master_parameters['respi_constants.R_lt']['value'] * self.C_cw),
             -1 / (self.C_cw * self.master_parameters['respi_constants.R_lt']['value']), 0, 0, 0]
        ])

        # Define B_mechanical matrix
        self.B_mechanical = np.array([
            [1 / ((1.021 * 1.5) * self.master_parameters['respi_constants.C_l']['value']), 0, 0],
            [0, 1, 0],
            [0, 1, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])

    def _initialize_heart_period(self):
        """Initialize heart period parameters."""
        self.HR = self.master_parameters['misc_constants.HR']['value']
        self.update_heart_period(self.HR)

    def _initialize_data_buffers(self):
        """Initialize data buffers for sliding windows."""
        self.P_store = deque([0.0] * self.window_size, maxlen=self.window_size)
        self.HR_store = deque([self.HR] * self.window_size, maxlen=self.window_size)
        window_b = int(5 / self.dt)

        self.avg_buffers = {
            "HR": deque(maxlen=window_b),
            "SaO2": deque(maxlen=window_b),
            "RR": deque(maxlen=window_b),
            "MAP": deque(maxlen=window_b),
            "etCO2": deque(maxlen=window_b)
        }



    def initialize_state(self):
        """Set initial state variables based on loaded parameters."""
        # Initialize state variables for the combined model
        state = np.zeros(35)  # Increased from 33 to accommodate new chemoreceptor states
        # Initialize blood volumes based on unstressed volumes
        TBV = self.master_parameters['misc_constants.TBV']['value']
        state[:10] = TBV * (self.uvolume / np.sum(self.uvolume)) 
        # Initialize mechanical states (indices 10 to 14)
        state[10:15] = np.zeros(5)  # Adjust initial values if necessary
        # Initialize other states as before
        state[15] = self.master_parameters['initial_conditions.FD_O2']['value'] / self.master_parameters['initial_conditions.conFrac']['value']  # FD_O2 
        state[16] = self.master_parameters['initial_conditions.FD_CO2']['value'] / self.master_parameters['initial_conditions.conFrac']['value']    # FD_CO2
        state[17] = self.master_parameters['initial_conditions.p_a_CO2']['value']
        state[18] = self.master_parameters['initial_conditions.p_a_O2']['value']
        state[19] = self.master_parameters['initial_conditions.c_Stis_CO2']['value']

        # Use pre-calculated derived initial concentrations for c_Scap_CO2 and c_Scap_O2
        state[20] = self.master_parameters['initial_conditions.c_Scap_CO2']['value']

        state[21] = self.master_parameters['initial_conditions.c_Stis_O2']['value']
        
        state[22] = self.master_parameters['initial_conditions.c_Scap_O2']['value']
        state[23] = self.master_parameters['respiratory_control_params.Delta_RR_c']['value']
        state[24] = self.master_parameters['respiratory_control_params.Delta_Pmus_c']['value']
        state[25] = self.master_parameters['initial_conditions.Pmus']['value']  # Pmus
        state[26] = 0   # Delta_HR_c
        state[27] = 0   # Delta_R_c
        state[28] = 0   # Delta_UV_c
        state[29] = self.master_parameters['initial_conditions.Pset']['value']    # Pbaro (initial arterial pressure)
        state[30] = 0     # dHRv
        state[31] = 0     # dHRh
        state[32] = 0     # dRs
        state[33] = 0     # dRR_chemo (chemoreceptor respiratory rate control)
        state[34] = 0     # dPmus_chemo (chemoreceptor Pmus control)
        return state

    def update_heart_period(self, HR):
        """Calculate heart period parameters."""
        self.HR = HR
        self.HP = 60 / HR
        self.Tas = 0.03 + 0.09 * self.HP
        self.Tav = 0.01
        self.Tvs = 0.16 + 0.2 * self.HP

    def is_cycle_complete(self, t):
        """Check if current cardiac cycle is complete."""
        return (t - self.last_cycle_start) >= self.HP

    def calculate_elastances(self, t):
        """Calculate heart elastances based on the current phase in the heart cycle."""
        # Check if cycle is complete and update timing if needed
        if self.is_cycle_complete(t):
            self.last_cycle_start = t
        
        # Time since start of current cycle
        cycle_time = t - self.last_cycle_start
        
        # Heart period and timing parameters
        HP = self.HP           
        Tas = self.Tas         
        Tav = self.Tav         
        Tvs = self.Tvs         
        
        # Normalize current time within the heart period to a fraction (0 to 1)
        phase_fraction = cycle_time / HP

        # Calculate activation function for contraction (aaf) during the early phase
        if phase_fraction <= Tas / HP:
            # Scale phase so it goes from 0 to pi over the activation time
            aaf = np.sin(np.pi * phase_fraction / (Tas / HP))
        else:
            aaf = 0.0

        # Elastance during early contraction using ascending elastance parameters (indices 8 and 4)
        ela = self.elastance[0, 8] + (self.elastance[1, 8] - self.elastance[0, 8]) * aaf
        era = self.elastance[0, 4] + (self.elastance[1, 4] - self.elastance[0, 4]) * aaf

        # Calculate activation function for the second elastance component (vaf) during later contraction
        if phase_fraction <= (Tas + Tav) / HP:
            vaf = 0.0
        elif phase_fraction <= (Tas + Tav + Tvs) / HP:
            # Scale the phase for the second activation phase
            vaf = np.sin(np.pi * (phase_fraction - (Tas + Tav) / HP) / (Tvs / HP))
        else:
            vaf = 0.0

        # Elastance during later contraction using descending elastance parameters (indices 9 and 5)
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
        # Compute respiratory rate
        RR = self.master_parameters['respiratory_control_params.RR_0']['value'] + y[23]  # Delta_RR_c is y[23]
        # Compute SaO2
        p_a_O2 = y[18]
        CaO2 = (self.master_parameters['params.K_O2']['value'] * np.power((1 - np.exp(-self.master_parameters['params.k_O2']['value'] * min(p_a_O2, 700))), 2)) * 100
        Sa_O2 = np.round(((CaO2 - p_a_O2 * 0.003 / 100) / (self.master_parameters['misc_constants.Hgb']['value'] * 1.34)) * 100)
        ##print(CaO2, p_a_O2, Sa_O2)
        # Store pressure value in the buffer
        self._update_pressure_buffer(P[0])

        #self.recent_MAP = self._update_pressure_buffer(P[0])
        return P, F, HR, Sa_O2, RR

    def ventilator_pressure(self, t):
        """Ventilator pressure as a function of time (simple square wave for demonstration)."""
        RR = self._ode_RR0  # Respiratory Rate (breaths per minute) - USE CACHED
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
    
    def compute_filtered_map(self,buffer):
        # fs = sample rate (Hz), cutoff = cutoff frequency (Hz)
        fs=1/self.dt
        cutoff=1
        if len(buffer) < 10:
            return np.mean(buffer)
        b, a = butter(N=2, Wn=cutoff / (0.5 * fs), btype='low')
        return filtfilt(b, a, buffer)[-1]

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
        Pbarodt, dHRv, dHRs, dRs = x[29], x[30], x[31], x[32]
        dRR_chemo, dPmus_chemo = x[33], x[34]  # New chemoreceptor states

        # ── 2) Pull “setpoint” parameters ─────────────────────────────────────────
        HR_n    = self._ode_HR_n
        R_n     = self._ode_R_n
        UV_n    = self._ode_UV_n
        RR0     = self._ode_RR0
        FI_O2   = self._ode_FI_O2
        FI_CO2  = self._ode_FI_CO2
        CO_nom  = self._ode_CO_nom
        shunt   = self._ode_shunt
        MV_mode = self._ode_MV_mode

        # Apply baroreflex deltas
        #HP = 60 / HR + dHRv + dHRh
        #HR = 60 / HP  # update HR for output tracking or logging if needed

        # Calculate desired HR from baroreflex and central control
        HR = 60/(60/HR_n + dHRv + dHRs) if self.use_reflex == True else HR_n
        #R_c = R_n + dRs if self.use_reflex == True else 1 
        R_c=1
        UV_c = UV_n
        # Apply chemoreceptor control to respiratory parameters
        RR = RR0 + dRR_chemo if self.use_reflex else RR0

        if RR < self._ode_RR0_min:
            RR = self._ode_RR0_min
        elif RR > self._ode_RR0_max:
            RR = self._ode_RR0_max
        #print(dRR_chemo)

        # Cap HR to physiological range
        if HR > self._ode_HR_n_max and self.AF == 0:
            HR = 200
        elif HR < self._ode_HR_n_min and self.AF == 0:
            HR = 30

        # Only update heart period if cardiac cycle is complete
        if self.is_cycle_complete(t):
            self.update_heart_period(HR)
        

        
        # Update respiratory rate (can change immediately)
        self.HR, self.RR = HR, RR

        # ── 3) Ventilator vs. Spontaneous Breathing ───────────────────────────────
        if MV_mode == 0:
            # spontaneous
            P_ao = 0
            Pmus_min = (
                self._ode_Pmus_0
                + Δ_Pmus_c
                + dPmus_chemo  # Add chemoreceptor contribution
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


        #print(np.sum(V),HR,P[0])
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
        # Use the new optimized buffer update function
        #self._update_pressure_buffer(P[0])


        # ── 8) Volume derivatives dV/dt ──────────────────────────────────────────
        dVdt = np.zeros(10)
        dVdt[0] = F[9] - F[0]
        for i in range(1,10):
            dVdt[i] = F[i-1] - F[i]
        # ── 7) Chemo‐ & Baroreflex updates via one call ──────────────────────────

        #print(P[0], F[0], self.elastance[0,0], Pbarodt, dHRv, dHRs, dRs)
        # Baroreflex control
        if self.use_reflex== True:
            dPbarodt, ddHRv, ddHRs, ddRs = self.baroreceptor_control(P[0], dVdt[0], self.elastance[0,0], Pbarodt, dHRv, dHRs, dRs)
        else:
            dPbarodt, ddHRv, ddHRs, ddRs = 0, 0, 0, 0
            
        # Chemoreceptor control
        if self.use_reflex == True:
            ddRR_chemo, ddPmus_chemo = self.chemoreceptor_control(p_a_CO2, p_a_O2, dRR_chemo, dPmus_chemo)
        else:
            ddRR_chemo, ddPmus_chemo = 0, 0

        # ── 9) Lung mechanics ─────────────────────────────────────────────────────
        R_lt = self._ode_R_lt
        C_cw = self.C_cw
        Ppl_dt = mech[0]/(R_lt*C_cw) - mech[1]/(C_cw*R_lt) + Pmus_dt
        dxdt_mech = (
            self.A_mechanical.dot(mech)
            + self.B_mechanical.dot([P_ao, Ppl_dt, Pmus_dt])
        )

        # ── 10) Gas exchange & alveolar ODEs ──────────────────────────────────────
        R_ml = 1.021 * 1.5
        Vdot_l = (P_ao - mech[0]) / R_ml
        R_bA   = self._ode_R_bA
        Vdot_A = (mech[2] - mech[3]) / R_bA

        p_D_CO2 = FD_CO2 * 713
        p_D_O2  = FD_O2  * 713
        FA_CO2  = p_a_CO2 / 713
        FA_O2   = p_a_O2  / 713

        K_CO2, k_CO2 = (
            self._ode_K_CO2,
            self._ode_k_CO2
        )
        c_a_CO2 = K_CO2 * p_a_CO2 + k_CO2

        if p_a_O2 > 700:
            p_a_O2 = 700
        K_O2, k_O2 = (
            self._ode_K_O2,
            self._ode_k_O2
        )
        c_a_O2 = K_O2 * (1 - np.exp(-k_O2 * p_a_O2))**2

        c_v_CO2, c_v_O2 = c_Scap_CO2, c_Scap_O2

        V_D = self._ode_V_D
        V_A = self._ode_V_A

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
        M_S_CO2 = self._ode_M_S_CO2
        D_S_CO2 = self._ode_D_S_CO2
        V_Stis_CO2 = self._ode_V_Stis_CO2
        V_Scap_CO2 = self._ode_V_Scap_CO2

        dc_Stis_CO2 = (M_S_CO2 - D_S_CO2 * (c_Stis_CO2 - c_Scap_CO2)) / V_Stis_CO2
        dc_Scap_CO2 = (
            CO_nom * 0.8 * (c_a_CO2 - c_Scap_CO2)
            + D_S_CO2 * (c_Stis_CO2 - c_Scap_CO2)
        ) / V_Scap_CO2

        M_S_O2 = self._ode_M_S_O2
        D_S_O2 = self._ode_D_S_O2
        V_Stis_O2 = self._ode_V_Stis_O2
        V_Scap_O2 = self._ode_V_Scap_O2

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
                0, 0,
                Pmus_dt,
                0, 0, 0,
                dPbarodt, ddHRv, ddHRs, ddRs,  # Baroreflex states
                ddRR_chemo, ddPmus_chemo        # NEW: Chemoreceptor states
            ]
        ])

        return dxdt

    def _update_pressure_buffer(self, new_pressure):
        """Efficiently update the pressure buffer."""
        # If the buffer is not full, simply append the new pressure
        if len(self.P_buffer) < self.window_size: # Use self.window_size
            self.P_buffer.append(new_pressure)
            self.P_buffer_sum += new_pressure  # Add to the sum when the buffer is not full
        else:
            # If the buffer is full, remove the oldest pressure and add the new one
            self.P_buffer_sum += new_pressure - self.P_buffer[0]  # Update the sum by adding the new value and removing the old one
            self.P_buffer.popleft()  # Remove the oldest pressure value
            self.P_buffer.append(new_pressure)  # Add the new pressure value

        #print(self.P_buffer_sum)  # Debugging print
    def baroreceptor_control(self, P, dVdt, elastance, Pbaro, dHRv, dHRh, dRs):
        
        # Use cached parameters for speed
        tz = self._baro_tz
        tp = self._baro_tp
        Fas_min = self._baro_Fas_min
        Fas_max = self._baro_Fas_max
        Ka = self._baro_Ka
        P_set = self._baro_P_set # This was cardio_control_params.ABP_n

        dPbarodt = (P + tz * (dVdt * elastance) - Pbaro) / tp
        Fas = (Fas_min + Fas_max * np.exp((Pbaro - P_set)/Ka)) / (1 + np.exp((Pbaro - P_set)/Ka))

        Fes_inf = self._baro_Fes_inf
        Fes_0 = self._baro_Fes_0
        Kes = self._baro_Kes
        Fev_0 = self._baro_Fev_0
        Fev_inf = self._baro_Fev_inf
        Kev = self._baro_Kev
        Fas_0 = self._baro_Fas_0

        Fes = Fes_inf + (Fes_0 - Fes_inf) * np.exp(-Kes * Fas)
        Fev = (Fev_0 + Fev_inf * np.exp((Fas - Fas_0) / Kev)) / (1 + np.exp((Fas - Fas_0) / Kev))

        self.Fes_delayed.append(max(Fes, 2.66))
        self.Fev_delayed.append(Fev)
        if len(self.Fes_delayed) > int(2/self.dt):
            self.Fes_delayed.pop(0)
        if len(self.Fev_delayed) > int(0.2/self.dt):
            self.Fev_delayed.pop(0)

        Gh = -0.13
        Ths = 2.0
        Gv = 0.09
        Thv = 1.5
        Grs = 0.45             # Baroreceptor gain resistance
        Trs = 6                     # Time constant for the resistance response to baroreceptor stimulation

        sFh = Gh * (np.log(self.Fes_delayed[0] - 2.65 + 1) - 1.1)
        sFv = Gv * (self.Fev_delayed[0] - 4.66)
        sFr = Grs * (np.log(self.Fes_delayed[-int(2/self.dt)]-2.65+1)-1.1)

        ddHRv = (sFv - dHRv) / Thv
        ddHRh = (sFh - dHRh) / Ths
        ddRs =  (sFr - dRs)/Trs

        return dPbarodt, ddHRv, ddHRh, ddRs

    def chemoreceptor_control(self, p_a_CO2, p_a_O2, dRR_chemo, dPmus_chemo):
        """
        Chemoreceptor control for respiratory rate and muscle pressure.
        Similar structure to baroreceptor_control but for respiratory variables.
        """
        # Update delayed buffers
        self.p_CO2_delayed.append(p_a_CO2)
        self.p_O2_delayed.append(p_a_O2)
        
        if len(self.p_CO2_delayed) > int(self._chemo_delay_CO2 / self.dt):
            self.p_CO2_delayed.pop(0)
        if len(self.p_O2_delayed) > int(self._chemo_delay_O2 / self.dt):
            self.p_O2_delayed.pop(0)
        
        #print(p_a_CO2, p_a_O2)
        # Get delayed values
        p_CO2_delayed = self.p_CO2_delayed[0]
        p_O2_delayed = self.p_O2_delayed[0]
        
        # CO2 chemoreceptor response (primary driver) - REDUCED SENSITIVITY
        CO2_error = p_CO2_delayed - self._chemo_p_CO2_set
        RR_drive_CO2 = self._chemo_G_CO2 * CO2_error * 0.5  # Additional 50% reduction
        
        # O2 chemoreceptor response (only active below threshold) - REDUCED SENSITIVITY
        if p_O2_delayed < self._chemo_p_O2_threshold:
            O2_error = self._chemo_p_O2_threshold - p_O2_delayed
            RR_drive_O2 = self._chemo_G_O2 * O2_error * -1 * 0.3  # Additional 70% reduction
        else:
            RR_drive_O2 = 0.0
        
        # Combined respiratory drive - FURTHER DAMPED
        total_RR_drive = RR_drive_CO2 + RR_drive_O2
        total_Pmus_drive = RR_drive_CO2 * 0.2 + RR_drive_O2 * 0.1  # Reduced Pmus scaling
        
        # First-order dynamics for chemoreceptor responses
        ddRR_chemo = (total_RR_drive - dRR_chemo) / self._chemo_tau_CO2
        ddPmus_chemo = (total_Pmus_drive - dPmus_chemo) / self._chemo_tau_O2
        
        return ddRR_chemo, ddPmus_chemo

    ## Start-stop calls
    def start_simulation(self):
        """Main simulation loop: integrates ODEs, processes events, and streams average data."""
        self.running = True
        last_print_time = self.t
        last_emit_time = self.t

        def _avg(buf, fallback):
            return np.mean(buf) if len(buf) > 0 else fallback

        while self.running:
            # Integrate ODEs
            t_span = [self.t, self.t + self.dt]
            t_eval = [self.t + self.dt]

            # Process events
            for processedEvent in list(self.events):
                if not self.processEvent(processedEvent):
                    self.events.remove(processedEvent)

            sol = solve_ivp(
                self.extended_state_space_equations,
                t_span,
                self.current_state,
                t_eval=t_eval,
                method='LSODA',
                rtol=1e-6,
                atol=1e-6
            )
            self.current_state = sol.y[:, -1]

            # Compute physiological variables
            P, F, HR, Sa_O2, RR = self.compute_variables(sol.t[-1], self.current_state)
           #print(Sa_O2)
            self.current_heart_rate = self.HR
            self.current_SaO2 = Sa_O2
            self.current_RR = self.RR

            self.current_ABP_1 = P[0]  # Link current_ABP_1 to P[3] (venous pressure compartment 3)

            # Store full-resolution waveform values
            self.P_store.append(P[0])
            self.HR_store.append(HR)

            # Store 1-second trend buffers
            self.avg_buffers["HR"].append(HR)
            self.avg_buffers["SaO2"].append(Sa_O2)
            self.avg_buffers["RR"].append(RR)
            self.avg_buffers["MAP"].append(P[0])
            #self.avg_buffers["MAP"].append(self.current_state[0])
            self.avg_buffers["etCO2"].append(self.current_state[17])

            # Estimate filtered MAP
            if len(self.P_store) == self.window_size:
                self.recent_MAP = self.compute_filtered_map(np.array(self.P_store))
            else:
                self.recent_MAP = np.mean(self.P_store) if len(self.P_store) > 0 else None

            if self.recent_MAP is None:
                self.recent_MAP = self.master_parameters['cardio_control_params.ABP_n']['value']

            # Emit average data every output_frequency seconds
            if self.t - last_emit_time >= self.output_frequency:
                avg_data = {
                    'time': self.start_timestamp + timedelta(seconds=np.round(self.t)),
                    'values': {
                        "heart_rate": np.round(_avg(self.avg_buffers["HR"], self.HR), 2),
                        "SaO2": np.round(_avg(self.avg_buffers["SaO2"], 97), 2),
                        "MAP": np.round(_avg(self.avg_buffers["MAP"], self.recent_MAP), 2),
                        "SAP": np.round(_avg(self.avg_buffers["MAP"], self.recent_MAP) + 30, 2),
                        "DAP": np.round(_avg(self.avg_buffers["MAP"], self.recent_MAP) - 15, 2),
                        "RR": np.round(_avg(self.avg_buffers["RR"], self.RR), 2),
                        "etCO2": np.round(_avg(self.avg_buffers["etCO2"], self.current_state[17]), 2)
                    }
                }

                self.redisClient.add_vital_sign(self.patient_id, avg_data)

                if not isinstance(self.data_epoch, list):
                    self.data_epoch = []

                self.data_epoch.append(avg_data)
                self.data_epoch = self.data_epoch[-self.data_points:]

                self.alarmModule.evaluate_data(curr_data=avg_data, historic_data=self.data_epoch)
                last_emit_time = self.t
                #print(avg_data)
            # Periodic logging (optional)
            if self.t - last_print_time >= 5:
                last_print_time = self.t

            # Advance time
            self.t += self.dt
            if self.sleep:
                time.sleep(self.dt)

    def addProcessedEvent(self, event):
        """Add a processed event to the event queue."""
        if event:
            self.events.append(event)
            print(f"Event added for patient {self.patient_id}: {event['eventType']}")
    def processEvent(self, eventContent):
        outcome = None
        error = None
        parameters_changed_by_event = [] # Keep track of parameters changed
        ode_parameters_changed = False # Flag to check if ODE specific parameters changed
        respi_constants_changed = False # NEW: Flag for respiratory constant changes
        
        ##print(f"Started processing: {eventContent}")
        """Process an event and update the model parameters."""
        if eventContent["eventType"] == "common": ## routine, no specific things -> Special events re-route to custom functions
            if eventContent["timeCategorical"] == "continuous" or eventContent["timeCategorical"] == "limited":
                ## ongoing event, check delta T 
                outcome = True ## keep the event 
                lastEmissionTime = eventContent["lastEmission"] if eventContent["lastEmission"] else 0
                if self.t - lastEmissionTime <= int(timedelta(**{eventContent["timeUnit"]: eventContent["timeInterval"]}).total_seconds()):
                    ## Not enough time has passed: Break early
                    return outcome # Return outcome directly
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
                newValue = None # Initialize newValue
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
                            error = f"Error: Unknown action {paramChange['action']} for parameter {paramName} in patient {self.patient_id}"
                            break # break from inner loop

                        ## get the actual value for the parameter in non-percentualness
                        splineFunction = self.pathologies.splineFunctions[paramName]
                        newValue = splineFunction(currPercent)
                        ## Update the parameter in the model
                        self.master_parameters[paramName]['value'] = newValue
                        parameters_changed_by_event.append(paramName)
                        print(f" === EVENT ACTION: param {paramName} updated to {self.master_parameters[paramName]['value']} ===")
                        if paramName.startswith(('cardio_control_params.', 'respiratory_control_params.', 'gas_exchange_params.', 'bloodflows.', 'misc_constants.', 'params.')):
                            ode_parameters_changed = True
                        if paramName.startswith('respi_constants.'):
                            respi_constants_changed = True

                    ## B] Type = 'absolute' -> set directly, no inference of spline
                    elif paramChange['type'] == 'absolute':
                        if paramChange['action'] == 'decay':
                            ## decay the parameter
                            self.master_parameters[paramName]['value'] += paramChange['value']
                            newValue = self.master_parameters[paramName]['value']
                        elif paramChange['action'] == 'set':
                            ## set parameter direcly
                            self.master_parameters[paramName]['value'] = paramChange['value']
                            newValue = self.master_parameters[paramName]['value']
                        else:
                            error = f"Error: Unknown action {paramChange['action']} for parameter {paramName} in patient {self.patient_id}"
                            break # break from inner loop
                        parameters_changed_by_event.append(paramName)
                        if paramName.startswith(('cardio_control_params.', 'respiratory_control_params.', 'gas_exchange_params.', 'bloodflows.', 'misc_constants.', 'params.')):
                            ode_parameters_changed = True
                        if paramName.startswith('respi_constants.'):
                            respi_constants_changed = True
                    else:
                        error = f"Error: Unknown type {paramChange['type']} for parameter {paramName} in patient {self.patient_id}"
                        break # break from inner loop
                
                else:
                    ## Parameter not in the master-parameters file, try manual CASE-MATCH
                    paramNameNew = None # Initialize to avoid UnboundLocalError
                    match paramName:
                        case "venous_compartiment_3":
                            paramNameNew = "self.current_state[3]"
                        ## match future parameters

                    if not paramNameNew:
                        error = f"Error: Unknown parameter {paramName} in patient {self.patient_id}"
                        break # break from inner loop
                    ## access the variable through the eval function
                    paramVariable = eval(paramNameNew) # Be cautious with eval
                    #print(f"+++ ParamChange = {paramChange} +++")
                    if paramChange['type'] == 'absolute':
                        if paramChange['action'] == 'decay':
                            paramVariable += paramValChange
                            exec(f"{paramNameNew} = {paramVariable}") # Be cautious with exec
                            #print(f" === EVENT ACTION: param {paramName} updated to {paramVariable} ===")
                        elif paramChange['action'] == 'set':
                            paramVariable = paramValChange 
                            exec(f"{paramNameNew} = {paramVariable}") # Be cautious with exec
                            #print(f" === EVENT ACTION: param {paramName} updated to {paramVariable} ===")
                        else:
                            error = f"Error: Unknown action {paramChange['action']} for parameter {paramName} in patient {self.patient_id}"
                            break # break from inner loop
                    else:
                        error = f"Error: Unknown type {paramChange['type']} for parameter {paramName} in patient {self.patient_id}"
                        break # break from inner loop
            
            if error: # If an error occurred in the loop
                print(error) # Or handle error more robustly
                return False # Indicate event processing failed or should be removed

            # If baroreflex parameters were changed by an event, re-cache them
            if any(p_name.startswith('baroreflex_params.') or p_name == 'cardio_control_params.ABP_n' for p_name in parameters_changed_by_event):
                self._cache_baroreflex_parameters()
                print("Re-cached baroreflex parameters due to event.")
                
            # If chemoreceptor parameters were changed by an event, re-cache them
            if any(p_name.startswith('chemoreceptor_params.') for p_name in parameters_changed_by_event):
                self._cache_chemoreceptor_parameters()
                print("Re-cached chemoreceptor parameters due to event.")

            # If ODE parameters were changed by an event (or implied by respi_constants change), re-cache them
            if ode_parameters_changed:
                self._cache_ode_parameters()
                print("Re-cached ODE parameters due to event.")

            ## recalculate cardiac elastances, resistance and uVolme if relevant parameters changed
            # This check can be made more specific if needed.
            # Note: compute_cardiac_parameters() updates self.elastance, self.resistance, self.uvolume.
            # If self.uvolume changes, self.initialize_state() might need to be called if initial volumes depend on it,
            # but that's a larger reset not currently handled here.
            if any(p_name.startswith('cardio.') for p_name in parameters_changed_by_event): # A more specific check for cardiac parameters
                self.compute_cardiac_parameters()
                print("Re-computed cardiac parameters due to event.")
        
        elif eventContent["eventType"] == "special": 
            ## Custom-made function for special events
            outcome = False
            print(f"Special events not yet incorporated, removing from list")
            ## to be filled later!

        return outcome 
     
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