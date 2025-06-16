import threading
import time
from datetime import datetime
import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from collections import deque

from digital_twin_model import DigitalTwinModel

class DummyRedisClient:
    def add_vital_sign(self, patient_id, data):
        pass
    def add_alarm_entry(self, patient_id, alarm_payload):
        print(f"Mock Redis: Alarm for {patient_id}: {alarm_payload.get('alarm_msg')}")

class ParameterControlWindow:
    def __init__(self, parent, twin_model):
        self.parent = parent
        self.twin_model = twin_model
        self.window = None
        self.current_temp = 37.2  # Add temperature storage
        
    def open_window(self):
        if self.window is not None:
            self.window.lift()
            return
            
        self.window = tk.Toplevel(self.parent)
        self.window.title("Parameter Controls")
        self.window.geometry("350x600")  # Increased height for temperature control
        self.window.configure(bg='#2b2b2b')
        self.window.protocol("WM_DELETE_WINDOW", self._on_close)
        
        # Title
        title_label = tk.Label(self.window, text="Physiological Parameters", 
                              font=('Helvetica', 16, 'bold'), fg='white', bg='#2b2b2b')
        title_label.pack(pady=10)
        
        # Get initial values from the model
        initial_fio2 = self._get_model_param('gas_exchange_params.FI_O2', 0.21) * 100
        initial_tbv = self._get_model_param('misc_constants.TBV', 5000.0)  # TBV is in ml
        initial_hr_n = self._get_model_param('cardio_control_params.HR_n', 72)
        initial_rr_0 = self._get_model_param('respiratory_control_params.RR_0', 15)
        
        # FiO2 Control
        self._create_parameter_control("FiO2 (%)", "fio2", 0, 100, initial_fio2, self._update_fio2)
        
        # TBV Control - now in ml  
        self._create_parameter_control("TBV (ml)", "tbv", 3000.0, 8000.0, initial_tbv, self._update_tbv)
        
        # HR_n Control
        self._create_parameter_control("HR Normal (bpm)", "hr_n", 40, 120, initial_hr_n, self._update_hr_n)
        
        # RR_0 Control
        self._create_parameter_control("RR Baseline (bpm)", "rr_0", 8, 40, initial_rr_0, self._update_rr_0)
        
        # Temperature Control
        self._create_parameter_control("Temperature (°C)", "temp", 34.0, 42.0, self.current_temp, self._update_temp)
        
        # Current values display
        values_frame = tk.Frame(self.window, bg='#2b2b2b')
        values_frame.pack(fill=tk.X, padx=20, pady=20)
        
        tk.Label(values_frame, text="Current Values:", font=('Helvetica', 12, 'bold'), 
                fg='white', bg='#2b2b2b').pack(anchor='w')
        
        self.current_values_label = tk.Label(values_frame, text="", font=('Helvetica', 10), 
                                           fg='#cccccc', bg='#2b2b2b', justify='left')
        self.current_values_label.pack(anchor='w', pady=5)
        
        # Reset button
        reset_button = tk.Button(self.window, text="Reset to Defaults", 
                                command=self._reset_parameters,
                                bg='#555555', fg='white', font=('Helvetica', 10))
        reset_button.pack(pady=10)
        
        self._update_current_values()
        
    def _get_model_param(self, param_name, default_value):
        """Get parameter value from the twin model."""
        try:
            return self.twin_model.master_parameters[param_name]['value']
        except (KeyError, AttributeError):
            return default_value
            
    def _create_event_for_parameter(self, param_name, new_value):
        """Create an event to update a parameter using the twin model's event system."""
        event = {
            "eventType": "common",
            "timeCategorical": "instant",
            "timeUnit": "seconds",
            "timeInterval": 0,
            "eventCount": 1,
            "lastEmission": None,
            "parameters": [
                {
                    "name": param_name,
                    "type": "absolute",
                    "action": "set",
                    "value": new_value
                }
            ]
        }
        self.twin_model.addProcessedEvent(event)
        
    def _create_parameter_control(self, label, param_name, min_val, max_val, default_val, callback):
        frame = tk.Frame(self.window, bg='#2b2b2b')
        frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Label
        tk.Label(frame, text=label, font=('Helvetica', 12), fg='white', bg='#2b2b2b').pack(anchor='w')
        
        # Scale and value display frame
        control_frame = tk.Frame(frame, bg='#2b2b2b')
        control_frame.pack(fill=tk.X, pady=5)
        
        # Value display
        value_var = tk.StringVar(value=f"{default_val:.1f}")
        value_label = tk.Label(control_frame, textvariable=value_var, font=('Helvetica', 11, 'bold'), 
                              fg='yellow', bg='#2b2b2b', width=8)
        value_label.pack(side=tk.RIGHT)
        
        # Scale
        scale = tk.Scale(control_frame, from_=min_val, to=max_val, orient=tk.HORIZONTAL,
                        bg='#2b2b2b', fg='white', highlightbackground='#2b2b2b',
                        activebackground='#444444', troughcolor='#444444',
                        command=lambda val, p=param_name, v=value_var, c=callback: self._on_scale_change(val, p, v, c))
        scale.set(default_val)
        scale.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        # Store reference
        setattr(self, f"{param_name}_var", value_var)
        setattr(self, f"{param_name}_scale", scale)
        
    def _on_scale_change(self, value, param_name, value_var, callback):
        float_val = float(value)
        value_var.set(f"{float_val:.1f}")
        callback(float_val)
        self._update_current_values()
        
    def _update_fio2(self, value):
        # Update FiO2 using the event system
        fio2_fraction = value / 100.0
        self._create_event_for_parameter('gas_exchange_params.FI_O2', fio2_fraction)
        print(f"FiO2 event created: {value}% (fraction: {fio2_fraction})")
        
    def _update_tbv(self, value):
        # Update TBV using the event system - value is already in ml
        self._create_event_for_parameter('misc_constants.TBV', value)
        print(f"TBV event created: {value}ml")
        
        # Also need to update blood volumes when TBV changes
        # This requires a special handling since it affects the state vector
        try:
            # Recalculate blood volumes based on new TBV and current uvolume distribution
            # value is already in ml, so use it directly
            new_volumes = value * (self.twin_model.uvolume / np.sum(self.twin_model.uvolume))
            # Update current blood volumes (indices 0-9 in the state vector)
            self.twin_model.current_state[:10] = new_volumes
            print(f"Blood volumes recalculated for new TBV: {value}ml")
        except Exception as e:
            print(f"Warning: Could not update blood volumes: {e}")
            
    def _update_hr_n(self, value):
        # Update HR_n using the event system
        self._create_event_for_parameter('cardio_control_params.HR_n', value)
        print(f"HR_n event created: {value} bpm")
            
    def _update_rr_0(self, value):
        # Update RR_0 using the event system
        self._create_event_for_parameter('respiratory_control_params.RR_0', value)
        print(f"RR_0 event created: {value} bpm")
            
    def _update_temp(self, value):
        # Update temperature (stored locally, not in model parameters)
        self.current_temp = value
        print(f"Temperature updated: {value}°C")
            
    def _reset_parameters(self):
        """Reset parameters to their default values using events."""
        # Reset FiO2 to 21%
        self.fio2_scale.set(21)
        self._update_fio2(21)
        
        # Reset TBV to 5000ml (5L)
        self.tbv_scale.set(5000.0)
        self._update_tbv(5000.0)
        
        # Reset HR_n to 72 bpm
        self.hr_n_scale.set(72)
        self._update_hr_n(72)
        
        # Reset RR_0 to 15 bpm
        self.rr_0_scale.set(15)
        self._update_rr_0(15)
        
        # Reset Temperature to 37.2°C
        self.temp_scale.set(37.2)
        self._update_temp(37.2)
        
        print("Parameter reset events created")
        
    def _update_current_values(self):
        if self.window is None:
            return
            
        # Get current values from twin model
        fio2_val = self._get_model_param('gas_exchange_params.FI_O2', 0.21) * 100
        tbv_val = self._get_model_param('misc_constants.TBV', 5000.0)  # TBV is in ml
        hr_n_val = self._get_model_param('cardio_control_params.HR_n', 72)
        rr_0_val = self._get_model_param('respiratory_control_params.RR_0', 15)
        
        values_text = f"FiO2: {fio2_val:.1f}%\nTBV: {tbv_val:.0f}ml\nHR_n: {hr_n_val:.0f} bpm\nRR_0: {rr_0_val:.0f} bpm\nTemp: {self.current_temp:.1f}°C"
        self.current_values_label.config(text=values_text)
        
        # Schedule next update
        if self.window:
            self.window.after(1000, self._update_current_values)
        
    def _on_close(self):
        self.window.destroy()
        self.window = None

class MonitorApp(tk.Tk):
    def __init__(self, twin_model):
        super().__init__()
        self.twin_model = twin_model
        self.title("Digital Twin ICU Monitor")
        self.configure(bg='black')
        self.geometry("1400x850")
        self._ecg_time = 0.0  # keeps track of global ECG time in seconds
        self._ecg_phase = 0.0  # in seconds, resets every beat

        self.running = False
        self.sim_thread = None
        self.gui_update_interval_ms = 100

        self.plot_window_seconds = 30
        self.data_points_per_plot = int(self.plot_window_seconds * (1000 / self.gui_update_interval_ms))

        self.time_data = deque(maxlen=self.data_points_per_plot)
        self.hr_data = deque(maxlen=self.data_points_per_plot)
        self.map_data = deque(maxlen=self.data_points_per_plot)
        self.spo2_data = deque(maxlen=self.data_points_per_plot)
        self.rr_data = deque(maxlen=self.data_points_per_plot)

        self.ecg_display_samples_per_update = 10
        self.ecg_total_display_points = self.data_points_per_plot * self.ecg_display_samples_per_update
        self.ecg_x_data = deque(maxlen=self.ecg_total_display_points)
        self.ecg_y_data = deque(maxlen=self.ecg_total_display_points)

        # Change EtCO₂ to use Pmus data with same sampling as ECG
        self.pmus_display_samples_per_update = self.ecg_display_samples_per_update
        self.pmus_total_display_points = self.data_points_per_plot * self.pmus_display_samples_per_update
        self.pmus_waveform_x_data = deque(maxlen=self.pmus_total_display_points)
        self.pmus_waveform_y_data = deque(maxlen=self.pmus_total_display_points)
        self._pmus_time = 0.0  # Add Pmus time tracking similar to ECG
        self.map_buffer = deque(maxlen=int(10 * (1000 / self.gui_update_interval_ms)))  # Store 10 seconds of map_val
        
        # Initialize parameter control window
        self.param_control = ParameterControlWindow(self, self.twin_model)
        
        self._setup_gui()
        self.protocol("WM_DELETE_WINDOW", self._on_closing)

    def _setup_gui(self):
        main_frame = tk.Frame(self, bg='black')
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Top patient info bar
        header = tk.Frame(main_frame, bg='#444444', height=40)
        header.pack(side=tk.TOP, fill=tk.X)
        self.patient_info_label = tk.Label(header, text="Bed 3   Doe, John   Adult   " + datetime.now().strftime("%H:%M"),
                                           font=('Helvetica', 14, 'bold'), fg='white', bg='#444444', anchor='w')
        self.patient_info_label.pack(side=tk.RIGHT, padx=10)

        control_frame = tk.Frame(main_frame, bg='black')
        control_frame.pack(side=tk.TOP, fill=tk.X)

        self.start_button = ttk.Button(control_frame, text="Start Simulation", command=self.start_simulation_action)
        self.start_button.pack(side=tk.LEFT, padx=5)

        self.stop_button = ttk.Button(control_frame, text="Stop Simulation", command=self.stop_simulation_action, state=tk.DISABLED)
        self.stop_button.pack(side=tk.LEFT, padx=5)
        
        # Add parameter control button
        self.param_button = ttk.Button(control_frame, text="Parameters", command=self.param_control.open_window)
        self.param_button.pack(side=tk.LEFT, padx=5)

        content_frame = tk.Frame(main_frame, bg='black')
        content_frame.pack(fill=tk.BOTH, expand=True)

        left_frame = tk.Frame(content_frame, bg='black')
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        right_frame = tk.Frame(content_frame, bg='black')
        right_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=20)

        self.fig = Figure(figsize=(10, 6), dpi=100, facecolor='black')
        self.canvas = FigureCanvasTkAgg(self.fig, master=left_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        signals = [
            ("ECG", self.ecg_x_data, self.ecg_y_data, 'lime', (-1.5, 1.5)),
            ("Art (mmHg)", self.time_data, self.map_data, 'red', (50, 160)),
            ("SpO₂ (%)", self.time_data, self.spo2_data, 'cyan', (140, 250)),
            ("EtCO₂ (mmHg)", self.pmus_waveform_x_data, self.pmus_waveform_y_data, 'yellow', (0, 60))  # Keep EtCO₂ label but use Pmus data
        ]

        self.axes, self.lines = {}, {}
        for i, (label, x_data, y_data, color, ylim) in enumerate(signals):
            ax = self.fig.add_subplot(len(signals), 1, i + 1)
            line, = ax.plot(x_data, y_data, color=color, lw=2)
            ax.set_ylim(ylim)
            ax.set_facecolor('black')
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.set_title(label, color='white', fontsize=10, pad=5, loc='left')
            self.axes[label] = ax
            self.lines[label] = line

        # Store reference to EtCO₂ fill for updating (keep original name)
        self.etco2_fill = None

        self.hr_var = tk.StringVar(value="--")
        self.map_var = tk.StringVar(value="--")
        self.map_dia_var = tk.StringVar(value="(--)")
        self.spo2_var = tk.StringVar(value="--")
        self.rr_var = tk.StringVar(value="--")
        self.etco2_var = tk.StringVar(value="--")
        self.temp_var = tk.StringVar(value="--")  # Add temperature variable

        def big_label(parent, textvar, color, size):
            return tk.Label(parent, textvariable=textvar, font=('Helvetica', size, 'bold'), fg=color, bg='black', anchor='w', width=6)

        def title_label(parent, text, color='white'):
            return tk.Label(parent, text=text, font=('Helvetica', 16), fg=color, bg='black', anchor='w')

        for label_text, var, color, size in [
            ("HR", self.hr_var, 'lime', 60),
            ("MAP", self.map_var, 'red', 60),
            ("", self.map_dia_var, 'red', 40),
            ("SpO₂", self.spo2_var, 'cyan', 60),
            ("RR", self.rr_var, 'white', 60),
            ("EtCO₂", self.etco2_var, 'yellow', 60),
            ("Temp", self.temp_var, 'orange', 60)  # Add temperature display
        ]:
            if label_text:
                title_label(right_frame, label_text).pack(anchor='w')
            big_label(right_frame, var, color, size).pack(anchor='w', pady=(0, 20))

        self.fig.tight_layout(pad=2)

    def start_simulation_action(self):
        if not self.running:
            self.running = True
            self.start_button.config(state=tk.DISABLED)
            self.stop_button.config(state=tk.NORMAL)
            self.twin_model.start_timestamp = datetime.now()
            self.twin_model.t = 0
            # Initialize time bases
            current_time = time.time()
            self._ecg_time = current_time
            self._pmus_time = current_time  # Initialize Pmus time
            self.sim_thread = threading.Thread(target=self.twin_model.start_simulation, daemon=True)
            self.sim_thread.start()
            self._update_gui()

    def stop_simulation_action(self):
        if self.running:
            self.running = False
            self.start_button.config(state=tk.NORMAL)
            self.stop_button.config(state=tk.DISABLED)
            self.twin_model.stop_simulation()

    def _generate_etco2_waveform(self, plateau_co2, rr_value, start_time):
        num_new_points = self.etco2_display_samples_per_update
        t_step = self.gui_update_interval_ms / 1000.0 / num_new_points
        new_times = np.linspace(start_time, start_time + (self.gui_update_interval_ms / 1000), num_new_points)
        waveform = np.zeros(num_new_points)

        if plateau_co2 is None or rr_value is None or rr_value <= 0:
            return new_times, waveform

        T_resp = 60.0 / rr_value
        T_insp = T_resp / 3
        T_exp = T_resp * 2 / 3
        baseline = 2.0

        for i, t_abs in enumerate(new_times):
            pos = (t_abs % T_resp)
            if pos < T_insp:
                waveform[i] = baseline
            else:
                phase = pos - T_insp
                if phase < 0.15 * T_exp:
                    waveform[i] = baseline
                elif phase < 0.5 * T_exp:
                    progress = (phase - 0.15 * T_exp) / (0.35 * T_exp)
                    waveform[i] = baseline + (plateau_co2 - baseline) * np.sin(progress * np.pi / 2)
                else:
                    waveform[i] = plateau_co2
        noise = np.random.normal(0, 0.2, num_new_points)
        return new_times, np.clip(waveform + noise, 0, 60)
    
    def _generate_ecg_waveform(self, t, heart_rate):
        if heart_rate is None or heart_rate < 30:
            heart_rate = 60

        cycle_duration = 60.0 / heart_rate
        p_wave_duration = 0.1 * cycle_duration
        qrs_duration = 0.1 * cycle_duration
        t_wave_duration = 0.2 * cycle_duration

        p_amp = 0.15
        qrs_amp = 1.2
        t_amp = 0.35

        ecg = np.zeros_like(t)
        t_in_cycle = t % cycle_duration

        # P-wave
        p_start = 0.1 * cycle_duration
        p_end = p_start + p_wave_duration
        p_mask = (t_in_cycle >= p_start) & (t_in_cycle < p_end)
        ecg[p_mask] += p_amp * np.sin(np.pi * (t_in_cycle[p_mask] - p_start) / p_wave_duration)

        # QRS complex
        qrs_start = 0.3 * cycle_duration
        qrs_end = qrs_start + qrs_duration
        qrs_mask = (t_in_cycle >= qrs_start) & (t_in_cycle < qrs_end)
        ecg[qrs_mask] += qrs_amp * np.sin(np.pi * (t_in_cycle[qrs_mask] - qrs_start) / qrs_duration)

        # T-wave
        t_start = 0.5 * cycle_duration
        t_end = t_start + t_wave_duration
        t_mask = (t_in_cycle >= t_start) & (t_in_cycle < t_end)
        ecg[t_mask] += t_amp * np.sin(np.pi * (t_in_cycle[t_mask] - t_start) / t_wave_duration)

        return ecg


    def _autoscale_axis(self, data, default_min, default_max, padding_factor=0.1):
        """Calculate autoscaled y-axis limits for given data"""
        valid_data = [val for val in data if not np.isnan(val)]
        if not valid_data:
            return default_min, default_max
        
        data_min, data_max = min(valid_data), max(valid_data)
        data_range = data_max - data_min
        
        if data_range == 0:
            # If all values are the same, add some padding around the value
            padding = max(10, abs(data_min) * 0.1)
            return data_min - padding, data_max + padding
        
        padding = data_range * padding_factor
        return data_min - padding, data_max + padding

    def _update_gui(self):
        if not self.running:
            return

        current_time = time.time()
        hr = getattr(self.twin_model, "current_heart_rate", None)
        map_val = getattr(self.twin_model, "recent_MAP", None)
        spo2 = getattr(self.twin_model, "current_SaO2", None)
        spo2_v2 = getattr(self.twin_model, "current_ABP_1", None)
        rr = getattr(self.twin_model, "current_RR", None)
        # Get Pmus value from the digital twin's state vector (index 25) for waveform
        pmus_val = self.twin_model.current_state[25] if hasattr(self.twin_model, "current_state") and len(self.twin_model.current_state) > 25 else None
        # Get EtCO₂ value from the digital twin's state vector (index 17) for numerical display
        etco2_val = self.twin_model.current_state[17] if hasattr(self.twin_model, "current_state") and len(self.twin_model.current_state) > 17 else None

        # Append map_val to the buffer
        self.map_buffer.append(map_val if map_val is not None else np.nan)

        # Calculate the mean of the last 10 seconds of map_val
        valid_map_values = [val for val in self.map_buffer if not np.isnan(val)]
        mean_map_val = np.mean(valid_map_values) if valid_map_values else None

        self.time_data.append(current_time)
        self.hr_data.append(hr or np.nan)
        self.map_data.append(map_val or np.nan)
        # Use MAP-derived SpO₂ values (back to original)
        #self.spo2_data_real.append(spo2 or np.nan)
        self.spo2_data.append(spo2_v2 or np.nan)
        self.rr_data.append(rr or np.nan)

        # Generate ECG based on internal continuous ECG timebase
        t_step = self.gui_update_interval_ms / 1000.0 / self.ecg_display_samples_per_update
        new_ecg_x = np.linspace(self._ecg_time, self._ecg_time + t_step * self.ecg_display_samples_per_update, self.ecg_display_samples_per_update)
        new_ecg_y = self._generate_ecg_waveform(new_ecg_x, 72)
        self._ecg_time += t_step * self.ecg_display_samples_per_update

        self.ecg_x_data.extend(new_ecg_x)
        self.ecg_y_data.extend(new_ecg_y)

        # Generate Pmus waveform using actual Pmus values from the digital twin
        pmus_t_step = self.gui_update_interval_ms / 1000.0 / self.pmus_display_samples_per_update
        new_pmus_x = np.linspace(self._pmus_time, self._pmus_time + pmus_t_step * self.pmus_display_samples_per_update, self.pmus_display_samples_per_update)
        # Use actual Pmus value for all samples in this update (real-time sampling)
        new_pmus_y = np.full(self.pmus_display_samples_per_update, pmus_val if pmus_val is not None else 0.0)
        self._pmus_time += pmus_t_step * self.pmus_display_samples_per_update
        
        self.pmus_waveform_x_data.extend(new_pmus_x)
        self.pmus_waveform_y_data.extend(new_pmus_y)

        # Update displayed values
        self.hr_var.set(f"{hr:.0f}" if hr else "--")
        self.map_var.set(f"{mean_map_val:.0f}" if mean_map_val is not None else "--")
        self.spo2_var.set(f"{spo2:.0f}" if spo2 else "--")
        self.rr_var.set(f"{rr:.0f}" if rr else "--")
        self.etco2_var.set(f"{etco2_val:.0f}" if etco2_val is not None else "--")
        self.temp_var.set(f"{self.param_control.current_temp:.1f}")  # Use temperature from parameter control

        if mean_map_val:
            sap, dap = (mean_map_val + 0.3 * mean_map_val), (mean_map_val - 0.2 * mean_map_val)
            self.map_dia_var.set(f"({sap:.0f}/{dap:.0f})")
        else:
            self.map_dia_var.set("(--)")

        self.lines["ECG"].set_data(list(self.ecg_x_data), list(self.ecg_y_data))
        self.axes["ECG"].set_xlim(self._ecg_time - self.plot_window_seconds, self._ecg_time)

        self.lines["Art (mmHg)"].set_data(list(self.time_data), list(self.spo2_data))
        self.axes["Art (mmHg)"].set_xlim(current_time - self.plot_window_seconds, current_time)
        # Add autoscaling for MAP
        map_min, map_max = self._autoscale_axis(list(self.spo2_data), 50, 160)
        self.axes["Art (mmHg)"].set_ylim(map_min, map_max)

        # Create derived SpO₂ plot data but with independent scaling
        #spo2_plot_data = [val * 2.2 for val in self.map_data]
        self.lines["SpO₂ (%)"].set_data(list(self.time_data), self.map_data)
        self.axes["SpO₂ (%)"].set_xlim(current_time - self.plot_window_seconds, current_time)
        # Independent autoscaling for SpO₂ with aggressive scaling to fill y-axis
        spo2_min, spo2_max = self._autoscale_axis(self.map_data, 140, 250, padding_factor=0.02)
        # Ensure minimum range for visual distinction
        if spo2_max - spo2_min < 10:
            mid_point = (spo2_max + spo2_min) / 2
            spo2_min = mid_point - 5
            spo2_max = mid_point + 5
        self.axes["SpO₂ (%)"].set_ylim(spo2_min, spo2_max)

        # Update EtCO₂ plot (keep original label but use Pmus data)
        self.lines["EtCO₂ (mmHg)"].set_data(list(self.pmus_waveform_x_data), list(self.pmus_waveform_y_data))
        self.axes["EtCO₂ (mmHg)"].set_xlim(self._pmus_time - self.plot_window_seconds, self._pmus_time)
        
        # Add dynamic y-axis scaling for Pmus data
        if len(self.pmus_waveform_y_data) > 0:
            pmus_min, pmus_max = self._autoscale_axis(list(self.pmus_waveform_y_data), -15, 5, padding_factor=0.1)
            # Ensure minimum range for visual distinction
            if pmus_max - pmus_min < 2:
                mid_point = (pmus_max + pmus_min) / 2
                pmus_min = mid_point - 1
                pmus_max = mid_point + 1
            self.axes["EtCO₂ (mmHg)"].set_ylim(pmus_min, pmus_max)
        else:
            self.axes["EtCO₂ (mmHg)"].set_ylim(-15, 5)  # Default range

        # Update EtCO₂ area fill (keep original name, fix fill direction)
        if self.etco2_fill is not None:
            self.etco2_fill.remove()
        if len(self.pmus_waveform_x_data) > 0 and len(self.pmus_waveform_y_data) > 0:
            # Get the minimum y-value to fill from the bottom of the curve
            min_y = min(self.pmus_waveform_y_data) if self.pmus_waveform_y_data else 0
            self.etco2_fill = self.axes["EtCO₂ (mmHg)"].fill_between(
                list(self.pmus_waveform_x_data), 
                min_y,  # Fill from bottom
                list(self.pmus_waveform_y_data),  # Fill to curve
                alpha=0.3, 
                color='yellow'
            )

        try:
            self.canvas.draw_idle()
        except tk.TclError:
            pass

        self.after(self.gui_update_interval_ms, self._update_gui)

    def _on_closing(self):
        if self.running:
            self.stop_simulation_action()
        self.destroy()

def main():
    twin = DigitalTwinModel(patient_id="demo_patient_gui", param_file="healthyFlat_leukemia.json", sleep=True, time_step=0.02)
    if not hasattr(twin, 'redisClient') or twin.redisClient is None:
        twin.redisClient = DummyRedisClient()
    app = MonitorApp(twin)
    app.mainloop()

if __name__ == "__main__":
    main()