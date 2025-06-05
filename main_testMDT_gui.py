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

        self.etco2_display_samples_per_update = self.ecg_display_samples_per_update
        self.etco2_total_display_points = self.data_points_per_plot * self.etco2_display_samples_per_update
        self.etco2_waveform_x_data = deque(maxlen=self.etco2_total_display_points)
        self.etco2_waveform_y_data = deque(maxlen=self.etco2_total_display_points)
        self.map_buffer = deque(maxlen=int(10 * (1000 / self.gui_update_interval_ms)))  # Store 10 seconds of map_val
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
            ("MAP (mmHg)", self.time_data, self.map_data, 'red', (50, 160)),
            ("SpO₂ (%)", self.time_data, self.spo2_data, 'cyan', (140, 250)),
            ("EtCO₂ (mmHg)", self.etco2_waveform_x_data, self.etco2_waveform_y_data, 'yellow', (0, 60))
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

        self.hr_var = tk.StringVar(value="--")
        self.map_var = tk.StringVar(value="--")
        self.map_dia_var = tk.StringVar(value="(--)")
        self.spo2_var = tk.StringVar(value="--")
        self.rr_var = tk.StringVar(value="--")
        self.etco2_var = tk.StringVar(value="--")

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
            ("EtCO₂", self.etco2_var, 'yellow', 60)
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
            self.sim_thread = threading.Thread(target=self.twin_model.start_simulation, daemon=True)
            self.sim_thread.start()
            self._update_gui()

    def stop_simulation_action(self):
        if self.running:
            self.running = False
            self.start_button.config(state=tk.NORMAL)
            self.stop_button.config(state=tk.DISABLED)
            self.twin_model.stop_simulation()

    def _generate_etco2_waveform(self, plateau_co2, rr_value):
        num_new_points = self.etco2_display_samples_per_update
        last_time = self.etco2_waveform_x_data[-1] if self.etco2_waveform_x_data else time.time()
        new_times = np.linspace(last_time, last_time + (self.gui_update_interval_ms / 1000), num_new_points)
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
        noise = np.random.normal(0, 0.5, num_new_points)
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


    def _update_gui(self):
        if not self.running:
            return

        current_time = time.time()
        hr = getattr(self.twin_model, "current_heart_rate", None)
        map_val = getattr(self.twin_model, "recent_MAP", None)
        spo2 = getattr(self.twin_model, "current_SaO2", None)
        spo2_v2 = getattr(self.twin_model, "current_ABP_1", None)
        rr = getattr(self.twin_model, "current_RR", None)
        etco2_val = self.twin_model.current_state[17] if hasattr(self.twin_model, "current_state") else None

        # Append map_val to the buffer
        self.map_buffer.append(map_val if map_val is not None else np.nan)

        # Calculate the mean of the last 10 seconds of map_val
        valid_map_values = [val for val in self.map_buffer if not np.isnan(val)]
        mean_map_val = np.mean(valid_map_values) if valid_map_values else None

        self.time_data.append(current_time)
        self.hr_data.append(hr or np.nan)
        self.map_data.append(map_val or np.nan)
        self.spo2_data.append(spo2_v2 or np.nan)
        self.rr_data.append(rr or np.nan)

        # Generate ECG based on internal continuous ECG timebase
        t_step = self.gui_update_interval_ms / 1000.0 / self.ecg_display_samples_per_update
        new_ecg_x = np.linspace(self._ecg_time, self._ecg_time + t_step * self.ecg_display_samples_per_update, self.ecg_display_samples_per_update)
        new_ecg_y = self._generate_ecg_waveform(new_ecg_x, 72)
        self._ecg_time += t_step * self.ecg_display_samples_per_update


        self.ecg_x_data.extend(new_ecg_x)
        self.ecg_y_data.extend(new_ecg_y)

        new_etco2_x, new_etco2_y = self._generate_etco2_waveform(etco2_val, rr)
        self.etco2_waveform_x_data.extend(new_etco2_x)
        self.etco2_waveform_y_data.extend(new_etco2_y)

        # Update displayed values
        self.hr_var.set(f"{hr:.0f}" if hr else "--")
        self.map_var.set(f"{mean_map_val:.0f}" if mean_map_val is not None else "--")  # Use mean_map_val
        self.spo2_var.set(f"{spo2:.0f}" if spo2 else "--")
        self.rr_var.set(f"{rr:.0f}" if rr else "--")
        self.etco2_var.set(f"{etco2_val:.0f}" if etco2_val is not None else "--")

        if mean_map_val:
            sap, dap = (mean_map_val + 0.3 * mean_map_val), (mean_map_val - 0.2 * mean_map_val)
            self.map_dia_var.set(f"({sap:.0f}/{dap:.0f})")
        else:
            self.map_dia_var.set("(--)")

        self.lines["ECG"].set_data(list(self.ecg_x_data), list(self.ecg_y_data))
        self.axes["ECG"].set_xlim(self._ecg_time - self.plot_window_seconds, self._ecg_time)

        self.lines["MAP (mmHg)"].set_data(list(self.time_data), list(self.map_data))
        self.axes["MAP (mmHg)"].set_xlim(current_time - self.plot_window_seconds, current_time)

        self.lines["SpO₂ (%)"].set_data(list(self.time_data), [val * 2.2 for val in self.map_data])
        self.axes["SpO₂ (%)"].set_xlim(current_time - self.plot_window_seconds, current_time)

        self.lines["EtCO₂ (mmHg)"].set_data(list(self.etco2_waveform_x_data), list(self.etco2_waveform_y_data))
        self.axes["EtCO₂ (mmHg)"].set_xlim(current_time - self.plot_window_seconds, current_time)


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
    twin = DigitalTwinModel(patient_id="demo_patient_gui", param_file="healthyFlat.json", sleep=True, time_step=0.02)
    if not hasattr(twin, 'redisClient') or twin.redisClient is None:
        twin.redisClient = DummyRedisClient()
    app = MonitorApp(twin)
    app.mainloop()

if __name__ == "__main__":
    main()