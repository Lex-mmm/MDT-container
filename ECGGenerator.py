from datetime import datetime
import numpy as np

class ECGGenerator:
    def __init__(self, samples_per_beat=500):
        self.samples_per_beat = samples_per_beat
        self.base_ecg = self._generate_base_waveform()
        self.last_time = 0.0
        self.buffer = np.array([])

    def _generate_base_waveform(self):
        t = np.linspace(0, 1, self.samples_per_beat)
        ecg = np.zeros_like(t)

        # Define wave parameters
        p_start, p_end = 0.1, 0.2
        qrs_start, qrs_end = 0.3, 0.35
        t_start, t_end = 0.5, 0.7

        p_amp, qrs_amp, t_amp = 0.15, 1.2, 0.35

        p_mask = (t >= p_start) & (t < p_end)
        ecg[p_mask] += p_amp * np.sin(np.pi * (t[p_mask] - p_start) / (p_end - p_start))

        qrs_mask = (t >= qrs_start) & (t < qrs_end)
        ecg[qrs_mask] += qrs_amp * np.sin(np.pi * (t[qrs_mask] - qrs_start) / (qrs_end - qrs_start))

        t_mask = (t >= t_start) & (t < t_end)
        ecg[t_mask] += t_amp * np.sin(np.pi * (t[t_mask] - t_start) / (t_end - t_start))

        return ecg

    def get_waveform(self, start_time, end_time, hr):
        if hr is None or hr < 30:
            hr = 60
        beat_duration = 60.0 / hr
        sample_interval = beat_duration / self.samples_per_beat

        num_samples = int((end_time - start_time) / sample_interval)
        times = np.linspace(start_time, end_time, num_samples)
        phases = ((times - self.last_time) % beat_duration) / beat_duration
        indices = (phases * self.samples_per_beat).astype(int) % self.samples_per_beat
        waveform = self.base_ecg[indices]
        self.last_time = end_time
        return times, waveform

