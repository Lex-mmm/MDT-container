from flask import Flask, request, jsonify
import os
import time
import requests
import numpy as np
from scipy.integrate import solve_ivp
from datetime import datetime, timedelta
import json

# Configuration class to handle JSON-based settings
class Config:
    properties_path: str | None = None

    @staticmethod
    def selectPropertyFile(path: str) -> None:
        Config.properties_path = path

    @staticmethod
    def get(path: str):
        with open(Config.properties_path, 'r') as file:
            config = json.load(file)
            vars = path.split('/')
            for var in vars:
                config = config[var]
            return config

# Random generator for identifiers and dates
class RandomGenerator:
    def __init__(self, seed: str) -> None:
        np.random.seed(seed)

    @staticmethod
    def generateRandomIdentifier(length: int) -> str:
        return ''.join(np.random.choice(list("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"), length))

    @staticmethod
    def generateRandomDate(start: datetime, end: datetime) -> datetime:
        return start + timedelta(
            seconds=np.random.randint(0, int((end - start).total_seconds()))
        )

# Patient setup
class Patient:
    def __init__(self) -> None:
        self.id: str = RandomGenerator.generateRandomIdentifier(7)
        self.bsn: int = RandomGenerator.generateRandomIdentifier(9)
        self.dob: datetime = datetime.strptime(Config.get('patient/dob'), "%Y-%m-%d")

# Flask app for communication
app = Flask(__name__)
Config.selectPropertyFile("patient_1.json")
RandomGenerator(seed=12345)

def update_heart_period(HR):
    """Calculate heart period parameters."""
    HP = 60 / HR
    Tas = 0.03 + 0.09 * HP
    Tav = 0.01
    Tvs = 0.16 + 0.2 * HP
    return HR, HP, Tas, Tav, Tvs

# Initialize heart period parameters
HR, HP, Tas, Tav, Tvs = update_heart_period(Config.get("cardio_control_params/HR_n"))

# Define core ODE system
def extended_state_space_equations(t, x):
    # Define parameters
    params = {
        "HR_n": Config.get("cardio_control_params/HR_n"),
        "R_n": Config.get("cardio_control_params/R_n"),
        "UV_n": Config.get("cardio_control_params/UV_n"),
        "FI_O2": Config.get("gas_exchange_params/FI_O2"),
        "FI_CO2": Config.get("gas_exchange_params/FI_CO2")
    }

    # Define the states for the cardiovascular and respiratory model
    V = x[:10]  # Cardiovascular volumes
    mechanical_states = x[10:15]
    FD_O2, FD_CO2, p_a_CO2, p_a_O2 = x[15], x[16], x[17], x[18]
    c_Stis_CO2, c_Scap_CO2, c_Stis_O2, c_Scap_O2 = x[19:23]
    Delta_RR_c = x[23]
    Delta_Pmus_c = x[24]
    Pmus = x[25]
    Delta_HR_c = x[26]
    Delta_R_c = x[27]
    Delta_UV_c = x[28]

    # Inputs and pressures
    ela, elv, era, erv = get_inputs(t, params["HR_n"])
    HR = params["HR_n"] - Delta_HR_c
    R_c = params["R_n"] - Delta_R_c
    UV_c = params["UV_n"] + Delta_UV_c

    # Calculate pressures and flows
    # Placeholders for calculations, e.g., P, F, dVdt based on cardiovascular model
    # Similar calculations for respiratory system omitted here for brevity
    dVdt = np.zeros(10)
    dFD_O2_dt, dFD_CO2_dt, dp_a_CO2, dp_a_O2 = 0, 0, 0, 0  # Replace with actual formulas

    # Mechanical state derivatives and cardiac control derivatives
    dxdt_mechanical = np.zeros(5)
    dDelta_Pmus_c, dDelta_RR_c, dDelta_HR_c, dDelta_R_c, dDelta_UV_c = 0, 0, 0, 0, 0  # Control system

    # Concatenate all derivatives
    dxdt = np.concatenate([dVdt, dxdt_mechanical, [dFD_O2_dt, dFD_CO2_dt, dp_a_CO2, dp_a_O2,
                                                  dDelta_RR_c, dDelta_Pmus_c, Pmus,
                                                  dDelta_HR_c, dDelta_R_c, dDelta_UV_c]])
    return dxdt

# Helper function for ventilator pressure control
def ventilator_pressure(t):
    RR = Config.get("respiratory_control_params/RR_0")
    PEEP = 5
    peak_pressure = 20
    T = 60 / RR
    cycle_time = t % T
    return peak_pressure if cycle_time < T / 2 else PEEP
def get_inputs(t, HR):
    """Get inputs for the cardiovascular system."""
    HR, HP, Tas, Tav, Tvs = update_heart_period(HR)
    ela, era, elv, erv = calculate_elastances(t, HR,Tas, Tav, Tvs, 0.01, elastance)
    return np.array([ela, elv, era, erv])


def generate_sinus_ecg_segment(t, heart_rate):
    # Calculate the duration of one cardiac cycle based on the heart rate
    cycle_duration = 60.0 / heart_rate
    segment_duration = t[-1] - t[0]

    # Proportional durations of the ECG components based on cycle duration
    p_wave_duration = 0.1 * cycle_duration
    p_wave_amplitude = 0.2
    qrs_duration = 0.1 * cycle_duration
    qrs_amplitude = 1.5
    t_wave_duration = 0.2 * cycle_duration
    t_wave_amplitude = 0.5

    # Initialize the ECG segment with zeros
    ecg_wave = np.zeros_like(t)

    # Calculate the start times within the cycle for P wave, QRS complex, and T wave
    p_wave_start = 0.1 * cycle_duration
    qrs_start = 0.3 * cycle_duration
    t_wave_start = 0.5 * cycle_duration

    # Modulo the time with the cycle duration to get the position within the cycle
    cycle_position = t % cycle_duration

    # P wave
    mask_p = (cycle_position >= p_wave_start) & (cycle_position < p_wave_start + p_wave_duration)
    ecg_wave[mask_p] = p_wave_amplitude * np.sin(np.pi * (cycle_position[mask_p] - p_wave_start) / p_wave_duration)

    # QRS complex
    mask_qrs = (cycle_position >= qrs_start) & (cycle_position < qrs_start + qrs_duration)
    ecg_wave[mask_qrs] = qrs_amplitude * np.sin(np.pi * (cycle_position[mask_qrs] - qrs_start) / qrs_duration)

    # T wave
    mask_t = (cycle_position >= t_wave_start) & (cycle_position < t_wave_start + t_wave_duration)
    ecg_wave[mask_t] = t_wave_amplitude * np.sin(np.pi * (cycle_position[mask_t] - t_wave_start) / t_wave_duration)

    # Add noise
    noise = np.random.normal(0, 0.02, t.shape)
    ecg_wave += noise

    return ecg_wave


def generate_vt_ecg_segment(t, heart_rate):
    # Calculate the duration of one cardiac cycle based on the heart rate
    heart_rate = 110
    cycle_duration = 60.0 / heart_rate
    segment_duration = t[-1] - t[0]

    # VT is characterized by wide QRS complexes and absent P and T waves
    qrs_duration = 0.2 * cycle_duration  # Wider QRS complex
    qrs_amplitude = 1.0  # Lower amplitude compared to normal QRS

    # Initialize the ECG segment with zeros
    ecg_wave = np.zeros_like(t)

    # Calculate the start time within the cycle for QRS complex
    qrs_start = 0.2 * cycle_duration

    # Modulo the time with the cycle duration to get the position within the cycle
    cycle_position = t % cycle_duration

    # QRS complex (simulating ventricular tachycardia morphology)
    mask_qrs = (cycle_position >= qrs_start) & (cycle_position < qrs_start + qrs_duration)
    ecg_wave[mask_qrs] = qrs_amplitude * np.sign(np.sin(np.pi * (cycle_position[mask_qrs] - qrs_start) / qrs_duration))

    # Add irregularity to the waveform to simulate VT variability
    irregularity = 0.1 * np.sin(2 * np.pi * cycle_position / (cycle_duration * 2))
    ecg_wave += irregularity

    # Add noise
    noise = np.random.normal(0, 0.05, t.shape)
    ecg_wave += noise

    return ecg_wave




def calculate_respiratory_curve(Pintra_t0, RRP, t):
    """Calculate the respiratory curve."""
    return Pintra_t0 + (1 + np.cos(2 * np.pi / RRP * t))

def calculate_elastances(t, HR, Tas, Tav, Tvs, T, elastance):
    """Calculate heart elastances."""
    HP = 60/HR

    ncc = (t % HP) / T

    if ncc <= round(Tas / T):
        aaf = np.sin(np.pi * ncc / (Tas / T))
    else:
        aaf = 0

    ela = elastance[0, 8] + (elastance[1, 8] - elastance[0, 8]) * aaf
    era = elastance[0, 4] + (elastance[1, 4] - elastance[0, 4]) * aaf

    if ncc <= round((Tas + Tav) / T):
        vaf = 0
    elif ncc <= round((Tas + Tav + Tvs) / T):
        vaf = np.sin(np.pi * (ncc - (Tas + Tav) / T) / (Tvs / T))
    else:
        vaf = 0

    elv = elastance[0, 9] + (elastance[1, 9] - elastance[0, 9]) * vaf
    erv = elastance[0, 5] + (elastance[1, 5] - elastance[0, 5]) * vaf

    return ela, era, elv, erv




class PatientSimulation:
    def __init__(self, patient_id, pathology):
        self.patient = Patient()
        self.patient_id = patient_id
        self.pathology = pathology

        # Initial conditions array
        self.state = np.zeros(29)

    def run_simulation(self):
        while True:
            sol = solve_ivp(extended_state_space_equations, [0, 1], self.state)
            self.state = sol.y[:, -1]
            self.check_conditions()
            time.sleep(1)

    def check_conditions(self):
        if self.state[0] > 6:
            event_data = {"patient_id": self.patient_id, "event": "CO2 level high"}
            self.send_event_to_main(event_data)

    def send_event_to_main(self, event_data):
        try:
            response = requests.post("http://main_controller:5000/patient_event", json=event_data)
            print(f"Sent event to main controller: {response.status_code}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to send event: {e}")

@app.route("/event", methods=["POST"])
def receive_event():
    data = request.json
    event = data.get("event")
    print(f"Patient {os.getenv('PATIENT_ID')}: Received event - {event}")
    return jsonify({"status": "event processed"}), 200

if __name__ == "__main__":
    patient_id = os.getenv("PATIENT_ID")
    pathology = os.getenv("PATHOLOGY")
    simulation = PatientSimulation(patient_id, pathology)
    import threading
    threading.Thread(target=simulation.run_simulation).start()
    app.run(host="0.0.0.0", port=5000)
