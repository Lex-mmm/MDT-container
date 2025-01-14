import json
import numpy as np
class Alarms:
    def __init__(self, patient_id):
        self.patient_id = patient_id
        self.alarms = []
        with open("Alarms/alarm_setting.json", "r") as f:
            self.alarm_data = json.load(f)

    def check_alarms(self, data, timestamp):
        alarms = []
        msg = None
        for key, value in data.items():
            match value:
                case _ if value > self.alarm_data[key]["max"]:
                    msg = f"Alarm (high): {key} at {timestamp} s."
                    break
                
                case _ if value < self.alarm_data[key]["min"]:
                    msg = f"Alarm (low): {key} at {timestamp} s."
                    break
        if msg:
            alarms.append(msg)
            return alarms
        else: return None
    


