import json
import numpy as np

class Inference:

    def __init__(self):
         
        with open("Inference/inference_setting.json", "r") as f:
            self.inference_data = json.load(f)
            self.status = {
            "heart_rate": str(self.inference_data["start"]["heart_rate"]),
            "RR": str(self.inference_data["start"]["RR"]),
            "SaO2": str(self.inference_data["start"]["SaO2"]),
            "MAP": str(self.inference_data["start"]["MAP"]),
            "etCO2": str(self.inference_data["start"]["etCO2"])
        }
        

    def check_inference(self, data):
        for key, value in data.items():
            value_swap = self.inference_data[key][ self.status[key] ]
            norm = np.random.normal(0, float(value_swap), 1)
            sd_alt = np.round(norm, 2)
            data[key] = float(value) + float(sd_alt)
            if key == "SpO2" and data["SpO2"] > 100:
                data["SpO2"] = 100
                ## manually limit saturation curve
        return data

    def change_inference(self, param, new_value):
        new_value = str(new_value)        
        self.status[param] = new_value
        return f"Inference level for {param} changed to {new_value}"


