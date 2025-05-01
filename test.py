import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d
import json

#from Pathologies.pathology import Pathology

def test(data):
    path = Pathology()
    disease = "myocardialInfarction"
    severity = 1

    with open("healthyFlat.json", "r") as f:
        masterParameters = json.load(f)
        
    masterParameters = path.solveEvent(disease, severity, masterParameters)
    print(f"Novel setting for Fi_O2 = {masterParameters['gas_exchange_params.FI_O2']['value']}")

def spline():
    y_points = np.array([0.21, 1])
    x_points = np.array([0, 100])
    cs = interp1d(x_points, y_points, kind='linear', fill_value="extrapolate")

    x = np.linspace(0, 100, 100)
    y = cs(x)
    plt.plot(x, y)
    plt.title("Linear Spline")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.show()

import datetime
from datetime import timedelta
datetime1 = datetime.datetime(2023, 10, 1, 12, 0)
datetime2 = datetime.datetime(2023, 10, 1, 12, 5)
## Timedelta in seconds

time = 10
unit = "seconds"
## convert variable to integer-based timedelta
print(int(timedelta(**{unit: time}).total_seconds()))
#timedelta = datetime2 - datetime1
#print(timedelta.total_seconds())

import random
from Events.Pathologies.pathology import Pathology
pathology = Pathology()
testing = pathology.CalcParamPercent(3, "cardio_parameters.elastance.max.9")
print(f"Testing: {testing}")


{'event': 'myocardialInfarction', 
            'eventSeverity': '5', 
            'eventType': 'common', 'timeCategorical': 'limited', 'lastEmission': 0, 'timeInterval': 0, 'timeUnit': 'seconds', 
            'eventCount': 1, 
            'parameters': [{'name': 'gas_exchange_params.FI_O2', 'value': 0.45, 'action': 'set', 'type': 'absolute'}]}