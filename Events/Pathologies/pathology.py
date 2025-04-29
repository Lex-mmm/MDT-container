import numpy as np
import json
import pandas as pd
from scipy.interpolate import CubicSpline, interp1d
from matplotlib import pyplot as plt


class Pathology:
    def __init__(self):
        with open("healthyFlat.json", "r") as f:
            splineParameters = json.load(f)
        
        self.splineFunctions = {}
        for parameter in splineParameters:
            self.splineFunctions[parameter] = self.splineSolver(splineParameters[parameter])
            ## Spline functionality returns the value for the given parameter
        
        with open("Pathologies/diseases.json", "r") as f:
            self.diseases = json.load(f)
    
    def splineSolver(self, parameterValues):
        ## Create a spline function for the parameter with min, max and mean value | Normal spline 
        ## Values range between -100 and 100 (percentual)
        min = parameterValues["min"]
        max = parameterValues["max"]
        mean = parameterValues["value"]
        if isinstance(mean, str):
            return None
        elif min == max == mean: ## no increasing sequence
            return None
        else:
            if min == mean:
                y_points = np.array([mean, max])
                x_points = np.array([0, 100])
                type = "Linear"
            elif max == mean:
                y_points = np.array([min, mean])
                x_points = np.array([-100, 0])
                type = "Linear"
            else:
                y_points = np.array([min, mean, max])
                x_points = np.array([-100, 0, 100])
                type = "Cubic"
            
            if type == "Linear":
                ## Create a linear spline function
                #print(f"Linear spline function for {parameterValues}")
                #print(f"x_points: {x_points}, y_points: {y_points}")
                cs = interp1d(x_points, y_points, kind='linear', fill_value="extrapolate")
            elif type == "Cubic":
                cs = CubicSpline(x_points, y_points, bc_type='natural')
            return cs

    def solveEvent(self, event, eventSeverity, masterParameters):
        ''' DEPRICATED: Use processPathology instead'''
        disease = self.diseases[event]
        if disease:
            diseaseParam = disease[f"severity_{eventSeverity}"]

            for paramName, paramValue in diseaseParam.items():
                ## Get the spline function for the parameter

                splineFunction = self.splineFunctions[paramName]
                #self.plot_spline(splineFunction, 0, 100)
                ## Get the value of the parameter
                val = splineFunction(paramValue)
                masterParameters[paramName]["value"] = val
            
            print(f"Novel setting for {paramName} = {masterParameters[paramName]['value']}")
            return masterParameters

            ## change parameters in the master parameter file

    def CalcParamPercent(self, paramContent, paramName):

        paramValue = paramContent["value"]
        paramMin = paramContent["min"]
        paramMax = paramContent["max"]

        paramFuncton = self.splineFunctions[paramName]
        ## get the
        x_values = np.linspace(paramMin, paramMax, 0.01)
        y_values = paramFuncton(x_values)

        # Find where the spline crosses the target Y value
        indices = np.where(np.isclose(y_values, paramValue, atol=1e-3))[0]
        listVal = x_values[indices].tolist()
        return listVal.mean()
        ## return the mean value of the list -> Average percentage for the current parameter, if multiple options

    
    def processPathology(self, eventName, eventSeverity):
        '''
        Input: 
            event: string, name of the event
            eventSeverity: int, severity of the event
        Output:
            processedEvent: events (n >=1), with following structure:
            {
                "event": eventName,                     String with name (e.g. myocardialInfarction)
                "eventSeverity": eventSeverity,         0 || 1 || 2 || 3 etc...
                "eventType": eventType,                 Disease || Therapy
                "timeCategorical": timeCategorical,     Continuous || Limited
                "lastEmission": lastEmission,           Last time of event emission (processed time), standard 0
                "timeInterval": timeInterval,           Time interval of the event (e.g. 0.5) for next processing
                "timeUnit": timeUnit,                   Time unit of the event (e.g. hours)
                "eventCount": eventCount,               Number of events, only used in case of limited timing: for single use, use n=1
                "parameters": {
                    name: {
                        value: [ Int || float ],
                        action: "set" || "decay"
                        type: "absolute" || "relative
                    },
                    .... 
            }
        '''
        ## Get the event from the diseases.json file
        if eventName not in self.diseases[eventName]:
            print(f"Event {eventName} not found in diseases.json")
            return None

        disease = self.diseases[eventName]
        if disease:

            diseaseParam = disease[f"severity_{eventSeverity}"]
            diseaseStartingCondition = diseaseParam["startingCondition"]
            diseaseDecayCondition = diseaseParam["decay"]

            ## process event for starting Condition:
            eventShell = {
                f"{eventName}_Start":{
                    "event": eventName,
                    "eventSeverity": eventSeverity,
                    "eventType": "disease",
                    "timeCategorical": 'limited',
                    "lastEmission": 0,
                    "timeInterval": 0,
                    "timeUnit": 's',
                    "eventCount": 1,
                    "parameters": {}
                },
                f"{eventName}_Decay":{
                    "event": eventName,
                    "eventSeverity": eventSeverity,
                    "eventType": "disease",
                    "timeCategorical": 'continuous',
                    "lastEmission": 0,
                    "timeInterval": diseaseDecayCondition[0][],
                    "timeUnit": 's',
                    "eventCount": 1,
                    "parameters": {}
                }
            }

