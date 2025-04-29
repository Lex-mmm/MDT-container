import numpy as np
import json
import pandas as pd
from scipy.interpolate import CubicSpline, interp1d
from matplotlib import pyplot as plt
import random
import statistics

class Pathology:
    def __init__(self):
        with open("healthyFlat.json", "r") as f:
            splineParameters = json.load(f)
        self.splineParameters = splineParameters    
        self.splineFunctions = {}
        for parameter in splineParameters:
            self.splineFunctions[parameter] = self.splineSolver(splineParameters[parameter])
            ## Spline functionality returns the value for the given parameter
        
        with open("Events/Pathologies/diseases.json", "r") as f:
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

    def CalcParamPercent(self, paramValue, paramName):
        curr_val = paramValue

        paramContent = self.splineParameters[paramName]
        paramMin = paramContent["min"]
        paramMax = paramContent["max"]
        paramMean = paramContent["value"]

        if paramMin == paramMean:
            searchMin = 0
            searchMax = 100
        elif paramMax == paramMean:
            searchMin = -100
            searchMax = 0
        else:
            searchMin = -100
            searchMax = 100

        print(f"min: {paramMin}, max: {paramMax}, value: {curr_val}")
        ## get the function for the parameter
        paramFuncton = self.splineFunctions[paramName]
        # create theoretical values for the parameter
        x_values = np.linspace(searchMin, searchMax, int(((searchMax - searchMin) / 0.1)+1) ) ## granularity
        ## get the values for the parameter
        y_values = paramFuncton(x_values)

        # Find where the spline crosses the target Y value
        indices = np.where(np.isclose(y_values, curr_val, atol=0.01))[0]

        ## get the indices of the values that are close to the current value
        listVal = x_values[indices].tolist()

        return statistics.mean(listVal) ## average percentage: No die-hard granularity

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

        # Get the event from the diseases.json file
        disease = self.diseases[eventName]
        if disease:
            diseaseParam = disease[f"severity_{eventSeverity}"]
            diseaseStartingCondition = diseaseParam["startingConditions"]
            diseaseDecayCondition = diseaseParam["decay"]

            # Process event for starting conditions
            startShell = {
                "event": eventName,
                "eventSeverity": eventSeverity,
                "eventType": "common",
                "timeCategorical": 'limited',
                "lastEmission": 0,
                "timeInterval": 0,
                "timeUnit": 'seconds',
                "eventCount": 1,  # Initial use only
                "parameters": []
            }

            for params in diseaseStartingCondition["parameters"]:
                paramFillShell = {
                    "name": params,
                    "value": diseaseStartingCondition["parameters"][params],
                    "action": "set",
                    "type": "relative"
                }
                startShell["parameters"].append(paramFillShell)

            print(f"Starting shell is: {startShell}")

            # Process event for decay conditions
            decayShell = {
                "event": eventName,
                "eventSeverity": eventSeverity,
                "eventType": "common",
                "timeCategorical": diseaseDecayCondition['rate']['type'],
                "lastEmission": None,
                "timeInterval": diseaseDecayCondition['rate']['timeValue'],
                "timeUnit": diseaseDecayCondition['rate']['timeUnit'],
                "eventCount": -10,  # Continuous function
                "parameters": []
            }

            for params in diseaseDecayCondition["parameters"]:
                paramData = diseaseDecayCondition["parameters"][params]
                paramFillShell = {
                    "name": params,
                    "value": round(random.uniform(paramData['min'], paramData['max']), 2),  # Two-decimal granularity
                    "action": "decay",
                    "type": "relative"
                }
                print(f"Adding shell based on {params}")
                decayShell["parameters"].append(paramFillShell)

            print(f"Decay shell is: {decayShell}")

            # Combine the processed events
            processedEvent = [startShell, decayShell]
            return processedEvent
        else:
            print(f"Event {eventName} not found in diseases.json")
            return None
