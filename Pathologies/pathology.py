import numpy as np
import json
import pandas as pd
from scipy.interpolate import CubicSpline


class Pathology:
    def __init__(self):
        with open("parameterSplines.json", "r") as f:
            splineParameters = json.load(f)
        
        self.splineFunctions = {}
        for parameter in splineParameters:
            self.splineFunctions[parameter] = self.splineSolver(splineParameters[parameter])
            ## Spline functionality returns the value for the given parameter
        
        with open("diseases.json", "r") as f:
            self.diseases = json.load(f)
    
    def splineSolver(self, parameterValues):
        ## Create a spline function for the parameter with min, max and mean value | Normal spline 
        ## Values range between -100 and 100 (percentual)
        x_points = np.array([parameterValues["min"], parameterValues["mean"], parameterValues["max"]])
        y_points = np.array([-100,0,100])

        cs = CubicSpline(x_points, y_points, bc_type='natural')
        return cs



    def solveEvent(self, event, eventSeverity, masterParameters):
        disease = self.diseases[event]
        if disease:
            diseaseParam = disease[f"{event}_{eventSeverity}"]

            for paramName, paramValue in diseaseParam.items():
                ## Get the spline function for the parameter
                splineFunction = self.splineFunctions[paramName]
                ## Get the value of the parameter
                diseaseParam[paramName] = splineFunction(paramValue)
            
            ## change parameters in the master parameter file
            for paramName, paramValue in diseaseParam.items():
                if paramName in masterParameters:
                    masterParameters[paramName] = paramValue
                else:
                    print(f"Parameter {paramName} not found in master parameters")

        
        
        