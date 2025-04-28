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
            print(f"Parameter {parameter} with values: {splineParameters[parameter]}")

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
        disease = self.diseases[event]
        if disease:
            diseaseParam = disease[f"severity_{eventSeverity}"]
            print(diseaseParam)
            for paramName, paramValue in diseaseParam.items():
                ## Get the spline function for the parameter
                print(paramName, paramValue)
                splineFunction = self.splineFunctions[paramName]
                #self.plot_spline(splineFunction, 0, 100)
                ## Get the value of the parameter
                val = splineFunction(paramValue)
                masterParameters[paramName]["value"] = val
            
            return masterParameters

            ## change parameters in the master parameter file

    def plot_spline(self, spline_function, x_min, x_max, num_points=100):
        """
        Plots a spline function over a specified range.

        Args:
            spline_function (callable): The spline function to plot.
            x_min (float): The minimum x value for the plot.
            x_max (float): The maximum x value for the plot.
            num_points (int): The number of points to use for plotting (default: 100).
        """
        # Generate x values
        x_values = np.linspace(x_min, x_max, num_points)
        # Evaluate the spline function for each x value
        y_values = spline_function(x_values)
        
        # Plot the spline
        plt.figure(figsize=(8, 6))
        plt.plot(x_values, y_values, label="Spline Function", color="blue")
        plt.title("Spline Function Plot")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(True)
        plt.legend()
        plt.show()

    # Example usage
    # Assuming `cs` is a spline function created using CubicSpline or interp1d
    # plot_spline(cs, x_min=0, x_max=10)