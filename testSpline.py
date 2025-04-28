# test script to generate splines

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import json

# Load your JSON (I'm pretending you already loaded it into a variable called `parameters`)
# Replace the current incorrect json.loads line with:
with open("parms_ranges.json", "r") as f:
    parameters = json.load(f)

# Dictionary to store spline functions
splines = {}

# Loop through all categories and parameters
for category, param_list in parameters.items():
    for param in param_list:
        name = f"{category}.{param['parameter']}"
        
        # Define x and y points
        x = np.array([-100, 0, 100])
        y = np.array([param['min'], param['mean'], param['max']])

        # Create the spline
        spline = CubicSpline(x, y)
        
        # Store in dictionary
        splines[name] = spline

# Now you can use splines like this:
example_param = 'cardio_parameters.elastance'
x_value = -30
y_value = splines[example_param](x_value)

print(f"At x={x_value}, {example_param} value is {y_value}")

# (Optional) plot a few splines
plt.figure(figsize=(10, 6))
for key in list(splines.keys())[:2]:  # Plot first 5 parameters as example
    x_fine = np.linspace(-100, 100, 500)
    plt.plot(x_fine, splines[key](x_fine), label=key)

plt.legend()
plt.grid(True)
plt.title("Example Splines from Parameters")
plt.xlabel("x [-100 to 100]")
plt.ylabel("parameter value")
plt.show()
