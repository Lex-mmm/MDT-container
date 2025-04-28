import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# Given points
x_points = np.array([-100, 0, 100])
y_points = np.array([2, 4, 8])

# Create a cubic spline interpolation
cs = CubicSpline(x_points, y_points, bc_type='natural')

# Generate x values for plotting
x_plot = np.linspace(-120, 120, 500)
y_plot = cs(x_plot)

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(x_plot, y_plot, label='Cubic Spline')
plt.scatter(x_points, y_points, color='red', zorder=5, label='Given points')
plt.title('Cubic Spline Interpolation')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()
