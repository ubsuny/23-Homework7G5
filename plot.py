#use import matplotlib.pyplot as plt instead of matplotlib notebook to achieve interactive plot.
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go

# Create a 3D plot using Matplotlib
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Assuming 'cluster' is an object with 'charges' and 'positions' attributes
charges = cluster.charges
x, y, z = cluster.positions[:, 0], cluster.positions[:, 1], cluster.positions[:, 2]

# Scatter plot with color-coded charges using the 'coolwarm' colormap
ax.scatter(x, y, z, c=charges, cmap='coolwarm')

# Set labels for each axis
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# Display the static Matplotlib plot
plt.show()

# Create an interactive 3D plot using Plotly
fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color=charges, colorscale='coolwarm'))])

# Set labels for each axis
fig.update_layout(scene=dict(xaxis_title='x', yaxis_title='y', zaxis_title='z'))

# Display the interactive Plotly plot
fig.show()

