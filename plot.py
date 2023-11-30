import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go

# Assuming 'cluster' is an object with 'charges' and 'positions' attributes
charges = cluster.charges
x, y, z = cluster.positions[:, 0], cluster.positions[:, 1], cluster.positions[:, 2]

# Create an interactive 3D plot using Plotly
fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color=charges, colorscale='Viridis'))])

# Set labels for each axis
fig.update_layout(scene=dict(xaxis_title='x', yaxis_title='y', zaxis_title='z'))

# Display the interactive Plotly plot
fig.show()


# Create an interactive 3D plot using Plotly
fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color=charges, colorscale='Viridis'))])

# Set labels for each axis
fig.update_layout(scene=dict(xaxis_title='x', yaxis_title='y', zaxis_title='z'))

# Display the interactive Plotly plot
fig.show()
