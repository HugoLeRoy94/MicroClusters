import pyvista as pv
import numpy as np

# Create grid data for the plot
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
x, y = np.meshgrid(x, y)
z = np.sin(np.sqrt(x**2 + y**2))

# Create a structured grid
grid = pv.StructuredGrid(x, y, z)

# Plot the grid
plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(grid, show_edges=True)

# Export to an interactive HTML file
plotter.export_html("pyvista_3d_plot.html")

print("3D plot saved as pyvista_3d_plot.html. You can open it in a web browser.")

