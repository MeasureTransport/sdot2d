import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt


xbnds = [0.0,1.0] # minimum and maximum x values
ybnds = [0.0,1.0] # minimum and maximum y values
Ns = [50,50]

bbox = ot.BoundingBox(xbnds[0],xbnds[1],ybnds[0],ybnds[1])

grid = ot.RegularGrid(bbox, Ns[0], Ns[1])

# Create the "discretized distribution"
densVals = np.ones(Ns)
dist = ot.DiscretizedDistribution(grid,densVals)

# Generate random points
numPoints = 10
points = np.random.rand(2,numPoints)
prices = np.ones(numPoints)

# Set up a Laguerre diagram
diag = ot.LaguerreDiagram(xbnds[0], xbnds[1], ybnds[0], ybnds[1], points, prices)

# Visualize the results, using the area of each cell as a color
areas = diag.Areas(dist)
print('Minimum area = ', np.min(areas))
print('Maximum area = ', np.max(areas))

# Plot the resulting centroidal Voronoi diagram
fig, axs = plt.subplots(ncols=2,figsize=(14,6))
ot.PlotDiagram(diag, axs[0], distribution=dist, cell_colors=areas)
axs[0].set_title('Weighted CVD')

axs[1].imshow(densVals.T,extent=[xbnds[0],xbnds[1],ybnds[0],ybnds[1]],origin='lower',alpha=0.8)
axs[1].set_title('Density')

plt.show()
