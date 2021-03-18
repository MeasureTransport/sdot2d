import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt

numPts = 100

xbnds = [0.0,1.0] # minimum and maximum x values
ybnds = [0.0,1.0] # minimum and maximum y values
Ns = [50,50]

bbox = ot.BoundingBox(xbnds[0],xbnds[1],ybnds[0],ybnds[1])

grid = ot.RegularGrid(bbox, Ns[0], Ns[1])

# Create the density, which is one inside a rectangular region and near zero outside
densVals = np.ones(Ns)
for i in range(Ns[0]):
    for j in range(Ns[1]):
        pt = grid.Center(i,j)
        if((pt[0]<0.8)&(pt[1]<0.6)&(pt[1]>0.4)):
            densVals[i,j] = 1.0
        else:
            densVals[i,j] = 1e-4



# Set the options for Lloyd's algorithm and the underlying optimization problem for the node prices.
opts = {'Lloyd Steps':200, 'Lloyd Tol':1e-6, 'Max Steps': 500}

# Construct the centroidal diagram with capacity constraints
dist = ot.DiscretizedDistribution(grid,densVals)
diag = ot.SemidiscreteOT.BuildCentroidal(dist, numPts, opts)

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
