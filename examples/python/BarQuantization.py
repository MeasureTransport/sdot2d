import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt

numPts = 100

xbnds = [0.0,1.0] # minimum and maximum x values
ybnds = [0.0,1.0] # minimum and maximum y values
Ns = [50,50]

bbox = ot.BoundingBox(xbnds[0],xbnds[1],ybnds[0],ybnds[1])

grid = ot.RegularGrid(bbox, Ns[0], Ns[1])
densVals = np.ones(Ns)
for i in range(Ns[0]):
    for j in range(Ns[1]):
        pt = grid.Center(i,j)
        if((pt[0]<0.8)&(pt[1]<0.6)&(pt[1]>0.4)):
            densVals[i,j] = 1.0
        else:
            densVals[i,j] = 1e-4

dist = ot.DiscretizedDistribution(grid,densVals)
diag = ot.SemidiscreteOT.BuildCentroidal(dist, numPts)

areas = diag.Areas(dist)

# Plot the resulting centroidal Voronoi diagram
fig, axs = plt.subplots(ncols=2,figsize=(14,6))
ot.PlotDiagram(diag, axs[0], distribution=dist, cell_colors=areas)
axs[0].set_title('Weighted CVD')

axs[1].imshow(densVals.T,extent=[xbnds[0],xbnds[1],ybnds[0],ybnds[1]],origin='lower',alpha=0.8)
axs[1].set_title('Density')

plt.show()
