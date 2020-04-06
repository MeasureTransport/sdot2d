import pysdot as ot
import numpy as np

from matplotlib.patches import Circle, Wedge, Polygon, Patch
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt


def PlotDiagram(lagDiag, ax):
    """ Plots the cells, cell centroids, and seeds points from a Laguerre diagram. """

    patches = []
    for cellInd in range(lagDiag.NumCells()):
        verts = lagDiag.GetCellVertices(cellInd)
        patches.append( Polygon(verts.T) )

    #colors = list(range(lagDiag.NumCells()))
    p = PatchCollection(patches, facecolor='gray',edgecolor='k', alpha=0.4)
    #p.set_array(np.array(colors))
    ax.add_collection(p)
    #fig.colorbar(p, ax=ax)

    centroids = lagDiag.Centroids()
    seeds = lagDiag.SeedPts()

    # Draw lines between each centroid and seed point
    for i in range(lagDiag.NumCells()):
        ax.plot([seeds[0,i], centroids[0,i]], [seeds[1,i], centroids[1,i]], 'r', linewidth=0.4)

    # Plot the cell centroids
    ax.plot(centroids[0,:],centroids[1,:],'.b')

    # Plot the seed points
    ax.plot(seeds[0,:], seeds[1,:],'.k')


    legend_elements = [Patch(facecolor='gray',edgecolor='k', alpha=0.4, label='Laguerre Cells'),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='b', markersize=10, label='Centroids'),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='k', markersize=10, label='Seeds')]

    ax.legend(handles=legend_elements, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=3, mode="expand", borderaxespad=0.)

    bbox = lagDiag.BoundBox()
    ax.set_xlim(bbox.xMin, bbox.xMax)
    ax.set_ylim(bbox.yMin, bbox.yMax)

    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')




numPts = 100

xbnds = [0.0,1.0] # minimum and maximum x values
ybnds = [0.0,1.0] # minimum and maximum y values
Ns = [10,10] # Number of cells in x and y direction

grid = ot.RegularGrid(xbnds[0],ybnds[0],xbnds[1],ybnds[1], Ns[0], Ns[1])

# Set the unnormalized density values for the "continuous" distribution
densVals = np.ones(Ns)
tgtDens = ot.DiscretizedDistribution(grid, densVals)

# Create random points in the domain of interest
pts = xbnds[0] + (xbnds[1]-xbnds[0])*np.random.rand(2,numPts)

# Prescribe uniform weights
probs = (1.0/numPts)*np.ones(numPts)

solver = ot.SemidiscreteOT(tgtDens, pts, probs)

initialPrices = np.ones(numPts)
optPrices, optObj = solver.Solve(initialPrices)

# Plot the optimal Laguerre diagram
fig, ax = plt.subplots()
PlotDiagram(solver.Diagram(), ax)
plt.savefig('RandomRectangle_LaguerreDiagram.pdf')
