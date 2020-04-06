import _pysdot as ot
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
