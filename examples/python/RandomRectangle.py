import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt

numPts = 3

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
ot.PlotDiagram(solver.Diagram(), ax)

plt.savefig('RandomRectangle_LaguerreDiagram.pdf')
plt.show()
