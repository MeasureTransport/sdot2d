import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt

opts = {'Max Steps': 500, 'GTol Abs':1e-1, 'XTol Abs':1e-2}


xbnds = [0.0,1.0] # minimum and maximum x values
ybnds = [0.0,1.0] # minimum and maximum y values
Ns = [40,40] # Number of cells in x and y direction

dx = (xbnds[1]-xbnds[0])/Ns[0]
dy = (ybnds[1]-ybnds[0])/Ns[1]
xs = 0.5*dx + dx*np.arange(Ns[0])
ys = 0.5*dy + dy*np.arange(Ns[1])

X,Y = np.meshgrid(xs,ys)

grid = ot.RegularGrid(xbnds[0],ybnds[0],xbnds[1],ybnds[1], Ns[0], Ns[1])

# Set the unnormalized density values for the "continuous" distribution
densVals = 0.0*np.ones(Ns)
radius = 0.1
densVals[(X-0.55)**2 + (Y-0.25)**2<radius**2] = 1.0
densVals[(X-0.55)**2 + (Y-0.75)**2<radius**2] = 1.0

tgtDens = ot.DiscretizedDistribution(grid, densVals)

# Create two points, one in each circle
numPts = 2

# Prescribe uniform weights
probs = (1.0/numPts)*np.ones(numPts)

x1s = np.linspace(0.001,0.7,500)
vals = np.zeros(x1s.shape)
derivs = np.zeros(x1s.shape)

for i in range(x1s.shape[0]):

    pts = np.array([[x1s[i],0.75],
                    [0.5, 0.5]])

    solver = ot.SemidiscreteOT(tgtDens, pts, probs)
    optPrices, vals[i] = solver.Solve(np.ones(numPts), opts)
    derivs[i] = solver.PointGradient()[0,0]

wassDists = -vals
# Plot the optimal Laguerre diagram
fig, ax = plt.subplots()
plt.plot(x1s,wassDists, label='$W_2$',linewidth=2)
plt.plot(x1s,derivs, label='Analytic',linewidth=2)
plt.plot(x1s[1:-1], (wassDists[2:]-wassDists[0:-2])/(x1s[2]-x1s[0]), label='Finite Difference',linewidth=1)
plt.legend()
plt.xlabel('X-Position of Point 1')

plt.show()
