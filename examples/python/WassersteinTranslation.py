"""
Plots the Wasserstein-2 distance as a function of translational distance.
"""


import pysdot as ot
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

xbnds = [0.0,1.0] # minimum and maximum x values
circle_radius = 0.1
ybnds = [0.5-2*circle_radius,0.5+2*circle_radius] # minimum and maximum y values
Ns = [150,int(150 * (ybnds[1]-ybnds[0])/(xbnds[1]-xbnds[0]))]

bbox = ot.BoundingBox(xbnds[0],xbnds[1],ybnds[0],ybnds[1])
grid = ot.RegularGrid(bbox, Ns[0], Ns[1])

def CreateCircleDensity(loc, radius):
    """ Returns a normalized density defined on the grid with a uniform circular
        region of nonzero density.

        ARGUMENTS:
            loc (np.array) : Center point of the circle.
            radius (float) : Radius of the circle.  Any grid cells with centroids inside this radius will be nonzero.
    """

    dens = np.zeros(Ns)
    for i in range(Ns[0]):
        for j in range(Ns[1]):
            pt = grid.Center(i,j)
            if( la.norm(pt-loc) < radius ):
                dens[i,j] = 1.0
            else:
                dens[i,j] = 0.0
    dens /= (np.pi*radius**2)
    return dens

class InitialState:
    """ Creates the initial density and optimal quantization of that density
        for a circle at position loc and radius.
    """
    def __init__(self, loc, radius, num_points, dens_vals0, ot_class,penalty=None):

        self.center = loc
        self.radius = radius
        self.dens_vals = dens_vals0

        self.dist = ot.DiscretizedDistribution(grid, self.dens_vals)

        self.probs = np.ones(num_points)*np.sum(grid.dx*grid.dy*dens_vals0)/num_points

        # Create an (balanced) optimal quantization of the initial circle
        opts = {'Lloyd Steps':200, 'Lloyd Tol':1e-5, 'GTol Abs':1e-9, 'Max Steps': 300}
        if (penalty is not None):
            opts['Penalty'] = penalty

        seed_pts = self._GenerateSeedPts(num_points)
        self.diag = ot_class.BuildCentroidal(self.dist, seed_pts, self.probs, opts)
        self.pts = self.diag.Centroids(self.dist)


    def _GenerateSeedPts(self, num_points):
        """ Generates random points in the domain that lie in the support (nonzero regions)
            of the density.  This is a good way to seed the centroidal solver.

            ARGUMENTS:
                num_point (int) : Number of points to generate.
                dens (np.array) : 2D array of density values used to define the support.

            RETURNS:
                np.array : a 2xN array of points
        """
        scale = np.array([xbnds[1]-xbnds[0], ybnds[1]-ybnds[0]])
        offset = np.array([xbnds[0], ybnds[0]])

        pts = np.zeros((2,num_points))
        for i in range(num_points):

            in_supp = False
            while(not in_supp):

                # Generate a random point in the domain
                pt = offset + scale * np.random.rand(2)

                if( la.norm(pt-self.center) < self.radius ):
                    in_supp = True
                    pts[:,i] = pt

        return pts





num_pts = 15
loc0 = np.array([circle_radius+0.05, 0.5])
dens_vals0 = CreateCircleDensity(loc0, circle_radius)

dxs = np.linspace(loc0[0],1+0.9*circle_radius,50)

vmin = -10
vmax = np.max(dens_vals0)

w2s = []
l2s = []
qrs = []
ghks = []


#######################
# Balanced
initial_state = InitialState(loc0, circle_radius, num_pts, dens_vals0, ot.SemidiscreteOT)
prices = initial_state.diag.Prices()

plt.figure(figsize=(10,10*(Ns[1]/Ns[0])))
#ot.PlotDiagram(initial_state.diag, distribution=initial_state.dist)
plt.imshow(dens_vals0.T, extent=[xbnds[0],xbnds[1],ybnds[0],ybnds[1]], cmap='Greys_r', vmin=vmin,vmax=vmax)
seed_pts = initial_state.diag.Centroids(initial_state.dist)
plt.plot(seed_pts[0,:],seed_pts[1,:],'.r', markersize=14)
plt.axis('off')
plt.savefig("InitialDensity.png", bbox_inches='tight')
plt.close()
plt.show()
quit()

i = 0
for dx in dxs:
    loc = np.array([dx, 0.5])

    dens_vals = CreateCircleDensity(loc, circle_radius)

    plt.figure(figsize=(10,10*(Ns[1]/Ns[0])))
    plt.imshow(dens_vals.T, cmap='Greys_r', vmin=vmin,vmax=vmax)
    plt.axis('off')
    plt.savefig("Density{:02d}.png".format(i), bbox_inches='tight')
    plt.close()

    l2s.append(0.01*np.sum((dens_vals-dens_vals0)**2)/(np.prod(Ns)))

    dens_vals *= (np.sum(dens_vals0)/np.sum(dens_vals))

    dist = ot.DiscretizedDistribution(grid, dens_vals)

    # Compute the transport between the first points and the
    solver = ot.SemidiscreteOT(dist, initial_state.pts, initial_state.probs)

    num_pts = initial_state.pts.shape[1]
    opts = {'Max Steps':500, 'GTol Abs':1e-6, 'FTol Abs':0.0, 'Accept Ratio':0.01, 'Shrink Ratio':0.01}

    prices, obj = solver.Solve(prices, opts)

    w2s.append(-2*obj)

    i+=1
    print('Wasserstein Distance: ', -obj)


#######################
# QR Unbalanced
penalty = 0.75
initial_state = InitialState(loc0, circle_radius, num_pts, dens_vals0, ot.SemidiscreteQR, penalty)
prices = initial_state.diag.Prices()


for dx in dxs:
    loc = np.array([dx, 0.5])

    dens_vals = CreateCircleDensity(loc, circle_radius)
    dist = ot.DiscretizedDistribution(grid, dens_vals)

    # Compute the transport between the first points and the
    solver = ot.SemidiscreteQR(dist, initial_state.pts, initial_state.probs,penalty)

    num_pts = initial_state.pts.shape[1]
    opts = {'Max Steps':500, 'GTol Abs':1e-6, 'FTol Abs':0.0, 'Accept Ratio':0.01, 'Shrink Ratio':0.01}

    prices, obj = solver.Solve(prices, opts)

    qrs.append(-2*obj)
    print('Wasserstein Distance: ', -obj)


#######################
# GHK Unbalanced
initial_state = InitialState(loc0, circle_radius, num_pts, dens_vals0, ot.SemidiscreteGHK,penalty)
prices = initial_state.diag.Prices()


for dx in dxs:
    loc = np.array([dx, 0.5])

    dens_vals = CreateCircleDensity(loc, circle_radius)
    dist = ot.DiscretizedDistribution(grid, dens_vals)

    # Compute the transport between the first points and the
    solver = ot.SemidiscreteGHK(dist, initial_state.pts, initial_state.probs, penalty)

    num_pts = initial_state.pts.shape[1]
    opts = {'Max Steps':500, 'GTol Abs':1e-6, 'FTol Abs':0.0, 'Accept Ratio':0.01, 'Shrink Ratio':0.01}

    prices, obj = solver.Solve(prices, opts)

    ghks.append(-2*obj)
    print('Wasserstein Distance: ', -obj)


plt.figure()
plt.plot(dxs,w2s, linewidth=2, label='W2')
plt.plot(dxs,qrs, linewidth=2, label='QR')
plt.plot(dxs,ghks, linewidth=2, label='GHK')
plt.plot(dxs, l2s, linewidth=2, label='L2')

plt.ylabel('Distance', fontsize=14)
plt.xlabel('Translation',fontsize=14)
plt.legend()

plt.savefig('MetricComparison.png')
plt.show()


for i in range(len(dxs)):

    plt.figure()
    plt.plot(dxs,w2s, linewidth=2, label='W2')
    plt.plot(dxs, l2s, linewidth=2, label='L2')

    ymin,ymax = plt.ylim()
    plt.plot([dxs[i],dxs[i]],[ymin,ymax], '--k')
    plt.ylim(ymin,ymax)

    plt.xlim(dxs[0],1.0-circle_radius)

    plt.ylabel('Distance', fontsize=14)
    plt.xlabel('Translation',fontsize=14)
    plt.legend(loc='upper left')

    plt.savefig('L2Comparison_{:02d}.png'.format(i),bbox_inches='tight')
    plt.close()


    plt.figure()
    plt.plot(dxs,w2s, linewidth=2, label='W2')
    plt.plot(dxs, l2s, linewidth=2, label='L2')
    plt.plot(dxs,qrs, linewidth=2, label='QR')
    plt.plot(dxs,ghks, linewidth=2, label='GHK')

    ymin,ymax = plt.ylim()
    plt.plot([dxs[i],dxs[i]],[ymin,ymax], '--k')
    plt.ylim(ymin,ymax)

    plt.ylabel('Distance', fontsize=14)
    plt.xlabel('Translation',fontsize=14)
    plt.legend(loc='upper left')

    plt.savefig('MetricComparison_{:02d}.png'.format(i),bbox_inches='tight')
    plt.close()
