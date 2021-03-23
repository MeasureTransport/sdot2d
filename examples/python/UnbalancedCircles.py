import pysdot as ot
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

"""
# Description:
Computes

"""


xbnds = [0.0,1.0] # minimum and maximum x values
ybnds = [0.0,1.0] # minimum and maximum y values
Ns = [50,50]
circle_radius = 0.1

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

    dens /= (grid.dx*grid.dy*np.sum(dens))

    return dens

class InitialState:
    """ Creates the initial density and optimal quantization of that density
        for a circle at position loc and radius.
    """
    def __init__(self, loc, radius, num_points):

        self.center = loc
        self.dens_vals = CreateCircleDensity(loc, radius)
        self.dist = ot.DiscretizedDistribution(grid, self.dens_vals)

        self.probs = (1.0/num_points)*np.ones(num_points)

        # Create an (balanced) optimal quantization of the initial circle
        opts = {'Lloyd Steps':200, 'Lloyd Tol':1e-5, 'GTol Abs':1e-8, 'Max Steps': 300}
        seed_pts = self._GenerateSeedPts(num_points)
        self.diag = ot.SemidiscreteOT.BuildCentroidal(self.dist, seed_pts, self.probs, opts)
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

                # Get the indices of the point
                xind = grid.LeftNode(pt[0])
                yind = grid.BottomNode(pt[1])

                # Check to see if the point is in a region of positive density
                if(self.dens_vals[xind,yind]>1e-3):
                    pts[:,i] = pt
                    in_supp = True

        return pts

def BalancedCase(initial_state):

    # Set the second position
    loc = np.array([0.25, 0.7])
    dens_vals = CreateCircleDensity(loc, circle_radius)
    dist = ot.DiscretizedDistribution(grid, dens_vals)

    # Compute the transport between the first points and the
    solver = ot.QuadraticRegularizedSDOT(dist, initial_state.pts, initial_state.probs, 1.0)

    num_pts = initial_state.pts.shape[1]
    opts = {'Max Steps':500, 'GTol Abs':1e-8, 'FTol Abs':0.0, 'Accept Ratio':0.01, 'Shrink Ratio':0.01}

    prices = np.ones(num_pts)
    prices, obj = solver.Solve(prices, opts)

    # Estimate the velocity
    diag = solver.Diagram()
    true_velocity = loc-initial_state.center
    ot_velocity = solver.MarginalCentroids()-initial_state.pts
    print('From:')
    print(initial_state.pts)
    print('To:')
    print(solver.MarginalCentroids())
    print('True Velocity:\n',true_velocity)
    print('OT Velocity:\n',ot_velocity)


    # Use the results to
    fig, ax = plt.subplots()
    ot.PlotDiagram(solver.Diagram(), ax, distribution=dist)

    plt.show()

def IntensityChange(initial_state):

    # Set the second position
    loc = np.array([0.25, 0.7])
    dens_vals = CreateCircleDensity(loc, circle_radius)
    dens_vals *= 0.7
    dist = ot.DiscretizedDistribution(grid, dens_vals)


    # Compute the transport between the first points and the
    solver = ot.QuadraticRegularizedSDOT(dist, initial_state.pts, initial_state.probs, 1.0)

    num_pts = initial_state.pts.shape[1]
    opts = {'Max Steps':500, 'GTol Abs':1e-8, 'FTol Abs':0.0, 'Accept Ratio':0.01, 'Shrink Ratio':0.01}

    prices = np.ones(num_pts)
    prices, obj = solver.Solve(prices, opts)

    # Estimate the velocity
    diag = solver.Diagram()
    true_velocity = loc-initial_state.center
    ot_velocity = solver.MarginalCentroids()-initial_state.pts
    print('From:')
    print(initial_state.pts)
    print('To:')
    print(solver.MarginalCentroids())
    print('True Velocity:\n',true_velocity)
    print('OT Velocity:\n',ot_velocity)

    # Use the results to
    fig, ax = plt.subplots()
    ot.PlotDiagram(solver.Diagram(), ax, distribution=dist)

    plt.show()


def PartiallyLeftDomain(initial_state):

    # Set the second position
    loc = np.array([0.75*circle_radius, 0.7])
    dens_vals = CreateCircleDensity(loc, circle_radius)

    # Scale the density so the values are the same as the initial state
    nnz0 = np.sum(initial_state.dens_vals>0)
    nnz1 = np.sum(dens_vals>0)
    dens_vals *= (nnz1/nnz0)

    dist = ot.DiscretizedDistribution(grid, dens_vals)

    # Compute the transport between the first points and the
    solver = ot.QuadraticRegularizedSDOT(dist, initial_state.pts, initial_state.probs, 5)

    num_pts = initial_state.pts.shape[1]
    opts = {'Max Steps':500, 'GTol Abs':1e-8, 'FTol Abs':0.0, 'Accept Ratio':0.01, 'Shrink Ratio':0.01}

    prices = np.ones(num_pts)
    prices, obj = solver.Solve(prices, opts)

    # Estimate the velocity
    diag = solver.Diagram()
    true_velocity = loc-initial_state.center
    ot_velocity = solver.MarginalCentroids()-initial_state.pts
    print('From:')
    print(initial_state.pts)
    print('To:')
    print(solver.MarginalCentroids())
    print('True Velocity:\n',true_velocity)
    print('OT Velocity:\n',ot_velocity)

    # Use the results to
    fig, ax = plt.subplots()
    ot.PlotDiagram(solver.Diagram(), ax, distribution=dist)

    plt.show()


loc0 = np.array([0.6,0.3])
initial_state = InitialState(loc0, circle_radius, 3)
print('\n\n================================')
print('Balanced Mass with Unbalanced Solver')
print('----------------------------------\n')
BalancedCase(initial_state)

print('\n\n================================')
print('Intensity Change')
print('----------------------------------\n')
IntensityChange(initial_state)

print('\n\n================================')
print('Partially Occluded')
print('----------------------------------\n')
PartiallyLeftDomain(initial_state)

print('\n\n')
