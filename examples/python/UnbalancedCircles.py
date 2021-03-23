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


def GenerateSeedPts(num_points, dens):
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
            if(dens[xind,yind]>1e-3):
                pts[:,i] = pt
                in_supp = True

    return pts

def BalancedCase():

    # Set the initial circle position
    loc0 = np.array([0.6,0.3])
    dens_vals0 = CreateCircleDensity(loc0, circle_radius)
    dist0 = ot.DiscretizedDistribution(grid, dens_vals0)

    # Set the second position
    loc1 = np.array([0.25, 0.7])
    dens_vals1 = CreateCircleDensity(loc1, circle_radius)
    dist1 = ot.DiscretizedDistribution(grid, dens_vals1)

    # Generate some random initial points inside the initial circle
    numPts = 20
    probs = (1.0/numPts)*np.ones(numPts)

    # Create an (balanced) optimal quantization of the initial circle
    opts = {'Lloyd Steps':200, 'Lloyd Tol':1e-5, 'GTol Abs':1e-8, 'Max Steps': 300}
    seed_pts = GenerateSeedPts(numPts, dens_vals0)
    diag = ot.SemidiscreteOT.BuildCentroidal(dist0, seed_pts, probs, opts)
    pts = diag.Centroids(dist0)

    # Compute the transport between the first points and the
    solver = ot.SemidiscreteOT(dist1, pts, probs)

    opts = {'Max Steps':500, 'GTol Abs':1e-8, 'FTol Abs':0.0}
    h = loc0-loc1
    newPrices = np.ones(numPts)#diag.Prices() -  2*np.sum(np.tile(h.reshape((2,1)),(1,numPts))*pts,axis=0)# Theorem 2.3 of James Ronan's thesis
    optPrices, optObj = solver.Solve(newPrices, opts)

    # Estimate the velocity
    diag2 = solver.Diagram()
    true_velocity = loc1-loc0
    ot_velocity = diag2.Centroids(dist1)-pts
    print('From:')
    print(pts)
    print('To:')
    print(diag2.Centroids(dist1))
    print('True Velocity:\n',true_velocity)
    print('OT Velocity:\n',ot_velocity)


    # Use the results to
    fig, ax = plt.subplots()
    ot.PlotDiagram(solver.Diagram(), ax, distribution=dist1)

    plt.show()

def IntensityChange():

    # Set the initial circle position
    loc0 = np.array([0.6,0.3])
    dens_vals0 = CreateCircleDensity(loc0, circle_radius)
    dist0 = ot.DiscretizedDistribution(grid, dens_vals0)

    # Set the second position
    loc1 = np.array([0.25, 0.7])
    dens_vals1 = CreateCircleDensity(loc1, circle_radius)
    dens_vals1 *= 0.7
    dist1 = ot.DiscretizedDistribution(grid, dens_vals1)

    # Generate some random initial points inside the initial circle
    numPts = 4
    probs = (1.0/numPts)*np.ones(numPts)

    # Create an (balanced) optimal quantization of the initial circle
    opts = {'Lloyd Steps':200, 'Lloyd Tol':1e-5, 'GTol Abs':1e-8, 'Max Steps': 300}
    seed_pts = GenerateSeedPts(numPts, dens_vals0)
    diag = ot.SemidiscreteOT.BuildCentroidal(dist0, seed_pts, probs, opts)
    pts = diag.Centroids(dist0)

    # Compute the transport between the first points and the
    solver = ot.QuadraticRegularizedSDOT(dist1, pts, probs, 0.5)

    opts = {'Max Steps':500, 'GTol Abs':3e-2, 'FTol Abs':0.0}
    h = loc0-loc1
    newPrices = np.ones(numPts)#diag.Prices() -  2*np.sum(np.tile(h.reshape((2,1)),(1,numPts))*pts,axis=0)# Theorem 2.3 of James Ronan's thesis
    optPrices, optObj = solver.Solve(newPrices, opts)

    # Estimate the velocity
    diag2 = solver.Diagram()
    true_velocity = loc1-loc0
    ot_velocity = diag2.Centroids(dist1)-pts
    print('From:')
    print(pts)
    print('To:')
    print(diag2.Centroids(dist1))
    print('True Velocity:\n',true_velocity)
    print('OT Velocity:\n',ot_velocity)


    # Use the results to
    fig, ax = plt.subplots()
    ot.PlotDiagram(solver.Diagram(), ax, distribution=dist1)

    plt.show()


def PartiallyLeftDomain():

    # Set the initial circle position
    loc0 = np.array([0.6,0.3])
    dens_vals0 = CreateCircleDensity(loc0, circle_radius)
    dist0 = ot.DiscretizedDistribution(grid, dens_vals0)
    nnz0 = np.sum(dens_vals0>0)

    # Set the second position
    loc1 = np.array([0.75*circle_radius, 0.7])
    dens_vals1 = CreateCircleDensity(loc1, circle_radius)
    nnz1 = np.sum(dens_vals1>0)

    area_frac = nnz1/nnz0
    dens_vals1 *= area_frac

    print(np.max(dens_vals1),np.max(dens_vals0))

    dist1 = ot.DiscretizedDistribution(grid, dens_vals1)

    # Generate some random initial points inside the initial circle
    numPts = 4
    probs = (1.0/numPts)*np.ones(numPts)

    # Create an (balanced) optimal quantization of the initial circle
    opts = {'Lloyd Steps':200, 'Lloyd Tol':1e-5, 'GTol Abs':1e-8, 'Max Steps': 300}
    seed_pts = GenerateSeedPts(numPts, dens_vals0)
    diag = ot.SemidiscreteOT.BuildCentroidal(dist0, seed_pts, probs, opts)
    pts = diag.Centroids(dist0)

    # Compute the transport between the first points and the
    solver = ot.QuadraticRegularizedSDOT(dist1, pts, probs, 1.0)

    opts = {'Max Steps':500, 'GTol Abs':3e-2, 'FTol Abs':0.0}
    h = loc0-loc1
    newPrices = np.ones(numPts)#diag.Prices() -  2*np.sum(np.tile(h.reshape((2,1)),(1,numPts))*pts,axis=0)# Theorem 2.3 of James Ronan's thesis
    optPrices, optObj = solver.Solve(newPrices, opts)

    # Estimate the velocity
    diag2 = solver.Diagram()
    true_velocity = loc1-loc0
    ot_velocity = diag2.Centroids(dist1)-pts
    print('From:')
    print(pts)
    print('To:')
    print(diag2.Centroids(dist1))
    print('True Velocity:\n',true_velocity)
    print('OT Velocity:\n',ot_velocity)


    # Use the results to
    fig, ax = plt.subplots()
    ot.PlotDiagram(solver.Diagram(), ax, distribution=dist1)

    plt.show()

PartiallyLeftDomain()
