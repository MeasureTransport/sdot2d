
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include <vector>
#include <memory>

#include "SDOT/SemidiscreteOT.h"
#include "SDOT/PolygonRasterize.h"
#include "SDOT/RegularGrid.h"
#include "SDOT/DiscretizedDistribution.h"

using namespace sdot;


void AddCircle(Eigen::MatrixXd &dens, double x, double y, double r, double dx, double dy)
{
  for(int j=0; j<dens.cols(); ++j){
    double yj = double(j)*dy;//  + 0.5*dy;

    for(int i=0; i<dens.rows(); ++i){
      double xi = double(i)*dx;// + 0.5*dx;

      if((xi-x)*(xi-x) + (yj-y)*(yj-y) < r*r)
        dens(i,j) = 1.0;
    }
  }
}


int main(int argc, char* argv[])
{
  int N = 20;

  Eigen::Matrix2Xd domain(2,4);
  domain << 0.0, 2.0, 2.0, 0.0,
            0.0, 0.0, 2.0, 2.0;

  // Construct the continuous distribution
  auto grid = std::make_shared<RegularGrid>(domain(0,0),domain(1,0), domain(0,2), domain(1,2), N, N);

  // Unnormalized density.  Will be normalized in DiscretizedDistribution constructor
  Eigen::MatrixXd density = Eigen::MatrixXd::Zero(grid->NumCells(0), grid->NumCells(1));

  double radius = 0.1001;
  AddCircle(density, 0.4,0.8, radius, grid->dx, grid->dy);
  AddCircle(density, 1.4,0.4, radius, grid->dx, grid->dy);

  std::cout << "Density = \n" << density << std::endl;

  Eigen::MatrixXd pts(2,6);
  pts << 1.26996324, 1.11996324, 1.11996324, 0.83003676, 0.68003676, 0.68003676,
         1.64491023, 1.73151277, 1.55830769, 0.60508977, 0.69169231, 0.51848723;
  unsigned int numPts = pts.cols();

  auto dist = std::make_shared<DiscretizedDistribution>(grid, density);

  // Evalaute the SDOT objective
  Eigen::VectorXd discrProbs = Eigen::VectorXd::Ones(pts.cols());
  discrProbs /= pts.cols();

  auto sdot = std::make_shared<SemidiscreteOT>(dist, pts, discrProbs);

  Eigen::VectorXd optPrices;
  double optVal;
  std::tie(optPrices,optVal) = sdot->Solve(Eigen::VectorXd::Ones(numPts));
  std::cout << "Optimal prices = " << optPrices.transpose() << std::endl;

  std::shared_ptr<LaguerreDiagram> lagDiag = sdot->Diagram();

  std::cout << "Laguerre cells = " << std::endl;
  for(int polyInd=0; polyInd<numPts; ++polyInd){
    std::cout << "Working on " <<  polyInd << std::endl;
    std::shared_ptr<PolygonRasterizeIter::Polygon_2> poly = lagDiag->GetCell(polyInd)->ToCGAL();

    if(poly->size()>0){
      auto vertIt = poly->vertices_begin();
      std::cout << "[[" << vertIt->x() << "," << vertIt->y() << "]";
      vertIt++;
      for(;  vertIt != poly->vertices_end(); ++vertIt){
        std::cout << ", [" << vertIt->x() << "," << vertIt->y() << "]";
      }
      std::cout << "]" << std::endl;
    }
  }

  return 0;
}
