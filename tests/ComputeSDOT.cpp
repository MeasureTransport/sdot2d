
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include <vector>
#include <memory>

#include "SDOT/SemiDiscreteOT.h"
#include "SDOT/PolygonRasterize.h"
#include "SDOT/RegularGrid.h"
#include "SDOT/DiscretizedDistribution.h"

using namespace sdot;


void LaguerreDiagramTest()
{
  int numPts = 2;

  Eigen::VectorXd prices = Eigen::VectorXd::Ones(numPts);
  prices << 1, 0.873;//, 0.978674;//, 0.851588;//, 0.996466;

  Eigen::Matrix2Xd pts(2,numPts);
  pts << 0.823295, 0.329554,
         0.604897, 0.536459;
  //pts = Eigen::Matrix2Xd::Random(2,numPts).cwiseAbs();

  Eigen::VectorXd discrProbs(numPts);
  discrProbs = (1.0/numPts)*Eigen::VectorXd::Ones(numPts);

  LaguerreDiagram lagDiag(0, 1, 0, 1, pts, prices);

  LaguerreDiagram::Point_2 startPt = std::get<1>(lagDiag.InternalEdges(0).at(0));
  LaguerreDiagram::Point_2 endPt = std::get<2>(lagDiag.InternalEdges(0).at(0));

  std::cout << "Dividing edge = (" << CGAL::to_double(startPt.x()) << "," << CGAL::to_double(startPt.y()) << ") -> (" << CGAL::to_double(endPt.x()) << "," << CGAL::to_double(endPt.y()) << ")" << std::endl;
  LaguerreDiagram::Point_2 midPt = startPt + 0.5*(endPt-startPt);

  double dist1 = std::sqrt( CGAL::to_double(CGAL::squared_distance(startPt,midPt)));
  double dist2 = std::sqrt( CGAL::to_double(CGAL::squared_distance(endPt,midPt)));

  assert( std::abs( (dist1 - prices(0)) - (dist2 - prices(0)))<1e-14 );
}


int main(int argc, char* argv[])
{
  LaguerreDiagramTest();

  int numPts = 100;

  Eigen::VectorXd prices = Eigen::VectorXd::Ones(numPts);

  Eigen::Matrix2Xd pts = Eigen::Matrix2Xd::Random(2,numPts).cwiseAbs();

  std::cout << "Points = \n";
  std::cout << "[[" << pts(0,0) << "," << pts(1,0) << "]";
  for(int i=1; i<numPts; ++i){
    std::cout << ", [" << pts(0,i) << "," << pts(1,i) << "]";
  }
  std::cout << "]" << std::endl;


  Eigen::VectorXd discrProbs(numPts);
  discrProbs << (1.0/numPts)*Eigen::VectorXd::Ones(numPts);

  Eigen::Matrix2Xd domain(2,4);
  domain << 0.0, 1.0, 1.0, 0.0,
            0.0, 0.0, 1.0, 1.0;

  // Construct the continuous distribution
  auto grid = std::make_shared<RegularGrid>(domain(0,0),domain(1,0), domain(0,2), domain(1,2), 10,10);

  // Unnormalized density.  Will be normalized in DiscretizedDistribution constructor
  Eigen::MatrixXd density = Eigen::MatrixXd::Ones(grid->NumCells(0), grid->NumCells(1));

  auto dist = std::make_shared<DiscretizedDistribution>(grid, density);

  // Evalaute the SDOT objective
  auto sdot = std::make_shared<SemidiscreteOT>(dist, pts, discrProbs);

  Eigen::VectorXd optPrices;
  double optVal;
  std::tie(optPrices,optVal) = sdot->Solve(Eigen::VectorXd::Ones(numPts));
  std::cout << "Optimal prices = " << optPrices.transpose() << std::endl;

  std::shared_ptr<LaguerreDiagram> lagDiag = sdot->Diagram();

  std::cout << "Laguerre cells = " << std::endl;
  for(int polyInd=0; polyInd<numPts; ++polyInd){
    std::shared_ptr<PolygonRasterizeIter::Polygon_2> poly = lagDiag->GetCell(polyInd)->ToCGAL();

    auto vertIt = poly->vertices_begin();
    std::cout << "[[" << vertIt->x() << "," << vertIt->y() << "]";
    vertIt++;
    for(;  vertIt != poly->vertices_end(); ++vertIt){
      std::cout << ", [" << vertIt->x() << "," << vertIt->y() << "]";
    }
    std::cout << "]" << std::endl;
  }

  return 0;
}
