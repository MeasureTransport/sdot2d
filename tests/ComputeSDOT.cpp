
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

  int numPts = 1000;

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
    std::shared_ptr<PolygonRasterizeIter::Polygon_2> poly = lagDiag->GetCell(polyInd);

    auto vertIt = poly->vertices_begin();
    std::cout << "[[" << vertIt->x() << "," << vertIt->y() << "]";
    vertIt++;
    for(;  vertIt != poly->vertices_end(); ++vertIt){
      std::cout << ", [" << vertIt->x() << "," << vertIt->y() << "]";
    }
    std::cout << "]" << std::endl;
  }

  exit(0);


  double obj, fdObj;
  Eigen::VectorXd grad, fdGrad, tempGrad;
  Eigen::SparseMatrix<double> hess, tempHess;
  Eigen::MatrixXd fdHess(numPts,numPts);

  const double stepSize = 1.0;
  for(int i=0; i<1; ++i){
    std::cout << "\n\nIteration " << i << std::endl;
    std::cout << "  Base Prices = " << prices.transpose() << std::endl;
    std::tie(obj, grad, hess) = sdot->Objective(prices);
    // std::cout << "  Objective = " << obj << std::endl;
    // std::cout << "  Gradient = " << grad.transpose() << std::endl;
    // // exit(0);

    // double tempObj1, tempObj2;
    // double fdstep = 1e-4;
    // Eigen::VectorXd tempGrad1, tempGrad2, fdGrad(numPts);
    // Eigen::SparseMatrix<double> tempHess;
    // for(unsigned int d=0; d<numPts; ++d){
    //   Eigen::VectorXd newPrices = prices;
    //   newPrices(d) = prices(d) + 0.5*fdstep;
    //   std::cout << "\n  FD Prices = " << newPrices.transpose() << std::endl;
    //   std::tie(tempObj1, tempGrad1, tempHess) = sdot->Objective(newPrices);
    //   std::cout << "\n  FD Prices = " << newPrices.transpose() << std::endl;
    //   newPrices(d) = prices(d) - 0.5*fdstep;
    //   std::tie(tempObj2, tempGrad2, tempHess) = sdot->Objective(newPrices);
    //
    //   fdGrad(d) = (tempObj1 - tempObj2)/fdstep;
    //   fdHess.col(d) = (tempGrad1-tempGrad2)/fdstep;
    // }
    //
    //
    // std::cout << "\nAn Gradient = " << grad.transpose() << std::endl;
    // std::cout << "FD Gradient = " << fdGrad.transpose() << std::endl;
    //
    // std::cout << "An Hessian = \n" << hess*Eigen::MatrixXd::Identity(numPts,numPts) << std::endl;
    // std::cout << "Fd Hessian = \n" << fdHess << std::endl;
    //
    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver(-1.0*hess.block(1,1,numPts-1,numPts-1));
    // //
    // // //std::cout << "  Inverse Hessian = \n" << -1.0*solver.solve(Eigen::MatrixXd::Identity(numPts,numPts)) << std::endl;
    // //
    // Eigen::VectorXd step = -solver.solve(-1.0*grad.tail(numPts-1));
    // std::cout << "  Newton Step = " << step.transpose() << std::endl;
    //
    //
    // prices.tail(numPts-1) += stepSize*step;
  }

  //std::cout << "Objective = " << obj << std::endl;

  return 0;
}
