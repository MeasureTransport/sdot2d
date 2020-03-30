
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

int main(int argc, char* argv[])
{
  int numPts = 3;
  Eigen::VectorXd prices(numPts);
  prices << 1.0, 1.0, 1.0;

  Eigen::Matrix2Xd pts(2,numPts);
  pts << 0.1, 0.2, 0.1,
         0.1, 0.1, 0.2;

  Eigen::VectorXd discrProbs(numPts);
  discrProbs = (1.0/numPts)*Eigen::VectorXd::Ones(numPts);

  Eigen::Matrix2Xd domain(2,4);
  domain << 0.0, 0.2, 0.2, 0.0,
            0.0, 0.0, 0.2, 0.2;

  // Construct the continuous distribution
  auto grid = std::make_shared<RegularGrid>(domain(0,0),domain(1,0), domain(0,2), domain(1,2), 100, 100);

  // Unnormalized density.  Will be normalized in DiscretizedDistribution constructor
  Eigen::MatrixXd density = Eigen::MatrixXd::Ones(grid->NumCells(0), grid->NumCells(1));

  auto dist = std::make_shared<DiscretizedDistribution>(grid, density);

  // Evalaute the SDOT objective
  auto sdot = std::make_shared<SemidiscreteOT>(dist, pts, discrProbs);

  double obj, fdObj;
  Eigen::VectorXd grad, fdGrad, tempGrad;
  Eigen::SparseMatrix<double> hess, tempHess;
  Eigen::MatrixXd fdHess;

  const double stepSize = 0.2;
  for(int i=0; i<100; ++i){
    std::cout << "Iteration " << i << std::endl;
    std::cout << "  Prices = " << prices.transpose() << std::endl;
    std::tie(obj, grad, hess) = sdot->Objective(prices);
    std::cout << "  Objective = " << obj << std::endl;
    std::cout << "  Gradient = " << grad.transpose() << std::endl;

    // // Finite difference hessian
    // double fdstep = 0.0;
    // fdGrad = Eigen::VectorXd::Zero(numPts);
    // fdHess = Eigen::MatrixXd::Zero(numPts,numPts);
    //
    // for(unsigned int d=0; d<numPts; ++d){
    //   Eigen::VectorXd newPrices = prices;
    //   newPrices(d) += fdstep;
    //   std::tie(fdObj, tempGrad, tempHess) = sdot->Objective(newPrices);
    //   std::cout << "fdObj(" << d << ")=" << fdObj << std::endl;
    //   fdGrad(d) = (fdObj-obj)/fdstep;
    //   fdHess.col(d) = (tempGrad-grad)/fdstep;
    // }
    // //fdHess = 0.5*(fdHess  + fdHess.transpose()).eval();
    // std::cout << "Analytic gradient = " << grad.transpose() << std::endl;
    // std::cout << "FD gradient = " << fdGrad.transpose() << std::endl;
    //
    // std::cout << "Analytic Hessian = \n" << hess*Eigen::MatrixXd::Identity(numPts,numPts) <<std::endl;
    // std::cout << "FD Hessian = \n" << fdHess << std::endl;
    // exit(0);
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver(-1.0*hess);

    //std::cout << "  Inverse Hessian = \n" << -1.0*solver.solve(Eigen::MatrixXd::Identity(numPts,numPts)) << std::endl;

    Eigen::VectorXd step = stepSize*solver.solve(grad);
    std::cout << "  Step = " << step.transpose() << std::endl;

    prices.tail(numPts-1) += stepSize*step.tail(numPts-1);
  }

  std::cout << "Objective = " << obj << std::endl;

  return 0;
}
