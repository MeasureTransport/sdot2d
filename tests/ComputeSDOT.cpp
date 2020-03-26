
#include <Eigen/Core>

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
  auto grid = std::make_shared<RegularGrid>(domain(0,0),domain(1,0), domain(0,2), domain(1,2), 200, 200);
  Eigen::MatrixXd probs = (1.0/grid->NumCells()) * Eigen::MatrixXd::Identity(grid->NumCells(0), grid->NumCells(1));

  auto dist = std::make_shared<DiscretizedDistribution>(grid, probs);

  // Evalaute the SDOT objective
  auto sdot = std::make_shared<SemidiscreteOT>(dist, pts, discrProbs);

  double obj = sdot->Objective(prices);
  std::cout << "Objective = " << obj << std::endl;

  return 0;
}
