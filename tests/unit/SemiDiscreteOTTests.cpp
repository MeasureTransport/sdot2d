#include "SDOT/SemiDiscreteOT.h"

#include "SDOT/RegularGrid.h"
#include "SDOT/DiscretizedDistribution.h"
#include "SDOT/Distances/Wasserstein2.h"

#include <gtest/gtest.h>

using namespace sdot;
using namespace sdot::distances;

TEST(SemiDiscreteOT, Construction)
{
  Eigen::MatrixXd pts(2,2);
  pts << 0.25, 0.75,
         0.5, 0.5;

  Eigen::VectorXd probs(2);
  probs << 0.75,0.25;

  // Construct the continuous distribution
  auto grid = std::make_shared<RegularGrid>(0.0, 0.0, 1.0, 1.0, 1, 1);
  auto dist = std::make_shared<DiscretizedDistribution>(grid, Eigen::MatrixXd::Ones(1,1));

  SemidiscreteOT<Wasserstein2> solver(dist, pts, probs);
  solver.Solve(Eigen::VectorXd::Ones(pts.cols()));

  auto diag = solver.Diagram();
  assert(diag != nullptr);

  Eigen::Matrix2Xd verts = diag->GetCellVertices(0);


  EXPECT_NEAR(probs(0), diag->CellArea(0, dist), 1e-3);
  EXPECT_NEAR(probs(1), diag->CellArea(1, dist), 1e-3);
}
