#include "SDOT/SemiDiscreteOT.h"

#include "SDOT/RegularGrid.h"
#include "SDOT/DiscretizedDistribution.h"
#include "SDOT/Distances/Wasserstein2.h"
#include "SDOT/Distances/QuadraticRegularization.h"

#include <gtest/gtest.h>

using namespace sdot;
using namespace sdot::distances;

TEST(SemiDiscreteOT, Objective)
{
  Eigen::MatrixXd pts(2,2);
  pts << 0.25, 0.75,
         0.5, 0.5;

  Eigen::VectorXd probs(2);
  probs << 0.75,0.25;

  // Construct the continuous distribution
  BoundingBox bbox(0.0, 1.0, 0.0, 1.0);
  int N = 1;
  auto grid = std::make_shared<RegularGrid>(bbox, N, N);
  auto dist = std::make_shared<DiscretizedDistribution>(grid, Eigen::MatrixXd::Ones(N,N));

  SemidiscreteOT<Wasserstein2> solver(dist, pts, probs);

  // Set the prices and compute the Laguerre diagram
  Eigen::VectorXd prices = 1.25*Eigen::VectorXd::Ones(pts.cols());
  LaguerreDiagram diag(bbox, pts, prices);


  // Compute the dual wass objective
  double obj;
  Eigen::VectorXd grad;
  std::tie(obj,grad) = solver.ComputeGradient(prices, diag);

  Eigen::SparseMatrix<double> hess = solver.ComputeHessian(prices, diag);

  // Truth obtained by analytically evaluating integrals
  double trueObj = 2.0*(1.0 - 0.94791666666666);
  EXPECT_NEAR(trueObj, obj, 1e-8);

  double obj2;
  Eigen::VectorXd grad2;
  double stepSize = 1e-5;
  Eigen::VectorXd stepDir = Eigen::VectorXd::Random(2);
  stepDir /= stepDir.norm();

  Eigen::VectorXd prices2 = prices + stepSize*stepDir;
  LaguerreDiagram diag2(bbox, pts, prices2);

  std::tie(obj2, grad2) = solver.ComputeGradient(prices2, diag2);

  EXPECT_NEAR((obj2-obj)/stepSize, grad.dot(stepDir), 1e-5);

  Eigen::VectorXd fdTruth = (grad2-grad)/stepSize;
  Eigen::VectorXd hessApp = hess * stepDir;
  EXPECT_NEAR(fdTruth(0), hessApp(0), 1e-5);
  EXPECT_NEAR(fdTruth(1), hessApp(1), 1e-5);

}


TEST(SemiDiscreteOT, Objective_QuadraticRegularization)
{
  Eigen::MatrixXd pts(2,2);
  pts << 0.25, 0.75,
         0.5, 0.5;

  Eigen::VectorXd probs(2);
  probs << 0.75,0.25;

  // Construct the continuous distribution
  BoundingBox bbox(0.0, 1.0, 0.0, 1.0);
  int N = 100;
  auto grid = std::make_shared<RegularGrid>(bbox, N, N);
  auto dist = std::make_shared<DiscretizedDistribution>(grid, Eigen::MatrixXd::Ones(N,N));

  SemidiscreteOT<QuadraticRegularization> solver(dist, pts, probs);

  // Set the prices and compute the Laguerre diagram
  Eigen::VectorXd prices = 1.25*Eigen::VectorXd::Ones(pts.cols());
  LaguerreDiagram diag(bbox, pts, prices);


  // Compute the dual wass objective
  double obj;
  Eigen::VectorXd grad;
  std::tie(obj,grad) = solver.ComputeGradient(prices, diag);

  double obj2, objMid;
  Eigen::VectorXd grad2, gradMid;
  double stepSize = 1e-5;
  Eigen::VectorXd stepDir = Eigen::VectorXd::Random(2);
  stepDir /= stepDir.norm();

  Eigen::VectorXd prices2 = prices + stepSize*stepDir;
  LaguerreDiagram diag2(bbox, pts, prices2);

  std::tie(obj2, grad2) = solver.ComputeGradient(prices2, diag2);

  Eigen::VectorXd pricesMid = prices + 0.5*stepSize*stepDir;
  LaguerreDiagram diagMid(bbox, pts, pricesMid);

  std::tie(objMid,gradMid) = solver.ComputeGradient(pricesMid, diagMid);
  Eigen::SparseMatrix<double> hess = solver.ComputeHessian(pricesMid, diagMid);

  EXPECT_NEAR((obj2-obj)/stepSize, gradMid.dot(stepDir), 1e-8);

  Eigen::VectorXd fdTruth = (grad2-grad)/stepSize;
  Eigen::VectorXd hessApp = hess * stepDir;
  EXPECT_NEAR(fdTruth(0), hessApp(0), 1e-3);
  EXPECT_NEAR(fdTruth(1), hessApp(1), 1e-3);
}

TEST(SemiDiscreteOT, Construction)
{
  Eigen::MatrixXd pts(2,2);
  pts << 0.25, 0.75,
         0.5, 0.5;

  Eigen::VectorXd probs(2);
  probs << 0.75,0.25;

  // Construct the continuous distribution
  auto grid = std::make_shared<RegularGrid>(0.0, 0.0, 1.0, 1.0, 1, 1);
  auto dist = std::make_shared<DiscretizedDistribution>(grid, 5.0*Eigen::MatrixXd::Ones(1,1));

  SemidiscreteOT<Wasserstein2> solver(dist, pts, probs);
  solver.Solve(Eigen::VectorXd::Ones(pts.cols()));

  auto diag = solver.Diagram();
  assert(diag != nullptr);

  Eigen::Matrix2Xd verts = diag->GetCellVertices(0);


  EXPECT_NEAR(probs(0), diag->CellArea(0, dist), 1e-3);
  EXPECT_NEAR(probs(1), diag->CellArea(1, dist), 1e-3);
}

TEST(SemiDiscreteOT, Construction_QuadraticRegularization)
{
  Eigen::MatrixXd pts(2,2);
  pts << 0.25, 0.75,
         0.5, 0.5;

  Eigen::VectorXd probs(2);
  probs << 0.75,0.25;

  // Construct the continuous distribution
  auto grid = std::make_shared<RegularGrid>(0.0, 0.0, 1.0, 1.0, 1, 1);
  auto dist = std::make_shared<DiscretizedDistribution>(grid, Eigen::MatrixXd::Ones(1,1));

  SemidiscreteOT<QuadraticRegularization> solver(dist, pts, probs);
  solver.Solve(2.0*Eigen::VectorXd::Ones(pts.cols()));

  auto diag = solver.Diagram();
  assert(diag != nullptr);

}
