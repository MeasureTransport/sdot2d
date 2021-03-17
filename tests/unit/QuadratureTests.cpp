#include "SDOT/Distances/LineQuadrature.h"
#include "SDOT/Distances/RectangularQuadrature.h"
#include "SDOT/Distances/TriangularQuadrature.h"
#include "SDOT/Distances/QuadraticRegularization.h"

#include <gtest/gtest.h>

using namespace sdot::distances;


TEST(Quadrature, Line){

  Eigen::VectorXd pts, wts;
  std::tie(pts,wts) = LineQuadrature::Get(5);

  EXPECT_NEAR(1.0, wts.sum(),  1e-8);

  double scale = 2.0;
  Eigen::VectorXd vals = (scale*pts-0.25*Eigen::VectorXd::Ones(wts.size())).array().pow(5);
  double integral = scale*vals.dot(wts);

  double trueIntegral = (std::pow(1.75,6.0) - std::pow(-0.25,6.0))/6.0;
  EXPECT_NEAR(trueIntegral,integral,1e-8);
}

TEST(Quadrature, Square){

  Eigen::Matrix2Xd pts;
  Eigen::VectorXd wts;
  std::tie(pts,wts) = RectangularQuadrature::Get(4);

  EXPECT_NEAR(1.0, wts.sum(),  1e-8);

  Eigen::VectorXd vals(wts.size());
  for(int  i=0; i<wts.size();++i){
    vals(i) = pts.col(i).squaredNorm();
  }

  double integral = vals.dot(wts);
  EXPECT_NEAR(2.0/3.0,integral,1e-8);
}


TEST(Quadrature, Triangle){

  Eigen::Matrix2Xd pts;
  Eigen::VectorXd wts;
  std::tie(pts,wts) = TriangularQuadrature::Get(7);

  EXPECT_NEAR(0.5, wts.sum(),  1e-8);

  Eigen::VectorXd vals(wts.size());
  for(int  i=0; i<wts.size();++i){
    vals(i) = pts.col(i).squaredNorm();
  }

  double integral = vals.dot(wts);
  EXPECT_NEAR(1.0/6.0, integral, 1e-8);
}

TEST(Quadrature, GenericTriangle){

  auto f = [](double z) { return 1.0; };

  Eigen::Vector2d xi{0.0,0.0};
  Eigen::Vector2d pt1{0.0,0.0}, pt2{1.0,0.0}, pt3{0.0,1.0};

  double res = QuadraticRegularization::GenericTriangularIntegral(f, 0.0, xi, pt1, pt2, pt3);
  EXPECT_NEAR(0.5, res, 1e-8);

  auto f2 = [](double z) { return -z; };
  res = QuadraticRegularization::GenericTriangularIntegral(f2, 0.0, xi, pt1, pt2, pt3);
  EXPECT_NEAR(1.0/6.0, res, 1e-8);
}

TEST(Quadrature, GenericRectangularIntegral){

  auto f = [](double z) { return 1.0; };

  Eigen::Vector2d xi{0.0,0.0};
  Eigen::Vector2d pt1{0.0,0.0}, pt2{2.0,2.0};

  double res = QuadraticRegularization::GenericRectangularIntegral(f, 0.0, xi, pt1, pt2);
  EXPECT_NEAR(4.0, res, 1e-8);

  auto f2 = [](double z) { return -z; };
  res = QuadraticRegularization::GenericRectangularIntegral(f2, 0.0, xi, pt1, pt2);
  EXPECT_NEAR(32.0/3.0, res, 1e-8);
}
