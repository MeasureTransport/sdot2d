#include "SDOT/Distances/Wasserstein2.h"

#include <gtest/gtest.h>

using namespace sdot;
using namespace sdot::distances;

TEST(Distances, Wasserstein2)
{
  Eigen::Vector2d bottomLeft{0.0,0.0};
  Eigen::Vector2d bottomRight{1.0,0.0};
  Eigen::Vector2d topRight{1.0,1.0};
  Eigen::Vector2d topLeft{0.0,1.0};

  Eigen::Vector2d xi{0.25,0.5};
  double wi = 0.5;

  double val = Wasserstein2::RectangularIntegral(wi, xi, bottomLeft, topRight);
  double truth = wi - ( std::pow(1.0-xi[0],3.0) + std::pow(xi[0],3.0) + std::pow(1.0-xi[1],3.0) + std::pow(xi[1],3.0))/6.0;
  EXPECT_DOUBLE_EQ(truth, val);

  val = Wasserstein2::TriangularIntegral(wi, xi, bottomLeft, bottomRight, topRight);
  truth = 0.5*wi - 0.25*xi[0]*xi[0] + xi[0]/3.0 -0.25*xi[1]*xi[1] +xi[1]/6.0 - 1.0/6.0;
  EXPECT_DOUBLE_EQ(truth, val);

  val = Wasserstein2::TriangularIntegral(wi, xi, bottomLeft, bottomRight, topLeft);
  truth = 0.5*wi - 0.25*xi[0]*xi[0] + xi[0]/6.0 -0.25*xi[1]*xi[1] +xi[1]/6.0 - 1.0/12.0;
  EXPECT_DOUBLE_EQ(truth, val);

  val = Wasserstein2::TriangularIntegralDeriv(wi,xi, bottomLeft, bottomRight, topLeft);
  EXPECT_DOUBLE_EQ(0.5, val);

  val = Wasserstein2::LineIntegralDeriv(wi, xi, bottomLeft, topRight);
  EXPECT_DOUBLE_EQ(std::sqrt(2), val);
}
