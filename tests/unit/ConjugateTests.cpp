#include "SDOT/ConjugateFunctions/BalancedConjugate.h"

#include <gtest/gtest.h>

using namespace sdot;

TEST(ConjugateFunctions, BalancedConjugate)
{
  Eigen::Vector2d bottomLeft{0.0,0.0};
  Eigen::Vector2d bottomRight{1.0,0.0};
  Eigen::Vector2d topRight{1.0,1.0};
  Eigen::Vector2d topLeft{0.0,1.0};

  Eigen::Vector2d xi{0.25,0.5};
  double wi = 0.5;

  double val = BalancedConjugate::RectangularIntegral(wi, xi, bottomLeft, topRight);
  double truth = wi - ( std::pow(1.0-xi[0],3.0) + std::pow(xi[0],3.0) + std::pow(1.0-xi[1],3.0) + std::pow(xi[1],3.0))/6.0;
  EXPECT_DOUBLE_EQ(truth, val);

  val = BalancedConjugate::TriangularIntegral(wi, xi, bottomLeft, bottomRight, topRight);
  truth = 0.5*wi - 0.25*xi[0]*xi[0] + xi[0]/3.0 -0.25*xi[1]*xi[1] +xi[1]/6.0 - 1.0/6.0;
  EXPECT_DOUBLE_EQ(truth, val);

  val = BalancedConjugate::TriangularIntegral(wi, xi, bottomLeft, bottomRight, topLeft);
  truth = 0.5*wi - 0.25*xi[0]*xi[0] + xi[0]/6.0 -0.25*xi[1]*xi[1] +xi[1]/6.0 - 1.0/12.0;
  EXPECT_DOUBLE_EQ(truth, val);

  val = BalancedConjugate::TriangularIntegralDeriv(wi,xi, bottomLeft, bottomRight, topLeft);
  EXPECT_DOUBLE_EQ(0.5, val);
}
