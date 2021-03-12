#include "SDOT/Integrands/ConstantIntegrand.h"
#include "SDOT/Integrands/TransportIntegrand.h"

#include <gtest/gtest.h>

using namespace sdot;

TEST(Integrands, Constant) {

    // Set the corners of a unit square
    Eigen::Vector2d bottomLeft{0.0,0.0};
    Eigen::Vector2d bottomRight{1.0,0.0};
    Eigen::Vector2d topRight{1.0,1.0};
    Eigen::Vector2d topLeft{0.0,1.0};
    Eigen::Vector2d center{0.5,0.5};

    double constVal = 2.0;
    ConstantIntegrand test(constVal);


    // The integral of the bottom left triangle
    double val = test.TriangularIntegral(bottomLeft, bottomRight, topLeft);
    EXPECT_DOUBLE_EQ(0.5*constVal, val);

    val = test.TriangularIntegral(bottomLeft, bottomRight, topRight);
    EXPECT_DOUBLE_EQ(0.5*constVal, val);

    val = test.TriangularIntegral(bottomLeft, bottomRight, center);
    EXPECT_DOUBLE_EQ(0.25*constVal, val);

    val = test.TriangularIntegral(center, topLeft, topRight);
    EXPECT_DOUBLE_EQ(0.25*constVal, val);

    val = test.RectangularIntegral(center, topRight);
    EXPECT_DOUBLE_EQ(0.25*constVal, val);

    val = test.RectangularIntegral(bottomLeft, center);
    EXPECT_DOUBLE_EQ(0.25*constVal, val);

    val = test.RectangularIntegral(bottomLeft, topRight);
    EXPECT_DOUBLE_EQ(constVal, val);
}


TEST(Integrands, Transport) {

  // Set the corners of a unit square
  Eigen::Vector2d bottomLeft{0.0,0.0};
  Eigen::Vector2d bottomRight{1.0,0.0};
  Eigen::Vector2d topRight{1.0,1.0};
  Eigen::Vector2d topLeft{0.0,1.0};
  Eigen::Vector2d center{0.5,0.5};


  TransportIntegrand test(center);

  double val = test.RectangularIntegral(bottomLeft, topRight);
  EXPECT_DOUBLE_EQ(1.0/12.0, val);

  val = test.TriangularIntegral(bottomLeft, bottomRight, topRight);
  EXPECT_DOUBLE_EQ(1.0/24.0, val);

  val = test.TriangularIntegral(bottomLeft, center, bottomRight);
  EXPECT_DOUBLE_EQ(1.0/48.0, val);
}
