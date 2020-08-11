#include "SDOT/Integrands/ConstantIntegrand.h"

using namespace sdot;

ConstantIntegrand::ConstantIntegrand(double constVal) : val(constVal){};


double ConstantIntegrand::TriangularIntegral(Eigen::Vector2d const& pt1,
                                             Eigen::Vector2d const& pt2,
                                             Eigen::Vector2d const& pt3)
{
  return 0.5*val*std::abs((pt2[0]*pt1[1]-pt1[0]*pt2[1])+(pt3[0]*pt2[1]-pt2[0]*pt3[1])+(pt1[0]*pt3[1]-pt3[0]*pt1[1]));
}


double ConstantIntegrand::RectangularIntegral(Eigen::Vector2d const& bottomLeft,
                                              Eigen::Vector2d const& topRight)
{
  return val * (topRight[0]-bottomLeft[0])*(topRight[1]-bottomLeft[1]);
}
