#include "SDOT/Integrands/TransportIntegrand.h"

using namespace sdot;

TransportIntegrand::TransportIntegrand(Eigen::Vector2d const& pIn) : p(pIn){}

double TransportIntegrand::TriangularIntegral(Eigen::Vector2d const& pt1,
                                              Eigen::Vector2d const& pt2,
                                              Eigen::Vector2d const& pt3)
{
  double triInt = (1.0/2.0)*std::pow(p[0], 2) - 1.0/3.0*p[0]*pt1[0] - 1.0/3.0*p[0]*pt2[0] - 1.0/3.0*p[0]*pt3[0] + (1.0/2.0)*std::pow(p[1], 2) - 1.0/3.0*p[1]*pt1[1] - 1.0/3.0*p[1]*pt2[1] - 1.0/3.0*p[1]*pt3[1] + (1.0/12.0)*std::pow(pt1[0], 2) + (1.0/12.0)*pt1[0]*pt2[0] + (1.0/12.0)*pt1[0]*pt3[0] + (1.0/12.0)*std::pow(pt2[0], 2) + (1.0/12.0)*pt2[0]*pt3[0] + (1.0/12.0)*std::pow(pt3[0], 2) + (1.0/12.0)*std::pow(pt1[1], 2) + (1.0/12.0)*pt1[1]*pt2[1] + (1.0/12.0)*pt1[1]*pt3[1] + (1.0/12.0)*std::pow(pt2[1], 2) + (1.0/12.0)*pt2[1]*pt3[1] + (1.0/12.0)*std::pow(pt3[1], 2);
  triInt *= 0.5*((pt2[0]-pt1[0])*(pt3[1]-pt1[1]) - (pt3[0]-pt1[0])*(pt2[1]-pt1[1]));

  return std::abs(triInt);
}

double TransportIntegrand::RectangularIntegral(Eigen::Vector2d const& bottomLeft,
                                               Eigen::Vector2d const& topRight)
{
  double rectInt = (0.5/3.0)*(topRight[1]-bottomLeft[1])*(std::pow(topRight[0]-p[0],3.0)-std::pow(bottomLeft[0]-p[0],3.0))
                 + (0.5/3.0)*(topRight[0]-bottomLeft[0])*(std::pow(topRight[1]-p[1],3.0)-std::pow(bottomLeft[1]-p[1],3.0));

  return std::abs(rectInt);
}
