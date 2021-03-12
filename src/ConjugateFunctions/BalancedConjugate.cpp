#include "SDOT/ConjugateFunctions/BalancedConjugate.h"

using namespace sdot;



double BalancedConjugate::Evaluate(double z){
  return z;
};

double BalancedConjugate::Derivative(double z){
  return 1;
};

double BalancedConjugate::TriangularIntegral(double wi,
                                             Eigen::Vector2d const& xi,
                                             Eigen::Vector2d const& pt1,
                                             Eigen::Vector2d const& pt2,
                                             Eigen::Vector2d const& pt3)
{
  double triInt = (1.0/2.0)*std::pow(xi[0], 2) - 1.0/3.0*xi[0]*pt1[0] - 1.0/3.0*xi[0]*pt2[0] - 1.0/3.0*xi[0]*pt3[0] + (1.0/2.0)*std::pow(xi[1], 2) - 1.0/3.0*xi[1]*pt1[1] - 1.0/3.0*xi[1]*pt2[1] - 1.0/3.0*xi[1]*pt3[1] + (1.0/12.0)*std::pow(pt1[0], 2) + (1.0/12.0)*pt1[0]*pt2[0] + (1.0/12.0)*pt1[0]*pt3[0] + (1.0/12.0)*std::pow(pt2[0], 2) + (1.0/12.0)*pt2[0]*pt3[0] + (1.0/12.0)*std::pow(pt3[0], 2) + (1.0/12.0)*std::pow(pt1[1], 2) + (1.0/12.0)*pt1[1]*pt2[1] + (1.0/12.0)*pt1[1]*pt3[1] + (1.0/12.0)*std::pow(pt2[1], 2) + (1.0/12.0)*pt2[1]*pt3[1] + (1.0/12.0)*std::pow(pt3[1], 2);
  triInt *= 0.5*((pt2[0]-pt1[0])*(pt3[1]-pt1[1]) - (pt3[0]-pt1[0])*(pt2[1]-pt1[1]));

  return -1.0*std::abs(triInt) + wi*TriangularIntegralDeriv(wi,xi,pt1,pt2,pt3);
}

double BalancedConjugate::RectangularIntegral(double wi,
                                              Eigen::Vector2d const& xi,
                                              Eigen::Vector2d const& bottomLeft,
                                              Eigen::Vector2d const& topRight)
{
  double rectInt = (0.5/3.0)*(topRight[1]-bottomLeft[1])*(std::pow(topRight[0]-xi[0],3.0)-std::pow(bottomLeft[0]-xi[0],3.0))
                 + (0.5/3.0)*(topRight[0]-bottomLeft[0])*(std::pow(topRight[1]-xi[1],3.0)-std::pow(bottomLeft[1]-xi[1],3.0));

  return -1.0*std::abs(rectInt) + wi*RectangularIntegralDeriv(wi,xi,bottomLeft, topRight);
}

double BalancedConjugate::TriangularIntegralDeriv(double wi,
                                                  Eigen::Vector2d const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  Eigen::Vector2d const& pt3)
{
  return 0.5*std::abs((pt2[0]*pt1[1]-pt1[0]*pt2[1])+(pt3[0]*pt2[1]-pt2[0]*pt3[1])+(pt1[0]*pt3[1]-pt3[0]*pt1[1]));
}

double BalancedConjugate::RectangularIntegralDeriv(double wi,
                                                   Eigen::Vector2d const& xi,
                                                   Eigen::Vector2d const& bottomLeft,
                                                   Eigen::Vector2d const& topRight)
{
  return std::abs( (topRight[0]-bottomLeft[0])*(topRight[1]-bottomLeft[1]) );
}

double BalancedConjugate::TriangularIntegralDeriv2(double wi,
                                                   Eigen::Vector2d const& xi,
                                                   Eigen::Vector2d const& pt1,
                                                   Eigen::Vector2d const& pt2,
                                                   Eigen::Vector2d const& pt3)
{
  return 0;
}

double BalancedConjugate::RectangularIntegralDeriv2(double wi,
                                                    Eigen::Vector2d const& xi,
                                                    Eigen::Vector2d const& lowerLeft,
                                                    Eigen::Vector2d const& upperRight)
{
  return 0;
}
