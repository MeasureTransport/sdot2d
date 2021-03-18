#include "SDOT/Distances/Wasserstein2.h"

#include "SDOT/Distances/LineQuadrature.h"
#include "SDOT/Distances/TriangularQuadrature.h"
#include "SDOT/Distances/RectangularQuadrature.h"

using namespace sdot;
using namespace sdot::distances;


double Wasserstein2::Evaluate(double z, double penaltyCoeff){
  return z;
};

double Wasserstein2::Derivative(double z, double penaltyCoeff){
  return 1;
};

double Wasserstein2::Derivative2(double z, double penaltyCoeff){
  return 0.0;
};

double Wasserstein2::TriangularIntegral(double wi,
                                        Eigen::Ref<const Eigen::Vector2d> const& xi,
                                        Eigen::Vector2d const& pt1,
                                        Eigen::Vector2d const& pt2,
                                        Eigen::Vector2d const& pt3,
                                        double penaltyCoeff)
{

  double triInt = (1.0/2.0)*std::pow(xi[0], 2) - 1.0/3.0*xi[0]*pt1[0] - 1.0/3.0*xi[0]*pt2[0] - 1.0/3.0*xi[0]*pt3[0] + (1.0/2.0)*std::pow(xi[1], 2) - 1.0/3.0*xi[1]*pt1[1] - 1.0/3.0*xi[1]*pt2[1] - 1.0/3.0*xi[1]*pt3[1] + (1.0/12.0)*std::pow(pt1[0], 2) + (1.0/12.0)*pt1[0]*pt2[0] + (1.0/12.0)*pt1[0]*pt3[0] + (1.0/12.0)*std::pow(pt2[0], 2) + (1.0/12.0)*pt2[0]*pt3[0] + (1.0/12.0)*std::pow(pt3[0], 2) + (1.0/12.0)*std::pow(pt1[1], 2) + (1.0/12.0)*pt1[1]*pt2[1] + (1.0/12.0)*pt1[1]*pt3[1] + (1.0/12.0)*std::pow(pt2[1], 2) + (1.0/12.0)*pt2[1]*pt3[1] + (1.0/12.0)*std::pow(pt3[1], 2);
  triInt *= 0.5*((pt2[0]-pt1[0])*(pt3[1]-pt1[1]) - (pt3[0]-pt1[0])*(pt2[1]-pt1[1]));

  return -2.0*std::abs(triInt) + wi*TriangularIntegralDeriv(wi,xi,pt1,pt2,pt3);
}

double Wasserstein2::RectangularIntegral(double wi,
                                         Eigen::Ref<const Eigen::Vector2d> const& xi,
                                         Eigen::Vector2d const& bottomLeft,
                                         Eigen::Vector2d const& topRight,
                                         double penaltyCoeff)
{
  double rectInt = (1.0/3.0)*(topRight[1]-bottomLeft[1])*(std::pow(topRight[0]-xi[0],3.0)-std::pow(bottomLeft[0]-xi[0],3.0))
                 + (1.0/3.0)*(topRight[0]-bottomLeft[0])*(std::pow(topRight[1]-xi[1],3.0)-std::pow(bottomLeft[1]-xi[1],3.0));

  return -1.0*std::abs(rectInt) + wi*RectangularIntegralDeriv(wi,xi,bottomLeft, topRight);
}

double Wasserstein2::TriangularIntegralDeriv(double wi,
                                             Eigen::Ref<const Eigen::Vector2d> const& xi,
                                             Eigen::Vector2d const& pt1,
                                             Eigen::Vector2d const& pt2,
                                             Eigen::Vector2d const& pt3,
                                             double penaltyCoeff)
{
  return 0.5*std::abs((pt2[0]*pt1[1]-pt1[0]*pt2[1])+(pt3[0]*pt2[1]-pt2[0]*pt3[1])+(pt1[0]*pt3[1]-pt3[0]*pt1[1]));
}

double Wasserstein2::RectangularIntegralDeriv(double wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& bottomLeft,
                                              Eigen::Vector2d const& topRight,
                                              double penaltyCoeff)
{
  return std::abs( (topRight[0]-bottomLeft[0])*(topRight[1]-bottomLeft[1]) );
}

double Wasserstein2::LineIntegralDeriv(double                 wi,
                                       Eigen::Ref<const Eigen::Vector2d> const& xi,
                                       Eigen::Vector2d const& pt1,
                                       Eigen::Vector2d const& pt2,
                                       double penaltyCoeff)
{
    return (pt1-pt2).norm();
}


double Wasserstein2::TriangularIntegralDeriv2(double                 wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& pt1,
                                              Eigen::Vector2d const& pt2,
                                              Eigen::Vector2d const& pt3,
                                              double penaltyCoeff)
{
  return 0;
}

double Wasserstein2::RectangularIntegralDeriv2(double                 wi,
                                               Eigen::Ref<const Eigen::Vector2d> const& xi,
                                               Eigen::Vector2d const& lowerLeft,
                                               Eigen::Vector2d const& upperRight,
                                               double penaltyCoeff)
{
  return 0;
}


Eigen::Vector2d Wasserstein2::TriangularIntegralPointGrad(double wi,
                                                 Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                 Eigen::Vector2d const& pt1,
                                                 Eigen::Vector2d const& pt2,
                                                 Eigen::Vector2d const& pt3,
                                                 double penaltyCoeff)
{
  double det   = (pt2(0)-pt1(0))*(pt3(1)-pt1(1)) - (pt3(0)-pt1(0))*(pt2(1)-pt1(1));
  return 2.0*(0.5*(pt2-xi) - (1.0/3.0)*(pt2-pt1)+(1.0/6.0)*(pt3-pt1))*det;
}


Eigen::Vector2d Wasserstein2::RectangularIntegralPointGrad(double wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  double penaltyCoeff)
{
  Eigen::Vector2d out(2);
  out << 0.5*(pt2(1)-pt1(1))*(std::pow(pt2(0)-xi(0),2.0)-std::pow(pt1(0)-xi(0),2.0)),
         0.5*(pt2(0)-pt1(0))*(std::pow(pt2(1)-xi(1),2.0)-std::pow(pt1(1)-xi(1),2.0));
  return 2.0*out;
}


Eigen::Matrix2d Wasserstein2::TriangularIntegralPointHessDiag(double wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  Eigen::Vector2d const& pt3,
                                                  double penaltyCoeff)
{
  return 2.0*Eigen::Matrix2d::Identity()*Wasserstein2::TriangularIntegralDeriv(wi, xi, pt1, pt2, pt3, penaltyCoeff);
}

Eigen::Matrix2d Wasserstein2::RectangularIntegralPointHessDiag(double wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  double penaltyCoeff)
{
  return 2.0*Eigen::Matrix2d::Identity()*Wasserstein2::RectangularIntegralDeriv(wi, xi, pt1, pt2, penaltyCoeff);
}

Eigen::Matrix2d Wasserstein2::LineIntegralPointHess(double                 wi,
                                Eigen::Ref<const Eigen::Vector2d> const& xi,
                                Eigen::Ref<const Eigen::Vector2d> const& xj,
                                Eigen::Vector2d const& pt1,
                                Eigen::Vector2d const& pt2,
                                double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return Eigen::Matrix2d(Wasserstein2::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff)*(x-xj)*(x-xi).transpose());
    };

  return -2.0*LineQuadrature::Integrate(func, pt1, pt2);
}
