#include "SDOT/Distances/QuadraticRegularization.h"

#include "SDOT/Distances/LineQuadrature.h"
#include "SDOT/Distances/TriangularQuadrature.h"
#include "SDOT/Distances/RectangularQuadrature.h"

using namespace sdot;
using namespace sdot::distances;


double QuadraticRegularization::Evaluate(double z, double penaltyCoeff){
  if(z > -2*penaltyCoeff){
    return (0.25/penaltyCoeff)*z*z+ z;
  }else{
    return -penaltyCoeff;
  }
};

double QuadraticRegularization::Derivative(double z, double penaltyCoeff){
  if(z > -2*penaltyCoeff){
    return (0.5/penaltyCoeff)*z + 1.0;
  }else{
    return 0;
  }
};

double QuadraticRegularization::Derivative2(double z, double penaltyCoeff){
  if(z > -2*penaltyCoeff){
    return 0.5/penaltyCoeff;
  }else{
    return 0;
  }
};

double QuadraticRegularization::TriangularIntegral(double wi,
                                        Eigen::Ref<const Eigen::Vector2d> const& xi,
                                        Eigen::Vector2d const& pt1,
                                        Eigen::Vector2d const& pt2,
                                        Eigen::Vector2d const& pt3,
                                        double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return QuadraticRegularization::Evaluate(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return TriangularQuadrature::Integrate(func,pt1,pt2,pt3);
}

double QuadraticRegularization::RectangularIntegral(double wi,
                                         Eigen::Ref<const Eigen::Vector2d> const& xi,
                                         Eigen::Vector2d const& bottomLeft,
                                         Eigen::Vector2d const& topRight,
                                         double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return QuadraticRegularization::Evaluate(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return RectangularQuadrature::Integrate(func,bottomLeft,topRight);
}

double QuadraticRegularization::TriangularIntegralDeriv(double wi,
                                             Eigen::Ref<const Eigen::Vector2d> const& xi,
                                             Eigen::Vector2d const& pt1,
                                             Eigen::Vector2d const& pt2,
                                             Eigen::Vector2d const& pt3,
                                             double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return QuadraticRegularization::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return TriangularQuadrature::Integrate(func,pt1,pt2,pt3);
}

double QuadraticRegularization::RectangularIntegralDeriv(double wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& bottomLeft,
                                              Eigen::Vector2d const& topRight,
                                              double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return QuadraticRegularization::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return RectangularQuadrature::Integrate(func,bottomLeft,topRight);
}

double QuadraticRegularization::LineIntegralDeriv(double                 wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return QuadraticRegularization::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return LineQuadrature::Integrate(func, pt1,pt2);
}


double QuadraticRegularization::TriangularIntegralDeriv2(double                 wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& pt1,
                                              Eigen::Vector2d const& pt2,
                                              Eigen::Vector2d const& pt3,
                                              double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return QuadraticRegularization::Derivative2(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return TriangularQuadrature::Integrate(func,pt1,pt2,pt3);
}

double QuadraticRegularization::RectangularIntegralDeriv2(double                 wi,
                                               Eigen::Ref<const Eigen::Vector2d> const& xi,
                                               Eigen::Vector2d const& lowerLeft,
                                               Eigen::Vector2d const& upperRight,
                                               double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return QuadraticRegularization::Derivative2(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return RectangularQuadrature::Integrate(func,lowerLeft,upperRight);
}


Eigen::Vector2d QuadraticRegularization::TriangularIntegralPointGrad(double wi,
                                          Eigen::Ref<const Eigen::Vector2d> const& xi,
                                          Eigen::Vector2d const& pt1,
                                          Eigen::Vector2d const& pt2,
                                          Eigen::Vector2d const& pt3,
                                          double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return (2.0*(x-xi).eval()*QuadraticRegularization::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff)).eval();
    };

  return TriangularQuadrature::Integrate(func, pt1,pt2,pt3);
}



Eigen::Vector2d QuadraticRegularization::RectangularIntegralPointGrad(double wi,
                                           Eigen::Ref<const Eigen::Vector2d> const& xi,
                                           Eigen::Vector2d const& pt1,
                                           Eigen::Vector2d const& pt2,
                                           double penaltyCoeff)
{
 auto func = [&](Eigen::Vector2d const& x)
   {
     return (2.0*(x-xi)*QuadraticRegularization::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff)).eval();
   };

 return RectangularQuadrature::Integrate(func, pt1,pt2);
}

Eigen::Matrix2d QuadraticRegularization::TriangularIntegralPointHessDiag(double wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  Eigen::Vector2d const& pt3,
                                                  double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      double wc = wi - (x - xi).squaredNorm();
      return (QuadraticRegularization::Derivative2(wc, penaltyCoeff)*(x-xi)*(x-xi).transpose() - QuadraticRegularization::Derivative(wc, penaltyCoeff)*Eigen::Matrix2d::Identity() ).eval();
    };

  return -2.0*TriangularQuadrature::Integrate(func, pt1,pt2,pt3);
}

Eigen::Matrix2d QuadraticRegularization::RectangularIntegralPointHessDiag(double wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      double wc = wi - (x - xi).squaredNorm();
      return (QuadraticRegularization::Derivative2(wc, penaltyCoeff)*(x-xi)*(x-xi).transpose() - QuadraticRegularization::Derivative(wc, penaltyCoeff)*Eigen::Matrix2d::Identity() ).eval();
    };

  return -2.0*RectangularQuadrature::Integrate(func, pt1,pt2);
}



Eigen::Matrix2d QuadraticRegularization::LineIntegralPointHess(double                 wi,
                                Eigen::Ref<const Eigen::Vector2d> const& xi,
                                Eigen::Ref<const Eigen::Vector2d> const& xj,
                                Eigen::Vector2d const& pt1,
                                Eigen::Vector2d const& pt2,
                                double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return Eigen::Matrix2d(QuadraticRegularization::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff)*(x-xj)*(x-xi).transpose());
    };

  return -2.0*LineQuadrature::Integrate(func, pt1, pt2);
}
