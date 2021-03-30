#include "SDOT/Distances/QuadratureDistance.h"

#include "SDOT/Distances/QuadraticRegularizationFunctions.h"
#include "SDOT/Distances/GHKFunctions.h"

#include "SDOT/Distances/LineQuadrature.h"
#include "SDOT/Distances/TriangularQuadrature.h"
#include "SDOT/Distances/RectangularQuadrature.h"

using namespace sdot;
using namespace sdot::distances;


template<typename ConjugateFunctionType>
double QuadratureDistance<ConjugateFunctionType>::TriangularIntegral(double wi,
                                        Eigen::Ref<const Eigen::Vector2d> const& xi,
                                        Eigen::Vector2d const& pt1,
                                        Eigen::Vector2d const& pt2,
                                        Eigen::Vector2d const& pt3,
                                        double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return ConjugateFunctionType::Evaluate(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return TriangularQuadrature::Integrate(func,pt1,pt2,pt3);
}

template<typename ConjugateFunctionType>
double QuadratureDistance<ConjugateFunctionType>::RectangularIntegral(double wi,
                                         Eigen::Ref<const Eigen::Vector2d> const& xi,
                                         Eigen::Vector2d const& bottomLeft,
                                         Eigen::Vector2d const& topRight,
                                         double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return ConjugateFunctionType::Evaluate(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return RectangularQuadrature::Integrate(func,bottomLeft,topRight);
}

template<typename ConjugateFunctionType>
double QuadratureDistance<ConjugateFunctionType>::TriangularIntegralDeriv(double wi,
                                             Eigen::Ref<const Eigen::Vector2d> const& xi,
                                             Eigen::Vector2d const& pt1,
                                             Eigen::Vector2d const& pt2,
                                             Eigen::Vector2d const& pt3,
                                             double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return TriangularQuadrature::Integrate(func,pt1,pt2,pt3);
}

template<typename ConjugateFunctionType>
double QuadratureDistance<ConjugateFunctionType>::RectangularIntegralDeriv(double wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& bottomLeft,
                                              Eigen::Vector2d const& topRight,
                                              double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return RectangularQuadrature::Integrate(func,bottomLeft,topRight);
}

template<typename ConjugateFunctionType>
double QuadratureDistance<ConjugateFunctionType>::LineIntegralDeriv(double                 wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return LineQuadrature::Integrate(func, pt1,pt2);
}

template<typename ConjugateFunctionType>
double QuadratureDistance<ConjugateFunctionType>::TriangularIntegralDeriv2(double                 wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& pt1,
                                              Eigen::Vector2d const& pt2,
                                              Eigen::Vector2d const& pt3,
                                              double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return ConjugateFunctionType::Derivative2(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return TriangularQuadrature::Integrate(func,pt1,pt2,pt3);
}

template<typename ConjugateFunctionType>
double QuadratureDistance<ConjugateFunctionType>::RectangularIntegralDeriv2(double                 wi,
                                               Eigen::Ref<const Eigen::Vector2d> const& xi,
                                               Eigen::Vector2d const& lowerLeft,
                                               Eigen::Vector2d const& upperRight,
                                               double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return ConjugateFunctionType::Derivative2(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return RectangularQuadrature::Integrate(func,lowerLeft,upperRight);
}

template<typename ConjugateFunctionType>
Eigen::Vector2d QuadratureDistance<ConjugateFunctionType>::TriangularIntegralPointGrad(double wi,
                                          Eigen::Ref<const Eigen::Vector2d> const& xi,
                                          Eigen::Vector2d const& pt1,
                                          Eigen::Vector2d const& pt2,
                                          Eigen::Vector2d const& pt3,
                                          double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return (2.0*(x-xi).eval()*ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff)).eval();
    };

  return TriangularQuadrature::Integrate(func, pt1,pt2,pt3);
}


template<typename ConjugateFunctionType>
Eigen::Vector2d QuadratureDistance<ConjugateFunctionType>::RectangularIntegralPointGrad(double wi,
                                           Eigen::Ref<const Eigen::Vector2d> const& xi,
                                           Eigen::Vector2d const& pt1,
                                           Eigen::Vector2d const& pt2,
                                           double penaltyCoeff)
{
 auto func = [&](Eigen::Vector2d const& x)
   {
     return (2.0*(x-xi)*ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff)).eval();
   };

 return RectangularQuadrature::Integrate(func, pt1,pt2);
}


template<typename ConjugateFunctionType>
double QuadratureDistance<ConjugateFunctionType>::TriangularIntegralMarginalMass(double wi,
                                          Eigen::Ref<const Eigen::Vector2d> const& xi,
                                          Eigen::Vector2d const& pt1,
                                          Eigen::Vector2d const& pt2,
                                          Eigen::Vector2d const& pt3,
                                          double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

  return TriangularQuadrature::Integrate(func, pt1,pt2,pt3);
}


template<typename ConjugateFunctionType>
double QuadratureDistance<ConjugateFunctionType>::RectangularIntegralMarginalMass(double wi,
                                           Eigen::Ref<const Eigen::Vector2d> const& xi,
                                           Eigen::Vector2d const& pt1,
                                           Eigen::Vector2d const& pt2,
                                           double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff);
    };

   return RectangularQuadrature::Integrate(func, pt1,pt2);
}

template<typename ConjugateFunctionType>
Eigen::Vector2d QuadratureDistance<ConjugateFunctionType>::TriangularIntegralMarginalCentroid(double wi,
                                          Eigen::Ref<const Eigen::Vector2d> const& xi,
                                          Eigen::Vector2d const& pt1,
                                          Eigen::Vector2d const& pt2,
                                          Eigen::Vector2d const& pt3,
                                          double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return (x*ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff)).eval();
    };

  return TriangularQuadrature::Integrate(func, pt1,pt2,pt3);
}


template<typename ConjugateFunctionType>
Eigen::Vector2d QuadratureDistance<ConjugateFunctionType>::RectangularIntegralMarginalCentroid(double wi,
                                           Eigen::Ref<const Eigen::Vector2d> const& xi,
                                           Eigen::Vector2d const& pt1,
                                           Eigen::Vector2d const& pt2,
                                           double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return (x*ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff)).eval();
    };

   return RectangularQuadrature::Integrate(func, pt1,pt2);
}

template<typename ConjugateFunctionType>
Eigen::Matrix2d QuadratureDistance<ConjugateFunctionType>::TriangularIntegralPointHessDiag(double wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  Eigen::Vector2d const& pt3,
                                                  double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      double wc = wi - (x - xi).squaredNorm();
      return (ConjugateFunctionType::Derivative2(wc, penaltyCoeff)*(x-xi)*(x-xi).transpose() - ConjugateFunctionType::Derivative(wc, penaltyCoeff)*Eigen::Matrix2d::Identity() ).eval();
    };

  return -2.0*TriangularQuadrature::Integrate(func, pt1,pt2,pt3);
}

template<typename ConjugateFunctionType>
Eigen::Matrix2d QuadratureDistance<ConjugateFunctionType>::RectangularIntegralPointHessDiag(double wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      double wc = wi - (x - xi).squaredNorm();
      return (ConjugateFunctionType::Derivative2(wc, penaltyCoeff)*(x-xi)*(x-xi).transpose() - ConjugateFunctionType::Derivative(wc, penaltyCoeff)*Eigen::Matrix2d::Identity() ).eval();
    };

  return -2.0*RectangularQuadrature::Integrate(func, pt1,pt2);
}


template<typename ConjugateFunctionType>
Eigen::Matrix2d QuadratureDistance<ConjugateFunctionType>::LineIntegralPointHess(double                 wi,
                                Eigen::Ref<const Eigen::Vector2d> const& xi,
                                Eigen::Ref<const Eigen::Vector2d> const& xj,
                                Eigen::Vector2d const& pt1,
                                Eigen::Vector2d const& pt2,
                                double penaltyCoeff)
{
  auto func = [&](Eigen::Vector2d const& x)
    {
      return Eigen::Matrix2d(ConjugateFunctionType::Derivative(wi - (x - xi).squaredNorm(), penaltyCoeff)*(x-xj)*(x-xi).transpose());
    };

  return -2.0*LineQuadrature::Integrate(func, pt1, pt2);
}

// Explicitly instantiate the versions of this class
template class sdot::distances::QuadratureDistance<QuadraticRegularizationFunctions>;
template class sdot::distances::QuadratureDistance<GHKFunctions>;
