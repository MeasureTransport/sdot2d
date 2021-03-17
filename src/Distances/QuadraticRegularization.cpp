#include "SDOT/Distances/QuadraticRegularization.h"

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

  return QuadraticRegularization::GenericTriangularIntegral(func,pt1,pt2,pt3);
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

  return QuadraticRegularization::GenericRectangularIntegral(func,bottomLeft,topRight);
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

  return QuadraticRegularization::GenericTriangularIntegral(func,pt1,pt2,pt3);
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

  return QuadraticRegularization::GenericRectangularIntegral(func,bottomLeft,topRight);
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

  return QuadraticRegularization::GenericLineIntegral(func, pt1,pt2);
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

  return QuadraticRegularization::GenericTriangularIntegral(func,pt1,pt2,pt3);
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

  return QuadraticRegularization::GenericRectangularIntegral(func,lowerLeft,upperRight);
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

  return QuadraticRegularization::GenericTriangularIntegral(func, pt1,pt2,pt3);
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

 return QuadraticRegularization::GenericRectangularIntegral(func, pt1,pt2);
}
