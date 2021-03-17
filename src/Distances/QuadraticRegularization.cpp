#include "SDOT/Distances/QuadraticRegularization.h"

using namespace sdot;
using namespace sdot::distances;


double QuadraticRegularization::Evaluate(double z){
  if(z > -2){
    return 0.25*z*z+ z;
  }else{
    return -1.0;
  }
};

double QuadraticRegularization::Derivative(double z){
  if(z > -2){
    return 0.5*z + 1.0;
  }else{
    return 0;
  }
};

double QuadraticRegularization::Derivative2(double z){
  if(z > -2){
    return 0.5;
  }else{
    return 0;
  }
};

double QuadraticRegularization::TriangularIntegral(double wi,
                                        Eigen::Ref<const Eigen::Vector2d> const& xi,
                                        Eigen::Vector2d const& pt1,
                                        Eigen::Vector2d const& pt2,
                                        Eigen::Vector2d const& pt3)
{
  return QuadraticRegularization::GenericTriangularIntegral(QuadraticRegularization::Evaluate, wi,xi,pt1,pt2,pt3);
}

double QuadraticRegularization::RectangularIntegral(double wi,
                                         Eigen::Ref<const Eigen::Vector2d> const& xi,
                                         Eigen::Vector2d const& bottomLeft,
                                         Eigen::Vector2d const& topRight)
{
  return QuadraticRegularization::GenericRectangularIntegral(QuadraticRegularization::Evaluate, wi,xi,bottomLeft,topRight);
}

double QuadraticRegularization::TriangularIntegralDeriv(double wi,
                                             Eigen::Ref<const Eigen::Vector2d> const& xi,
                                             Eigen::Vector2d const& pt1,
                                             Eigen::Vector2d const& pt2,
                                             Eigen::Vector2d const& pt3)
{
  return QuadraticRegularization::GenericTriangularIntegral(QuadraticRegularization::Derivative, wi,xi,pt1,pt2,pt3);
}

double QuadraticRegularization::RectangularIntegralDeriv(double wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& bottomLeft,
                                              Eigen::Vector2d const& topRight)
{
  return QuadraticRegularization::GenericRectangularIntegral(QuadraticRegularization::Derivative, wi,xi,bottomLeft,topRight);
}

double QuadraticRegularization::LineIntegralDeriv(double                 wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2)
{
  return QuadraticRegularization::GenericLineIntegral(QuadraticRegularization::Derivative, wi,xi,pt1,pt2);
}


double QuadraticRegularization::TriangularIntegralDeriv2(double                 wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& pt1,
                                              Eigen::Vector2d const& pt2,
                                              Eigen::Vector2d const& pt3)
{
  return QuadraticRegularization::GenericTriangularIntegral(QuadraticRegularization::Derivative2, wi,xi,pt1,pt2,pt3);
}

double QuadraticRegularization::RectangularIntegralDeriv2(double                 wi,
                                               Eigen::Ref<const Eigen::Vector2d> const& xi,
                                               Eigen::Vector2d const& lowerLeft,
                                               Eigen::Vector2d const& upperRight)
{
  return QuadraticRegularization::GenericRectangularIntegral(QuadraticRegularization::Derivative2, wi,xi,lowerLeft,upperRight);
}
