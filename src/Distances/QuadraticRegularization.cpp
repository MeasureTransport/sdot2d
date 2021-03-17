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
  return QuadraticRegularization::GenericTriangularIntegral(QuadraticRegularization::Evaluate, wi,xi,pt1,pt2,pt3, penaltyCoeff);
}

double QuadraticRegularization::RectangularIntegral(double wi,
                                         Eigen::Ref<const Eigen::Vector2d> const& xi,
                                         Eigen::Vector2d const& bottomLeft,
                                         Eigen::Vector2d const& topRight,
                                         double penaltyCoeff)
{
  return QuadraticRegularization::GenericRectangularIntegral(QuadraticRegularization::Evaluate, wi,xi,bottomLeft,topRight, penaltyCoeff);
}

double QuadraticRegularization::TriangularIntegralDeriv(double wi,
                                             Eigen::Ref<const Eigen::Vector2d> const& xi,
                                             Eigen::Vector2d const& pt1,
                                             Eigen::Vector2d const& pt2,
                                             Eigen::Vector2d const& pt3,
                                             double penaltyCoeff)
{
  return QuadraticRegularization::GenericTriangularIntegral(QuadraticRegularization::Derivative, wi,xi,pt1,pt2,pt3, penaltyCoeff);
}

double QuadraticRegularization::RectangularIntegralDeriv(double wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& bottomLeft,
                                              Eigen::Vector2d const& topRight,
                                              double penaltyCoeff)
{
  return QuadraticRegularization::GenericRectangularIntegral(QuadraticRegularization::Derivative, wi,xi,bottomLeft,topRight, penaltyCoeff);
}

double QuadraticRegularization::LineIntegralDeriv(double                 wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  double penaltyCoeff)
{
  return QuadraticRegularization::GenericLineIntegral(QuadraticRegularization::Derivative, wi,xi,pt1,pt2, penaltyCoeff);
}


double QuadraticRegularization::TriangularIntegralDeriv2(double                 wi,
                                              Eigen::Ref<const Eigen::Vector2d> const& xi,
                                              Eigen::Vector2d const& pt1,
                                              Eigen::Vector2d const& pt2,
                                              Eigen::Vector2d const& pt3,
                                              double penaltyCoeff)
{
  return QuadraticRegularization::GenericTriangularIntegral(QuadraticRegularization::Derivative2, wi,xi,pt1,pt2,pt3, penaltyCoeff);
}

double QuadraticRegularization::RectangularIntegralDeriv2(double                 wi,
                                               Eigen::Ref<const Eigen::Vector2d> const& xi,
                                               Eigen::Vector2d const& lowerLeft,
                                               Eigen::Vector2d const& upperRight,
                                               double penaltyCoeff)
{
  return QuadraticRegularization::GenericRectangularIntegral(QuadraticRegularization::Derivative2, wi,xi,lowerLeft,upperRight, penaltyCoeff);
}
