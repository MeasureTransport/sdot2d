#ifndef DISTANCES_QUADRATUREDISTANCE_H
#define DISTANCES_QUADRATUREDISTANCE_H

#include <Eigen/Core>

namespace sdot {

namespace distances{

/**
@class QuadratureDistance
@ingroup Distances
@brief Provides an interface to integrating conjugate functions with quadrature.
@details Following the formulation of [Bourne et al.], entropic unbalanced optimal
transport can be framed in terms of integrals involving a conjugate function
\f$F^\ast(z)\f$ and its derivatives.  This template class can approximate all of
the integrals necessary for SDOT using quadrature.   The template argument is
another class which must implement three static functions:
- `double Evaluate(double z, double penaltyCoeff=1)`
- `double Derivative(double z, double penaltyCoeff=1)`
- `double Derivative2(double z, double penaltyCoeff=1)`

@seealso  QuadraticRegularizationFunctions
*/
template<typename ConjugateFunctionType>
class QuadratureDistance{
public:

  static double Evaluate(double z, double penaltyCoeff=1){return  ConjugateFunctionType::Evaluate(z,penaltyCoeff);};
  static double Derivative(double z, double penaltyCoeff=1){return  ConjugateFunctionType::Derivative(z,penaltyCoeff);};
  static double Derivative2(double z, double penaltyCoeff=1){return  ConjugateFunctionType::Derivative2(z,penaltyCoeff);};

  /**
  Returns
  \f[
  \int_{\Omega_{\text{tri}}} F^\ast(w_i - c(x,x_i)) dx,
  \f]
  where \f$\Omega_{\text{tri}}\f$ is a triangular region defined by three corner
  points.

    Uses numerical quadrature to approximate the integral.
  */
  static double TriangularIntegral(double                 wi,
                                   Eigen::Ref<const Eigen::Vector2d> const& xi,
                                   Eigen::Vector2d const& pt1,
                                   Eigen::Vector2d const& pt2,
                                   Eigen::Vector2d const& pt3,
                                   double penaltyCoeff=1);

   /**
   Returns
   \f[
   \int_{\Omega_{\text{rect}}} F^\ast(w_i - c(x,x_i)) dx,
   \f]
   where \f$\Omega_{\text{rect}}\f$ is a rectangular region defined by two of
   its corner points.  The corner on the bottom left (minx, miny) and the corner
   on the top right (maxx, maxy).

     Uses numerical quadrature to approximate the integral.
   */
  static double RectangularIntegral(double                wi,
                                    Eigen::Ref<const Eigen::Vector2d> const& xi,
                                    Eigen::Vector2d const& lowerLeft,
                                    Eigen::Vector2d const& upperRight,
                                    double penaltyCoeff=1);

  static double TriangularIntegralDeriv(double                 wi,
                                        Eigen::Ref<const Eigen::Vector2d> const& xi,
                                        Eigen::Vector2d const& pt1,
                                        Eigen::Vector2d const& pt2,
                                        Eigen::Vector2d const& pt3,
                                        double penaltyCoeff=1);

  static double RectangularIntegralDeriv(double                 wi,
                                         Eigen::Ref<const Eigen::Vector2d> const& xi,
                                         Eigen::Vector2d const& lowerLeft,
                                         Eigen::Vector2d const& upperRight,
                                         double penaltyCoeff=1);

  static double TriangularIntegralDeriv2(double                 wi,
                                         Eigen::Ref<const Eigen::Vector2d> const& xi,
                                         Eigen::Vector2d const& pt1,
                                         Eigen::Vector2d const& pt2,
                                         Eigen::Vector2d const& pt3,
                                         double penaltyCoeff=1);

  static double RectangularIntegralDeriv2(double                 wi,
                                          Eigen::Ref<const Eigen::Vector2d> const& xi,
                                          Eigen::Vector2d const& lowerLeft,
                                          Eigen::Vector2d const& upperRight,
                                          double penaltyCoeff=1);

  /**
  Returns
  \f[
  \int_{\partial \Omega} (F^\ast)^\prime(w_i - c(x,x_i)) dS,
  \f]
  over a line segment $\partial \Omega$ define by a starting point and an
  ending point.   Because \f$(F^\ast)^\prime=1\f$ for the Wasserstein-2 distance,
  the value of this line integral is simply the length of the line segment.
  */
  static double LineIntegralDeriv(double                 wi,
                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                  Eigen::Vector2d const& pt1,
                                  Eigen::Vector2d const& pt2,
                                  double penaltyCoeff=1);


  static Eigen::Vector2d TriangularIntegralMarginalCentroid(double wi,
                                            Eigen::Ref<const Eigen::Vector2d> const& xi,
                                            Eigen::Vector2d const& pt1,
                                            Eigen::Vector2d const& pt2,
                                            Eigen::Vector2d const& pt3,
                                            double penaltyCoeff=1);


  /**
   This function is used to compute centroids weighted by the marginal
   distribution of the optimal coupling.  It returns the integral
   of
   \f[
   \int_{\Omega} (F^\ast)^\prime(w_i - c(x,x_i)) x dx
   \f]
   over a rectangular region \f$\Omega\f$.
  */
  static Eigen::Vector2d RectangularIntegralMarginalCentroid(double wi,
                                             Eigen::Ref<const Eigen::Vector2d> const& xi,
                                             Eigen::Vector2d const& pt1,
                                             Eigen::Vector2d const& pt2,
                                             double penaltyCoeff=1);

 static double TriangularIntegralMarginalMass(double wi,
                                           Eigen::Ref<const Eigen::Vector2d> const& xi,
                                           Eigen::Vector2d const& pt1,
                                           Eigen::Vector2d const& pt2,
                                           Eigen::Vector2d const& pt3,
                                           double penaltyCoeff=1);


 /**
  This function is used to compute centroids weighted by the marginal
  distribution of the optimal coupling.  It returns the integral
  of
  \f[
  \int_{\Omega} (F^\ast)^\prime(w_i - c(x,x_i)) x dx
  \f]
  over a rectangular region \f$\Omega\f$.
 */
 static double RectangularIntegralMarginalMass(double wi,
                                            Eigen::Ref<const Eigen::Vector2d> const& xi,
                                            Eigen::Vector2d const& pt1,
                                            Eigen::Vector2d const& pt2,
                                            double penaltyCoeff=1);

  static Eigen::Vector2d TriangularIntegralPointGrad(double wi,
                                            Eigen::Ref<const Eigen::Vector2d> const& xi,
                                            Eigen::Vector2d const& pt1,
                                            Eigen::Vector2d const& pt2,
                                            Eigen::Vector2d const& pt3,
                                            double penaltyCoeff=1);



  /**
   This function is used to compute the PointGradient.  It returns the integral
   of
   \f[
   \int_{\Omega} 2(F^\ast)^\prime(w_i - c(x,x_i)) (x-x_i) dx
   \f]
   over a rectangular region \f$\Omega\f$.
  */
  static Eigen::Vector2d RectangularIntegralPointGrad(double wi,
                                             Eigen::Ref<const Eigen::Vector2d> const& xi,
                                             Eigen::Vector2d const& pt1,
                                             Eigen::Vector2d const& pt2,
                                             double penaltyCoeff=1);


  static Eigen::Matrix2d TriangularIntegralPointHessDiag(double wi,
                                                    Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                    Eigen::Vector2d const& pt1,
                                                    Eigen::Vector2d const& pt2,
                                                    Eigen::Vector2d const& pt3,
                                                    double penaltyCoeff=1);

  static Eigen::Matrix2d RectangularIntegralPointHessDiag(double wi,
                                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                                  Eigen::Vector2d const& pt1,
                                                  Eigen::Vector2d const& pt2,
                                                  double penaltyCoeff=1);

  static Eigen::Matrix2d LineIntegralPointHess(double                 wi,
                                  Eigen::Ref<const Eigen::Vector2d> const& xi,
                                  Eigen::Ref<const Eigen::Vector2d> const& xj,
                                  Eigen::Vector2d const& pt1,
                                  Eigen::Vector2d const& pt2,
                                  double penaltyCoeff=1);

}; // class Wasserstein2


}// namespace distances
} // namespace sdot

#endif
