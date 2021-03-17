#ifndef DISTANCES_QUADRATICREGULARIZATION_H
#define DISTANCES_QUADRATICREGULARIZATION_H

#include <Eigen/Core>

namespace sdot {

namespace distances{


/** @class QuadraticRegularization
    @ingroup Distances
    @details The Quadratic Regularization described in Example 2.8 of
    [Bourne et al., 2018](https://arxiv.org/pdf/1808.01962.pdf) uses the function
    \f[
    F(s) = (s-1)^2 ,
    \f]
    which results in a conjugate function \f$F^\ast(z)\f$ of the form
    \f[
    F^\ast(z) = \left\{ \begin{array}{ll} \frac{z^2}{4} + z & \text{if} z\geq -2\\ -1 & \text{otherwise}\endd{array}\right.
    \f]
*/
class QuadraticRegularization{
public:

  /** Evaluates the conjugate function \f$F^\ast(z)\f$ at a point \f$z\f$. */
  static double Evaluate(double z, double penaltyCoeff=1);

  /** Evaluates the derivative of the conjugate function \f$\frac{\partial}{\partial z}F^\ast(z)\f$ at a point \f$z\f$. */
  static double Derivative(double z, double penaltyCoeff=1);

  /** Evaluates the second derivative of the conjugate function \f$\frac{\partial}{\partial z}F^\ast(z)\f$ at a point \f$z\f$. */
  static double Derivative2(double z, double penaltyCoeff=1);


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
   \int_{\Omega} (F^\ast)^\prime(w_i - c(x,x_i)) (x-x_i) dx
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
