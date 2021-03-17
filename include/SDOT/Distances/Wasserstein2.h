#ifndef DISTANCES_WASSERSTEIN2_H
#define DISTANCES_WASSERSTEIN2_H

#include <Eigen/Core>

namespace sdot {

namespace distances{


/** @class Wasserstein2
    @ingroup Distances
    @details The standard Wasserstein-2 optimal transport problem is obtained
    for a cost function \f$c(x,y)=\frac{1}{2}\|x-y\|^2\f$ and a function
    \f[
    F(s) = \left\{ \begin{array}{ll} 0 & \text{if } s=1\\ \infty & \text{otherwise}\end{array} \right. ,
    \f]
    which results in a conjugate function \f$F^\ast(z)\f$ of the simple form
    \f[
    F^\ast(z) = z
    \f]
*/
class Wasserstein2{
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


}; // class Wasserstein2

}// namespace distances
} // namespace sdot

#endif
