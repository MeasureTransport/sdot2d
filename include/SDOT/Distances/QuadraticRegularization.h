#ifndef DISTANCES_QUADRATICREGULARIZATION_H
#define DISTANCES_QUADRATICREGULARIZATION_H

#include <type_traits>

#include "SDOT/Distances/TriangularQuadrature.h"
#include "SDOT/Distances/RectangularQuadrature.h"
#include "SDOT/Distances/LineQuadrature.h"

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

  /** Uses quadrature to compute an integral of \f$ f(w_i-c(x,x_i))\f$ over a triangle. */
  template<typename FunctionType>
  static typename std::invoke_result<FunctionType,Eigen::Vector2d>::type GenericTriangularIntegral(FunctionType f,
                                          Eigen::Vector2d const& pt1,
                                          Eigen::Vector2d const& pt2,
                                          Eigen::Vector2d const& pt3)
  {
      Eigen::Matrix2d A(2,2); // 2x2 matrix transforming reference coordinates to spatial coordinates

      A << pt2(0)-pt1(0), pt3(0)-pt1(0),
           pt2(1)-pt1(1), pt3(1)-pt1(1);

      double jacDet =  (pt2(0)-pt1(0))*(pt3(1)-pt1(1)) - (pt3(0)-pt1(0))*(pt2(1)-pt1(1));

      Eigen::Matrix2Xd quadPts;
      Eigen::VectorXd quadWts;
      std::tie(quadPts, quadWts) = TriangularQuadrature::Get(7);

      typename std::invoke_result<FunctionType,Eigen::Vector2d>::type output = quadWts(0)*f(pt1 + A*quadPts.col(0));

      for(int i=1; i<quadWts.size(); ++i){
        output += quadWts(i)*f(pt1 + A*quadPts.col(i)); // map pt in reference triangle to real coordinates
      }

      output *= jacDet;
      return output;
  };

  /** Uses quadrature to compute an integral of \f$ f(w_i-c(x,x_i))\f$ over a line segment. */
  template<typename FunctionType>
  static typename std::invoke_result<FunctionType,Eigen::Vector2d>::type GenericLineIntegral(FunctionType f,
                                    Eigen::Vector2d const& pt1,
                                    Eigen::Vector2d const& pt2)
  {
      Eigen::Vector2d diff = pt2-pt1;
      double segLength = diff.norm();


      Eigen::VectorXd quadPts, quadWts;
      std::tie(quadPts, quadWts) = LineQuadrature::Get(7);

      typename std::invoke_result<FunctionType,Eigen::Vector2d>::type output = quadWts(0)*f(pt1 + diff*quadPts(0));

      for(int i=1; i<quadWts.size(); ++i){
        output += quadWts(i)*f(pt1 + diff*quadPts(i)); // map pt in reference triangle to real coordinates
      }

      output *= segLength;
      return output;
  };

  /** Uses quadrature to compute an integral of \f$ f(w_i-c(x,x_i))\f$ over a triangle. */
  template<typename FunctionType>
  static typename std::invoke_result<FunctionType,Eigen::Vector2d>::type GenericRectangularIntegral(FunctionType f,
                                         Eigen::Vector2d const& pt1,
                                         Eigen::Vector2d const& pt2)
  {

      Eigen::Matrix2Xd quadPts;
      Eigen::VectorXd quadWts;
      std::tie(quadPts, quadWts) = RectangularQuadrature::Get(7);

      Eigen::Vector2d scale = pt2-pt1;

      typename std::invoke_result<FunctionType,Eigen::Vector2d>::type output = quadWts(0)*f(pt1 + scale.asDiagonal()*quadPts.col(0));


      for(int i=1; i<quadWts.size(); ++i){
        output += quadWts(i)*f(pt1 + scale.asDiagonal()*quadPts.col(i)); // map pt in reference triangle to real coordinates
      }

      output *= scale.prod();
      return output;
  };

}; // class Wasserstein2

}// namespace distances
} // namespace sdot

#endif
