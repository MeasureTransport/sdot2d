#ifndef BALANCEDCONJUGATE_H
#define BALANCEDCONJUGATE_H

#include <Eigen/Core>

namespace sdot {

/** @class BalancedConjugate
    @ingroup ConjugateFunctions
    @brief Defines the conjugate function \f$F^\ast\f$ for the standard balanced Wasserstein-2 case.
*/
class BalancedConjugate{
public:

  static double Evaluate(double z);

  static double Derivative(double z);

  static double TriangularIntegral(double wi,
                                   Eigen::Vector2d const& xi,
                                   Eigen::Vector2d const& pt1,
                                   Eigen::Vector2d const& pt2,
                                   Eigen::Vector2d const& pt3);

  static double RectangularIntegral(double wi,
                                   Eigen::Vector2d const& xi,
                                   Eigen::Vector2d const& lowerLeft,
                                   Eigen::Vector2d const& upperRight);

  static double TriangularIntegralDeriv(double wi,
                                        Eigen::Vector2d const& xi,
                                        Eigen::Vector2d const& pt1,
                                        Eigen::Vector2d const& pt2,
                                        Eigen::Vector2d const& pt3);

  static double RectangularIntegralDeriv(double wi,
                                         Eigen::Vector2d const& xi,
                                         Eigen::Vector2d const& lowerLeft,
                                         Eigen::Vector2d const& upperRight);

  static double TriangularIntegralDeriv2(double wi,
                                         Eigen::Vector2d const& xi,
                                         Eigen::Vector2d const& pt1,
                                         Eigen::Vector2d const& pt2,
                                         Eigen::Vector2d const& pt3);

  static double RectangularIntegralDeriv2(double wi,
                                          Eigen::Vector2d const& xi,
                                          Eigen::Vector2d const& lowerLeft,
                                          Eigen::Vector2d const& upperRight);

}; // class BalancedConjugate

} // namespace sdot

#endif
