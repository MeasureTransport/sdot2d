#ifndef DISTANCES_WASSERSTEIN2_H
#define DISTANCES_WASSERSTEIN2_H

#include <Eigen/Core>

namespace sdot {

namespace distances{


/** @class Wasserstein2
    @ingroup Distances
    @brief Defines the conjugate function \f$F^\ast\f$ for the standard balanced Wasserstein-2 case.
*/
class Wasserstein2{
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



}; // class Wasserstein2

}// namespace distances
} // namespace sdot

#endif
