#ifndef INTEGRAND_H
#define INTEGRAND_H

#include <Eigen/Core>

namespace sdot{

/**
@defgroup Integrands
*/

  /** @class Integrand
      @ingroup Integrands
      Provides an interface to functions that will be integrated
      over Laguerre cells.   Let $f(x)$ be a function defined over
      a two dimensional domain $\Omega_x$.   This class provides
      an interface to that function with the intent of computing
      $$
      \int_{A} f(x) dx
      $$
      for some region $A\subseteq \Omega_x$ that is either rectangular
      or triangular.
  */
  class Integrand {
  public:

    virtual ~Integrand() = default;

    /** Returns the integral of the integrand f(x) over a triangular region.
        @param[in] pt1 The first point in the triangle.
        @param[in] pt2 The second point in the triangle.
        @param[in] pt3 The last point in the triangle.
        @return double The integral $\int_A f(x)dx$ for the triangle $A$.
    */
    virtual double TriangularIntegral(Eigen::Vector2d const& pt1,
                                      Eigen::Vector2d const& pt2,
                                      Eigen::Vector2d const& pt3) = 0;

    /** Returns the integral of the integrand f(x) over an axis-aligned
        rectangular region.
        @param[in] bottomLeft The point at the bottom left of the rectangle.
        @param[in] topRight The point at the top right of the rectangle.
        @return double The integral $\int_A f(x)dx$ for the rectangle $A$.
    */
    virtual double RectangularIntegral(Eigen::Vector2d const& bottomLeft,
                                       Eigen::Vector2d const& topRight) = 0;

  };

};


#endif
