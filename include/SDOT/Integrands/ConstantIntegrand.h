#ifndef CONSTANTINTEGRAND_H
#define CONSTANTINTEGRAND_H

#include "SDOT/Integrands/Integrand.h"

namespace sdot{

  /** @class ConstantIntegrand
      Defines a constant integrand f(x).
  */
  class ConstantIntegrand : public Integrand {
  public:

    ConstantIntegrand(double constVal  = 1.0);

    virtual ~ConstantIntegrand() = default;

    virtual double TriangularIntegral(Eigen::Vector2d const& pt1,
                                      Eigen::Vector2d const& pt2,
                                      Eigen::Vector2d const& pt3) override;

    virtual double RectangularIntegral(Eigen::Vector2d const& bottomLeft,
                                       Eigen::Vector2d const& topRight) override;


  private:
    const double val;
  };

};


#endif
