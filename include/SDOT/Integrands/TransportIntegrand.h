#ifndef TRANSPORTINTEGRAND_H
#define TRANSPORTINTEGRAND_H

#include "SDOT/Integrands/Integrand.h"

namespace sdot {

  /**
   @class TransportIntegrand
   @ingroup Integrands
   Defines the transport cost to move mass from a point x to p.
   */
  class TransportIntegrand : public Integrand{
  public:

    TransportIntegrand(Eigen::Vector2d const& pIn);

    virtual ~TransportIntegrand() = default;

    virtual double TriangularIntegral(Eigen::Vector2d const& pt1,
                                      Eigen::Vector2d const& pt2,
                                      Eigen::Vector2d const& pt3) override;

    virtual double RectangularIntegral(Eigen::Vector2d const& bottomLeft,
                                       Eigen::Vector2d const& topRight) override;


  private:
    const Eigen::Vector2d p;

  };

} //namespace sdot

#endif
