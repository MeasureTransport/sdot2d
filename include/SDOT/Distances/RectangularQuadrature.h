#ifndef DISTANCES_RECTANGULARQUADRATURE_H
#define DISTANCES_RECTANGULARQUADRATURE_H

#include <Eigen/Core>

namespace sdot{
namespace distances{

  /**
  @class RectangularQuadrature
  @ingroup Quadrature

  Contains functions that return quadrature rules for the unit square [(0,0), (1,0), (1,0), (1,1)].

  The quadrature rules are built from tensor products of 1D Gauss-Legendre rules.

  Maximum degree is 5.
  */
  class RectangularQuadrature{

  public:

    static std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> Get(unsigned int order);

  };
}
}


#endif
