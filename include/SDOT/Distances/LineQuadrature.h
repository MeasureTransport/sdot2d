#ifndef DISTANCES_LINEQUADRATURE_H
#define DISTANCES_LINEQUADRATURE_H

#include <Eigen/Core>

namespace sdot{
namespace distances{

  /**
  @defgroup Quadrature
  @ingroup Distances
  @details Contains quadrature rules for integrating over lines, triangles, and rectangles.
  */

  /**
  @class LineQuadrature
  @ingroup Quadrature
  Contains functions that return scalar quadrature rules for the unit line [0,1]

  Uses Gauss-Legendre quadrature with a maximum degree of 7.

  */
  class LineQuadrature{

  public:

    static std::pair<Eigen::VectorXd, Eigen::VectorXd> Get(unsigned int degree);

    static std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre0();
    static std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre1();
    static std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre2();
    static std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre3();
    static std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre4();
    static std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre5();
    static std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre6();

  };
}
}


#endif
