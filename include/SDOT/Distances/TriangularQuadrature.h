#ifndef DISTANCES_TRIANGULARQUADRATURE_H
#define DISTANCES_TRIANGULARQUADRATURE_H

#include <Eigen/Core>

namespace sdot{
namespace distances{

  /**
  @class TriangularQuadrature
  @ingroup Quadrature
  
  Contains functions that return quadrature rules for the unit triangle [(0,0), (1,0), (1,1)]

  Maximum degree is 7.

  This adapts the "Strang" rules reported at https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html


  */
  class TriangularQuadrature{

  public:

    static std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> Get(unsigned int degree);

    static std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> Centroid();
    static std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> Strang1();
    static std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> Strang3();
    static std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> Strang5();
    static std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> Strang7();
    static std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> Strang9();
    static std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> Strang10();

  };
}
}


#endif
