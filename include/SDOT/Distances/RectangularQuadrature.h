#ifndef DISTANCES_RECTANGULARQUADRATURE_H
#define DISTANCES_RECTANGULARQUADRATURE_H

#include <type_traits>
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

    /**
    Integrates a function f over an axis-aligned rectangle with bottom left corner pt1 and top right corner pt2.
    */
    template<typename FunctionType>
    static typename std::invoke_result<FunctionType,Eigen::Vector2d>::type Integrate(FunctionType f,
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

  };
}
}


#endif
