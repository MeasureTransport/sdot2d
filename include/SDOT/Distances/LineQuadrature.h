#ifndef DISTANCES_LINEQUADRATURE_H
#define DISTANCES_LINEQUADRATURE_H

#include <type_traits>
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

    /**
    Uses a quadrature rule to integrate a function f along a line segment
    between points pt1 and pt2.
    */
    template<typename FunctionType>
    static typename std::invoke_result<FunctionType,Eigen::Vector2d>::type Integrate(FunctionType f,
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

  };
}
}


#endif
