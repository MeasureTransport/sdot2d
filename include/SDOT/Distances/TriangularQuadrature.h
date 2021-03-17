#ifndef DISTANCES_TRIANGULARQUADRATURE_H
#define DISTANCES_TRIANGULARQUADRATURE_H

#include <type_traits>
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

    /** Uses quadrature to compute the intgral of a function f over a triangle
    specified by corner points pt1, pt2, pt3 in counter-clockwise order.
    */
    template<typename FunctionType>
    static typename std::invoke_result<FunctionType,Eigen::Vector2d>::type Integrate(FunctionType f,
                                            Eigen::Vector2d const& pt1,
                                            Eigen::Vector2d const& pt2,
                                            Eigen::Vector2d const& pt3)
    {
        Eigen::Matrix2d A(2,2); // 2x2 matrix transforming reference coordinates to spatial coordinates

        A << pt2(0)-pt1(0), pt3(0)-pt1(0),
             pt2(1)-pt1(1), pt3(1)-pt1(1);

        double jacDet =  (pt2(0)-pt1(0))*(pt3(1)-pt1(1)) - (pt3(0)-pt1(0))*(pt2(1)-pt1(1));

        Eigen::Matrix2Xd quadPts;
        Eigen::VectorXd quadWts;
        std::tie(quadPts, quadWts) = TriangularQuadrature::Get(7);

        typename std::invoke_result<FunctionType,Eigen::Vector2d>::type output = quadWts(0)*f(pt1 + A*quadPts.col(0));

        for(int i=1; i<quadWts.size(); ++i){
          output += quadWts(i)*f(pt1 + A*quadPts.col(i)); // map pt in reference triangle to real coordinates
        }

        output *= jacDet;
        return output;
    };


  };
}
}


#endif
