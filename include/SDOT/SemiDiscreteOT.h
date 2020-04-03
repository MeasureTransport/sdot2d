#ifndef SEMIDISCRETEOT_H_
#define SEMIDISCRETEOT_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "SDOT/Distribution2d.h"
#include "SDOT/LaguerreDiagram.h"
#include "SDOT/PolygonRasterize.h"

namespace sdot{

  /** Class for solving semi-discrete optimal transport problems. */
  class SemidiscreteOT
  {

  public:

    /**
      @param[in] gridIn Regular 2d grid defining cells of "continuous" distribution.
      @param[in] gridProbs Matrix of probabilities for each cell.  Rows correspond to X locations, columns to Y locations.
      @param[in] discrPts Matrix of discrete points.
      @param[in] discrProbs Probabilities assigned to each discrete point.
    */
    SemidiscreteOT(std::shared_ptr<Distribution2d> const& distIn,
                   Eigen::Matrix2Xd                const& discrPtsIn,
                   Eigen::VectorXd                 const& discrProbsIn);

    /** Returns the Semi discrete OT objective, gradient, and Hessian.
        @param[in] prices A vector containing the current prices for each discrete point
        @return A tuple containing three things: The objective function, the gradient, and the Hessian.
    */
    std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>> Objective(Eigen::VectorXd const& prices) const;


    /** Solves the dual OT problem for the prices.  Uses a Trust region newton method.
    */
    std::pair<Eigen::VectorXd, double> Solve(Eigen::VectorXd const& prices0);

    /** Returns the Laguerre diagram that was constructed during Solve.  If the
    diagram hasn't  been constructed yet, the returned shared_ptr will be a nullptr.
    */
    std::shared_ptr<LaguerreDiagram> const& Diagram() const{return lagDiag;}

  private:

    /** Solves the trust region subproblem using a Steihaug-CG approach.
        @param[in] obj The current value of the objective.
        @param[in] grad A view of the gradient with  one of the components removed to make the solution unique.
        @param[in] hess A view of the Hessian with one component  removed.
        @param[in] trustRadius The current trust region radius.
        @return A vector holding the optimization step.
    */
    Eigen::VectorXd SolveSubProblem(double obj,
                                    Eigen::Ref<const Eigen::VectorXd> const& grad,
                                    Eigen::Ref<const Eigen::SparseMatrix<double>> const& hess,
                                    double trustRadius) const;


    std::shared_ptr<Distribution2d> dist;
    std::shared_ptr<RegularGrid> grid;

    Eigen::Matrix2Xd discrPts;
    Eigen::VectorXd  discrProbs;

    Eigen::Matrix2Xd domain; // <- Points defining a polygon surrounding the domain of interest

    // The most recently constructed Laguerre diagram
    std::shared_ptr<LaguerreDiagram> lagDiag;

    /** Given an existing Laguerre diagram, this function computes the objective
        and gradient by integrating over each Laguerre cell.
        @param[in] prices The current prices vector.  These must be the same
                   prices used to construct the Laguerre diagram.
        @param[in] lagDiag Laguerre diagram computed with the current prices
        @return Pair containing the objective and gradient of the objective.
    */
    std::pair<double, Eigen::VectorXd> ComputeGradient(Eigen::VectorXd const& prices,
                                                       LaguerreDiagram const& lagDiag) const;

    /** Given an existing Laguerre diagram, this function computes the Hessian by
        integrating over each pair of Laguerre cell boundaries.

    */
    Eigen::SparseMatrix<double> ComputeHessian(LaguerreDiagram const& lagDiag) const;

    double LineIntegral(LaguerreDiagram::Point_2 const& srcPt,
                        LaguerreDiagram::Point_2 const& tgtPt) const;

    static double SquareIntegral(double xmin, double xmax,
                                 double ymin, double ymax,
                                 double px,   double py);

    static double TriangleIntegral(double x1, double y1,
                                   double x2, double y2,
                                   double x3, double y3,
                                   double px, double py);

  };

} // namespace sdot

#endif // #ifndef SEMIDISCRETEOT_H_
