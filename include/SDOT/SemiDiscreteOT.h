#ifndef SEMIDISCRETEOT_H_
#define SEMIDISCRETEOT_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "SDOT/Distribution2d.h"
#include "SDOT/LaguerreDiagram.h"
#include "SDOT/PolygonRasterize.h"

#include "SDOT/OptionUtilities.h"

namespace sdot{

  /** Class for solving semi-discrete optimal transport problems. */
  template<typename ConjugateFunctionType>
  class SemidiscreteOT
  {

  public:

    /**
      @param[in] gridIn Regular 2d grid defining cells of "continuous" distribution.
      @param[in] gridProbs Matrix of probabilities for each cell.  Rows correspond to X locations, columns to Y locations.
      @param[in] discrPts Matrix of discrete points.
      @param[in] discrProbs Probabilities assigned to each discrete point.
      @param[in] unbalancedPenalty The weight on the penalties in the unbalanced setting
    */
    SemidiscreteOT(std::shared_ptr<Distribution2d> const& distIn,
                   Eigen::Matrix2Xd                const& discrPtsIn,
                   Eigen::VectorXd                 const& discrProbsIn,
                   double                                 unbalancedPenaltyIn=1);

    /** Returns the Semi discrete OT objective, gradient, and Hessian.
        @param[in] prices A vector containing the current prices for each discrete point
        @return A tuple containing three things: The objective function, the gradient, and the Hessian.
    */
    std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>> Objective(Eigen::VectorXd const& prices) const;

    /** In Xin et al., "Centroidal Power Diagrams with Capacity Constraints: Computation, Applications, and Extension",
        the authors wish to minimize the maximum of the Kantorovich dual with
        respect to the point locations.  We will call this the centroidal objective
        because it was used in the Xin et al. paper to construct centroidal
        power diagrams.    This function, PointGradient(), computes the gradient
        of the centroidal objective with respect to the point locations.  The
        form of the gradient is given in (9) of Xin et al:
        \f[
        \nabla_{x_i}  F(X,W^{\ast}(X)) = \int_{GD(W^{|ast})_i} \nabla_{x_i} d(x,x_i) dx.
        \f]
        In our case, \f$d(x,x_i)\f$ = \frac{1}{2} (x-x_i)^2, so
        \f[
        \nabla_{x_i} d(x,x_i) = x_i-x.
        \f]
        Thus,
        \fbegin{eqnarray*}
        \int_{GD(W^{|ast})_i} \nabla_{x_i} d(x,x_i) dx &=& x_i \int_{GD(W^{|ast})_i}  dx - \int_{GD(W^{|ast})_i} x dx\\
        &=& A_i ( x_i - x_{ci}),
        \fend{eqnarray*}
        where \f$A_i\f$ is the area of the \f$i^{th}\f$ Laguerre cell \f$i\f$ and
        \f$x_{ci}\f$ is the centroid of the cell.

        @param[in] The optimal Laguerre diagram we wish to use to compute the gradient wrt each x_i
        @return A matrix of derivatives wrt to each seed point in the Laguerre diagram.  output(0,i) contains d/dx_i and output(1,i) contains d/dy_i
    */
    Eigen::Matrix2Xd PointGradient(Eigen::VectorXd const& prices, LaguerreDiagram const& lagDiag) const;

    Eigen::Matrix2Xd PointGradient() const;

    /** Sets the weight on the marginal discrepancy penalties in the unbalanced
        setting.  Has no impact on  the usual Wasserstein -2 balanced  case.
    */
    void SetPenalty(double newPenalty){unbalancedPenalty=newPenalty;};

    /** Gets the weight on the marginal discrepancy penalties in the unbalanced
        setting.
    */
    double GetPenalty() const{return unbalancedPenalty;};

    /** Solves the dual OT problem for the prices.  Uses a Trust region newton method.
    */
    std::pair<Eigen::VectorXd, double> Solve(Eigen::VectorXd const& prices0,
                                             OptionList             opts = OptionList());

    /** Returns the Laguerre diagram that was constructed during Solve.  If the
    diagram hasn't  been constructed yet, the returned shared_ptr will be a nullptr.
    */
    std::shared_ptr<LaguerreDiagram> const& Diagram() const{return lagDiag;}

    /** Resets the locations of the discrete points. */
    void SetPoints(Eigen::Matrix2Xd const& newPts);

    static std::shared_ptr<LaguerreDiagram> BuildCentroidal(std::shared_ptr<Distribution2d> const& distIn,
                                                            Eigen::Matrix2Xd                const& initialPoints,
                                                            Eigen::VectorXd                 const& pointProbs,
                                                            OptionList                             opts = OptionList());

    static std::shared_ptr<LaguerreDiagram> BuildCentroidal(std::shared_ptr<Distribution2d> const& distIn,
                                                            Eigen::VectorXd                 const& pointProbs,
                                                            OptionList                             opts = OptionList());

    static std::shared_ptr<LaguerreDiagram> BuildCentroidal(std::shared_ptr<Distribution2d> const& distIn,
                                                            unsigned int                           numPts,
                                                            OptionList                             opts = OptionList());

    /** Given an existing Laguerre diagram, this function computes the objective
        and gradient by integrating over each Laguerre cell.
        @param[in] prices The current prices vector.  These must be the same
                   prices used to construct the Laguerre diagram.
        @param[in] lagDiag Laguerre diagram computed with the current prices
        @return Pair containing the objective and gradient of the objective.
    */
    std::pair<double, Eigen::VectorXd> ComputeGradient(Eigen::VectorXd const& prices,
                                                       LaguerreDiagram const& lagDiag) const;

    /** Computes the weighted centroids, e.g.,  of each Laguerre cell. */
    Eigen::Matrix2Xd WeightedCentroids(LaguerreDiagram const& lagDiag) const;

    /** Given an existing Laguerre diagram, this function computes the Hessian by
        integrating over each pair of Laguerre cell boundaries.

    */
    Eigen::SparseMatrix<double> ComputeHessian(Eigen::VectorXd const& prices,
                                               LaguerreDiagram const& lagDiag) const;
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

    /** Penalty on the marginal discrepancies for the unbalanced formulations. */
    double unbalancedPenalty=1;

    Eigen::Matrix2Xd domain; // <- Points defining a polygon surrounding the domain of interest

    // The most recently constructed Laguerre diagram
    std::shared_ptr<LaguerreDiagram> lagDiag;

    Eigen::VectorXd optPrices;


    /**
    Computes the integral
    \f[
    \int_{\partial \Omega} (F^\ast)^\prime(w_i-c(x,x_i)) u(x) dS
    \f]
    for a line segment \f$\partial \Omega\f$ between two cells in the Laguerre
    diagram.
    */
    double LineIntegral(double wi,
                        Eigen::Ref<const Eigen::Vector2d> const& xi,
                        LaguerreDiagram::Point_2 const& srcPt,
                        LaguerreDiagram::Point_2 const& tgtPt) const;

    /**
    Computes the integral
    \f[
    \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}} 1/2 (x-x_j)^2 + 1/2 (y-y_j)^2 dydx
    \f]
    over a rectangular region \f$\Omega= [x_{min},x_{max}] \times [y_{min},y_{max}]\f$.
    */
    static double SquareIntegral(double xmin, double xmax,
                                 double ymin, double ymax,
                                 double px,   double py);

    /**
    Computes the integral
    \f[
    \int_{\Omega} 1/2 (x-x_j)^2 + 1/2 (y-y_j)^2 dxdy
    \f]
    over a triangle defined by points (x1,y1), (x2,y2), and (x3,y3).
    */
    static double TriangleIntegral(double x1, double y1,
                                   double x2, double y2,
                                   double x3, double y3,
                                   double px, double py);

  };

} // namespace sdot

#endif // #ifndef SEMIDISCRETEOT_H_
