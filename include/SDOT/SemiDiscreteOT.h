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
    Eigen::Matrix2Xd PointGradient(Eigen::VectorXd const& prices,
                                   LaguerreDiagram const& lagDiag) const;

    Eigen::Matrix2Xd PointGradient() const;

    /**
    Computes the Hessian of the Kantorovich dual objective wrt to the point
    locations.   The components are order (x_1, y_1, x_2, y_2, ..., x_N, y_N),
    which corresponds to the column-stacked version matrix returned by PointGradient().
    */
    Eigen::SparseMatrix<double> PointHessian(Eigen::VectorXd const& prices,
                                             LaguerreDiagram const& lagDiag) const;

    Eigen::SparseMatrix<double> PointHessian() const;

    Eigen::Matrix2Xd LloydPointHessian(Eigen::VectorXd const& prices,
                                       LaguerreDiagram const& lagDiag) const;

    Eigen::Matrix2Xd LloydPointHessian() const;


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

    void CheckNormalization(){};

    template<typename FunctionType>
    typename std::invoke_result<FunctionType,Eigen::Vector2d,Eigen::Vector2d>::type LineIntegral(FunctionType f,
                                    LaguerreDiagram::Point_2 const& srcPt,
                                    LaguerreDiagram::Point_2 const& tgtPt) const
    {
      const double compTol = 1e-11;

      double xs = CGAL::to_double(srcPt.x());
      double ys = CGAL::to_double(srcPt.y());
      double xt = CGAL::to_double(tgtPt.x());
      double yt = CGAL::to_double(tgtPt.y());

      // Points on either end of a line segment
      Eigen::Vector2d startPt(2), endPt(2);

      // if the edge is vertical...
      if(std::abs(xt-xs)<compTol){

        startPt(0) = xs;
        endPt(0) = xs;

        // Assume we are working from bottom to top
        if(yt<ys){
          std::swap(xt,xs);
          std::swap(yt,ys);
        }

        unsigned int xInd = grid->LeftNode(xs);
        unsigned int yInd = grid->BottomNode(ys);
        const double maxY = CGAL::to_double(yt);

        double currY = CGAL::to_double( ys );
        double nextY = grid->TopNode(ys);

        // If the source point is on a boundary...
        if(std::abs(nextY-currY)<compTol)
          nextY += grid->dy;

        // Initial value
        startPt(1) = currY;
        endPt(1) = nextY;
        typename std::invoke_result<FunctionType,Eigen::Vector2d,Eigen::Vector2d>::type val = f(startPt, endPt)*dist->Density(xInd,yInd);

        yInd++;
        currY = nextY;
        nextY = currY+grid->dy;

        while(nextY<maxY-compTol){
          startPt(1) = currY;
          endPt(1) = nextY;

          val += f(startPt, endPt)*dist->Density(xInd,yInd);

          yInd++;
          currY = nextY;
          nextY = currY+grid->dy;
        }

        startPt(1) = currY;
        endPt(1) = maxY;

        val += f(startPt, endPt)*dist->Density(xInd,yInd);

        return val;

      // If the edge is horizontal
      }else if(std::abs(yt-ys)<compTol){

        startPt(1) = ys;
        endPt(1) = ys;

        // Assume we are working from left to right and swap direction if we're not
        if(xt<xs){
          std::swap(xt,xs);
          std::swap(yt,ys);
        }

        unsigned int xInd = grid->LeftNode(xs);
        unsigned int yInd = grid->BottomNode(ys);

        const double maxX = CGAL::to_double(xt);

        double currX = CGAL::to_double( xs );
        double nextX = grid->RightNode(xs);

        // If the source is on a boundary...
        if(std::abs(nextX-currX)<compTol)
          nextX += grid->dx;

        // Initial value
        startPt(0) = currX;
        endPt(0) = nextX;
        typename std::invoke_result<FunctionType,Eigen::Vector2d,Eigen::Vector2d>::type val = f(startPt, endPt)*dist->Density(xInd,yInd);
        xInd++;
        currX = nextX;
        nextX = currX+grid->dx;

        while(nextX<maxX-compTol){
          startPt(0) = currX;
          endPt(0) = nextX;

          val += f(startPt, endPt)*dist->Density(xInd,yInd);

          xInd++;
          currX = nextX;
          nextX = currX+grid->dx;
        }

        startPt(0) = currX;
        endPt(0) = maxX;

        val += f(startPt, endPt)*dist->Density(xInd,yInd);

        return val;

      // Otherwise there is a nonzero finite slope
      }else{

        // Assume we are working from left to right and swap direction if we're not
        if(xt<xs){
          std::swap(xt,xs);
          std::swap(yt,ys);
        }

        double dy = yt-ys;
        double dx = xt-xs;

        // The length of the source -> target line segment
        double segLenth = std::sqrt(dy*dy+dx*dx);

        /* We parameterize the line segment as ps + t*(pt-ps), where ps is the
           source point and target is the target node.  As we walk along the line
           segment and compute the integral,
           - currt holds the current position along the line
           - nextt_vert holds the next value of t where the line segment crosses a vertical grid line
           - nextt_horz holds the next value of t where the line segment crosses a horizontal grid line
           - nextt holds the minimum of nextt_vert and nextt_horz
        */
        double currt = 0.0;
        double nextt_vert, nextt_horz;
        double nextt;


        // Compute the slope of the line
        bool posSlope = dy>0;

        // Get the starting grid cells
        unsigned int xInd = grid->LeftNode(xs);
        unsigned int yInd = grid->BottomNode(ys);

        // Handle situations where the source starts on a boundary
        if(std::abs(yInd*grid->dy+grid->yMin - ys)<compTol){
          if(yt<ys-compTol){
            yInd--;
          }
        }

        nextt_vert = std::min(1.0, ( (xInd+1)*grid->dx + grid->xMin - xs) / dx);

        if(posSlope){
          nextt_horz = std::min(1.0, ( (yInd+1.0)*grid->dy + grid->yMin - ys) / dy);
        }else{
          nextt_horz = std::min(1.0, ( yInd*grid->dy + grid->yMin - ys) / dy);
        }

        nextt = std::min(nextt_horz,nextt_vert);

        // initial value
        startPt(0) = xs + dx*currt;
        startPt(1) = ys + dy*currt;
        endPt(0) = xs + dx*nextt;
        endPt(1) = ys + dy*nextt;
        typename std::invoke_result<FunctionType,Eigen::Vector2d,Eigen::Vector2d>::type val = f(startPt, endPt)*dist->Density(xInd,yInd);
        val -= val; // Make surethe value starts off at 0

        while(nextt<1.0-compTol){
          startPt(0) = xs + dx*currt;
          startPt(1) = ys + dy*currt;
          endPt(0) = xs + dx*nextt;
          endPt(1) = ys + dy*nextt;

          val += f(startPt, endPt)*dist->Density(xInd,yInd);

          // we leave out the top or bottom
          if(std::abs(nextt-nextt_horz)<compTol){

            if(posSlope){
              yInd++;
              nextt_horz = std::min(1.0, ((yInd+1)*grid->dy + grid->yMin-ys)/dy);
            }else{
              yInd--;
              nextt_horz = std::min(1.0, ((yInd)*grid->dy + grid->yMin-ys)/dy);
            }
          }

          // leave out the right (note that we could leave out the right and top/bottom if we leave a corner)
          if(std::abs(nextt-nextt_vert)<compTol){
            xInd++;
            nextt_vert = std::min(1.0,  ( (xInd+1)*grid->dx + grid->xMin - xs) / dx);
          }
          currt = nextt;
          nextt = std::min(nextt_horz,nextt_vert);
        }

        if((xInd<grid->NumCells(0))&&(yInd<grid->NumCells(1))){
          startPt(0) = xs + dx*currt;
          startPt(1) = ys + dy*currt;
          endPt(0) = xs + dx*nextt;
          endPt(1) = ys + dy*nextt;

          val += f(startPt, endPt)*dist->Density(xInd,yInd);
        }

        return val;
      }
    }

  };

} // namespace sdot

#endif // #ifndef SEMIDISCRETEOT_H_
