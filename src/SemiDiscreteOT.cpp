#include "SDOT/SemidiscreteOT.h"

#include "SDOT/Integrands/ConstantIntegrand.h"
#include "SDOT/Integrands/TransportIntegrand.h"

#include <CGAL/Kernel/global_functions.h>

#include <algorithm>
#include <chrono>

using namespace sdot;


SemidiscreteOT::SemidiscreteOT(std::shared_ptr<Distribution2d> const& distIn,
                               Eigen::Matrix2Xd                const& discrPtsIn,
                               Eigen::VectorXd                 const& discrProbsIn) : dist(distIn),
                                                                                      grid(distIn->Grid()),
                                                                                      discrPts(discrPtsIn),
                                                                                      discrProbs(discrProbsIn)
{
  assert(discrPtsIn.cols()==discrProbsIn.size());

  // Check to make sure all the points are inside the grid domain
  for(unsigned int i=0; i<discrPts.cols(); ++i){
    assert(discrPts(0,i)>=grid->xMin);
    assert(discrPts(0,i)<=grid->xMax);
    assert(discrPts(1,i)>=grid->yMin);
    assert(discrPts(1,i)<=grid->yMax);
  }
}

void SemidiscreteOT::SetPoints(Eigen::Matrix2Xd const& newPts){
  assert(newPts.cols()==discrProbs.size());

  // Check to make sure all the points are inside the grid domain
  for(unsigned int i=0; i<newPts.cols(); ++i){
    assert(newPts(0,i)>=grid->xMin);
    assert(newPts(0,i)<=grid->xMax);
    assert(newPts(1,i)>=grid->yMin);
    assert(newPts(1,i)<=grid->yMax);
  }

  discrPts = newPts;
}

std::tuple<double,Eigen::VectorXd, Eigen::SparseMatrix<double>> SemidiscreteOT::Objective(Eigen::VectorXd const& prices) const
{
    // Notes:
    //   - The cost c(x,y) is the squared distance between x and y
    //   - See (17) of https://arxiv.org/pdf/1710.02634.pdf

    const int numCells = discrPts.cols();
    assert(numCells==prices.size());

    // Construct the Laguerre diagram
    LaguerreDiagram lagDiag(grid->xMin, grid->xMax, grid->yMin, grid->yMax, discrPts, prices);

    double obj;
    Eigen::VectorXd grad;
    Eigen::SparseMatrix<double> hess;

    std::tie(obj,grad) = ComputeGradient(prices, lagDiag);
    hess = ComputeHessian(lagDiag);

    return std::make_tuple(obj,grad,hess);
}

Eigen::Matrix2Xd SemidiscreteOT::PointGradient() const
{
  return PointGradient(*lagDiag);
}

Eigen::Matrix2Xd SemidiscreteOT::PointGradient(LaguerreDiagram const& lagDiag) const
{
  Eigen::Matrix2Xd seeds = lagDiag.SeedPts();
  Eigen::VectorXd areas = lagDiag.Areas();
  Eigen::Matrix2Xd centroids = lagDiag.Centroids();

  return (seeds-centroids)*areas.asDiagonal();
}


std::pair<double,Eigen::VectorXd> SemidiscreteOT::ComputeGradient(Eigen::VectorXd const& prices,
                                                                  LaguerreDiagram const& lagDiag) const
{
 const int numCells = prices.size();

 // Holds the part of the objective for each cell in the Laguerre diagram
 Eigen::VectorXd objParts(numCells);
 Eigen::VectorXd gradient(numCells);

 Eigen::VectorXd probs(numCells);

 //Eigen::MatrixXd cellAreas = Eigen::MatrixXd::Zero(grid->Nx, grid->Ny);
#if defined(_OPENMP)
 #pragma omp parallel for
#endif
 for(int cellInd=0; cellInd<numCells; ++cellInd){

   auto area_integrand = std::make_shared<ConstantIntegrand>();

   auto trans_integrand = std::make_shared<TransportIntegrand>(discrPts.col(cellInd));

   double weightedArea  = lagDiag.IntegrateOverCell(cellInd, area_integrand, dist);

   objParts(cellInd) = prices(cellInd)*discrProbs(cellInd)
                     + lagDiag.IntegrateOverCell(cellInd, trans_integrand, dist)
                     - prices(cellInd)*weightedArea;

   gradient(cellInd) = discrProbs(cellInd) - weightedArea;
   probs(cellInd) = weightedArea;
 }

 double totalProb = probs.sum();
 if(std::abs(totalProb-1.0)>5e-7){

   std::cout << "Warning:  Total probability has an error of " << totalProb-1.0 << std::endl;
   // std::cout << "Prices = " << prices.transpose() << std::endl;
   // for(unsigned int cellInd=0; cellInd<lagDiag.NumCells(); ++cellInd){
   //   std::cout << "Cell " << cellInd << " has points " << std::endl;
   //   Eigen::MatrixXd pts = lagDiag.GetCellVertices(cellInd);
   //   for(unsigned int ptInd=0; ptInd<pts.cols(); ++ptInd){
   //      std::cout << "[" << pts(0,ptInd) << "," << pts(1,ptInd) << "], ";
   //   }
   //   std::cout << std::endl;
   // }
   throw std::runtime_error("Error in total probability.");

   //assert(std::abs(weightedArea-1.0)<1e-10);
 }


 return std::make_pair(objParts.sum(), gradient);
}

double SemidiscreteOT::LineIntegral(LaguerreDiagram::Point_2 const& srcPt,
                                    LaguerreDiagram::Point_2 const& tgtPt) const
{
  const double compTol = 1e-11;

  double xs = CGAL::to_double(srcPt.x());
  double ys = CGAL::to_double(srcPt.y());
  double xt = CGAL::to_double(tgtPt.x());
  double yt = CGAL::to_double(tgtPt.y());

  // if the edge is vertical...
  if(std::abs(xt-xs)<compTol){

    // Assume we are working from bottom to top
    if(yt<ys){
      std::swap(xt,xs);
      std::swap(yt,ys);
    }

    double val = 0.0;
    unsigned int xInd = grid->LeftNode(xs);
    unsigned int yInd = grid->BottomNode(ys);
    const double maxY = CGAL::to_double(yt);

    double currY = CGAL::to_double( ys );
    double nextY = grid->TopNode(ys);

    // If the source point is on a boundary...
    if(std::abs(nextY-currY)<compTol)
      nextY += grid->dy;

    while(nextY<maxY-compTol){
      val += (nextY-currY)*dist->Density(xInd,yInd);

      yInd++;
      currY = nextY;
      nextY = currY+grid->dy;
    }

    val += (maxY-currY)*dist->Density(xInd,yInd);

    return val;

  // If the edge is horizontal
  }else if(std::abs(yt-ys)<compTol){

    // Assume we are working from left to right and swap direction if we're not
    if(xt<xs){
      std::swap(xt,xs);
      std::swap(yt,ys);
    }

    double val = 0.0;
    unsigned int xInd = grid->LeftNode(xs);
    unsigned int yInd = grid->BottomNode(ys);

    const double maxX = CGAL::to_double(xt);

    double currX = CGAL::to_double( xs );
    double nextX = grid->RightNode(xs);

    // If the source is on a boundary...
    if(std::abs(nextX-currX)<compTol)
      nextX += grid->dx;

    while(nextX<maxX-compTol){
      val += (nextX-currX)*dist->Density(xInd,yInd);

      xInd++;
      currX = nextX;
      nextX = currX+grid->dx;
    }

    val += (maxX-currX)*dist->Density(xInd,yInd);

    return val;

  // Otherwise there is a nonzero finite slope
  }else{

    // Assume we are working from left to right and swap direction if we're not
    if(xt<xs){
      std::swap(xt,xs);
      std::swap(yt,ys);
    }

    double val = 0.0;

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


    while(nextt<1.0-compTol){
      val += (nextt-currt)*segLenth*dist->Density(xInd,yInd);

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

    if((xInd<grid->NumCells(0))&&(yInd<grid->NumCells(1)))
      val += (nextt-currt)*segLenth*dist->Density(xInd,yInd);

    return val;
  }
}

Eigen::SparseMatrix<double> SemidiscreteOT::ComputeHessian(LaguerreDiagram const& lagDiag) const
{
  const unsigned int numCells = discrPts.cols();
  typedef Eigen::Triplet<double> T;

  /* The diagonal entries of the Hessian are the negative sum of the off diagonals
     See equation 27 of https://arxiv.org/pdf/1710.02634.pdf
     This vector is used to keep track of this sum for each cell.
  */
  Eigen::VectorXd diagVals = Eigen::VectorXd::Zero(numCells);

  // Hold the i,j,val triplets defining the sparse Hessian
  std::vector<T> hessVals;

  double intVal;
  unsigned int cellInd2;
  LaguerreDiagram::Point_2 srcPt, tgtPt;

  for(unsigned int cellInd1=0; cellInd1<numCells; ++cellInd1){

    for(auto edgeTuple : lagDiag.InternalEdges(cellInd1)){
      std::tie(cellInd2, srcPt, tgtPt) = edgeTuple;

      // Compute the integral of the target density along the edge
      intVal = 0.25*LineIntegral(srcPt,tgtPt)/(discrPts.col(cellInd1)-discrPts.col(cellInd2)).norm();

      diagVals(cellInd1) -= intVal;

      hessVals.push_back(T(cellInd1,cellInd2,intVal));
    }
  }

  for(int i=0; i<numCells; ++i)
    hessVals.push_back(T(i,i,diagVals(i)));

  Eigen::SparseMatrix<double> hess(numCells,numCells);
  hess.setFromTriplets(hessVals.begin(), hessVals.end());

  return hess;
}

// double SemidiscreteOT::SquareIntegral(double xmin, double xmax,
//                                       double ymin, double ymax,
//                                       double px,   double py)
// {
//  double rectInt = (0.5/3.0)*(ymax-ymin)*(std::pow(xmax-px,3.0)-std::pow(xmin-px,3.0))
//                 + (0.5/3.0)*(xmax-xmin)*(std::pow(ymax-py,3.0)-std::pow(ymin-py,3.0));
//
//  return rectInt;
// }

// double SemidiscreteOT::TriangleIntegral(double x1, double y1,
//                                        double x2, double y2,
//                                        double x3, double y3,
//                                        double px, double py)
// {
//  double triInt = (1.0/2.0)*std::pow(px, 2) - 1.0/3.0*px*x1 - 1.0/3.0*px*x2 - 1.0/3.0*px*x3 + (1.0/2.0)*std::pow(py, 2) - 1.0/3.0*py*y1 - 1.0/3.0*py*y2 - 1.0/3.0*py*y3 + (1.0/12.0)*std::pow(x1, 2) + (1.0/12.0)*x1*x2 + (1.0/12.0)*x1*x3 + (1.0/12.0)*std::pow(x2, 2) + (1.0/12.0)*x2*x3 + (1.0/12.0)*std::pow(x3, 2) + (1.0/12.0)*std::pow(y1, 2) + (1.0/12.0)*y1*y2 + (1.0/12.0)*y1*y3 + (1.0/12.0)*std::pow(y2, 2) + (1.0/12.0)*y2*y3 + (1.0/12.0)*std::pow(y3, 2);
//  triInt *= 0.5*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
//
//  return triInt;
// }


std::pair<Eigen::VectorXd, double> SemidiscreteOT::Solve(Eigen::VectorXd                  const& prices0,
                                                         OptionList                              options)
{
  assert(prices0.size()==discrPts.cols());
  const unsigned int dim = prices0.size();

  unsigned int printLevel = GetOpt("Print Level", options, 3);

  // Trust region approach with a double dogleg step
  double trustRadius = GetOpt("Trust Radius", options, 1.0);
  const unsigned int maxEvals = GetOpt("Max Steps", options, 100.0);

  const double xtol_abs = GetOpt("XTol Abs", options, 1e-13*std::sqrt(double(dim)));
  const double gtol_abs = GetOpt("GTol Abs", options, 2e-4*std::sqrt(double(dim)));
  const double ftol_abs = GetOpt("FTol Abs", options, 1e-11);

  const double acceptRatio = GetOpt("Accept Ratio", options, 1e-4);//0.1;
  const double shrinkRatio = GetOpt("Shrink Ratio", options, 1e-4);//0.1;
  const double growRatio = GetOpt("Grow Ratio", options, 0.75);
  const double growRate = GetOpt("Grow Rate", options, 2.0);
  const double shrinkRate = GetOpt("Shrink Rate", options, 0.25);
  const double maxRadius = GetOpt("Max Radius", options, 10);



  double fval, newF, gradNorm, newGradNorm;
  Eigen::VectorXd grad, newGrad;
  Eigen::SparseMatrix<double> hess;

  Eigen::VectorXd x = prices0;
  Eigen::VectorXd newX(x);
  Eigen::VectorXd step = Eigen::VectorXd::Zero(dim);

  std::shared_ptr<LaguerreDiagram> newLagDiag;

  // Compute an initial gradient and Hessian
  lagDiag  = std::make_shared<LaguerreDiagram>(grid->xMin, grid->xMax, grid->yMin, grid->yMax, discrPts, x);
  assert(lagDiag!=nullptr);

  std::tie(fval, grad) = ComputeGradient(x, *lagDiag);

  fval *= -1.0;
  grad *= -1.0;
  gradNorm = grad.norm();

  if(printLevel>0){
    std::cout << "Using NewtonTrust optimizer..." << std::endl;
    std::cout << "  Iteration, TrustRadius,       ||g||,   ||g||/sqrt(dim)" << std::endl;
  }

  for(int it=0; it<maxEvals; ++it) {

    auto start = std::chrono::steady_clock::now();
    hess = ComputeHessian(*lagDiag);
    hess *= -1.0;

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> hess_time = end-start;

    if(printLevel>0){
      std::printf("  %9d, %11.2e,  %5.3e, %5.3e\n", it, trustRadius, gradNorm, gradNorm/std::sqrt(double(dim)));
    }

    if(gradNorm < gtol_abs){
      if(printLevel>0){
        std::printf("Terminating because gradient norm (%4.2e) is smaller than gtol_abs (%4.2e).\n", gradNorm, gtol_abs);
      }
      return std::make_pair(x,fval);
    }

    start = std::chrono::steady_clock::now();
    step.tail(dim-1) = SolveSubProblem(fval, grad.tail(dim-1),  hess.block(1,1,dim-1,dim-1), trustRadius);
    end = std::chrono::steady_clock::now();
    std::chrono::duration<double> sub_time = end-start;

    // Crude method (kind of like a line search) for ensuring the prices are positive
    newX = x+step;
    while(newX.minCoeff()<1e-15){
      std::cout << "            Shrinking trust region because of positivity constraint." << std::endl;
      step *= 0.5;
      newX = x+step;
    }

    // Try constructing the new Laguerre diagram.  If we can't then shrink the trust region size
    try{

      start = std::chrono::steady_clock::now();
      newLagDiag  = std::make_shared<LaguerreDiagram>(grid->xMin, grid->xMax, grid->yMin, grid->yMax, discrPts, newX);
      end = std::chrono::steady_clock::now();
      std::chrono::duration<double> lag_time = end-start;

      start = std::chrono::steady_clock::now();
      std::tie(newF, newGrad) = ComputeGradient(newX, *newLagDiag);
      end = std::chrono::steady_clock::now();
      std::chrono::duration<double> grad_time = end-start;

      // std::cout << "Times:" << std::endl;
      // std::cout << "    Hessian: "  << hess_time.count()  << std::endl;
      // std::cout << "    Subproblem: " << sub_time.count()  << std::endl;
      // std::cout << "    Laguerre Construction: " << lag_time.count()  << std::endl;
      // std::cout << "    Gradient: " << grad_time.count()  << std::endl;


      newF *= -1.0;
      newGrad *= -1.0;

      newGradNorm = newGrad.norm();

      // Use the quadratic submodel to predict the change in the gradient norm
      double modDelta = gradNorm - (grad + hess.selfadjointView<Eigen::Lower>()*step).norm();
      double trueDelta = gradNorm - newGradNorm;

      double rho = trueDelta/modDelta;
      // std::cout << "          step.dot(grad) = " << step.dot(grad) << std::endl;
      // std::cout << "          delta f = " << trueDelta << std::endl;
      // std::cout << "          modDelta = " << modDelta << std::endl;
      // std::cout << "          New prices = " << newX.transpose() << std::endl;
      // std::cout << "          rho = " << rho << std::endl;

      double stepNorm = step.norm();
      if(stepNorm < xtol_abs){
        if(printLevel>0){
          std::printf("Terminating because stepsize (%4.2e) is smaller than xtol_abs (%4.2e).\n", stepNorm, xtol_abs);
        }
        return std::make_pair(newX,newF);
      }

      // Update the position.  If the model is really bad, we'll just stay put
      if(rho>acceptRatio){

        if(std::abs(fval-newF)<ftol_abs){
          if(printLevel>0){
            std::printf("Terminating because change in objective (%4.2e) is smaller than ftol_abs (%4.2e).\n", fval-newF, ftol_abs);
          }
          return std::make_pair(newX,newF);
        }

        x = newX;
        fval = newF;
        lagDiag = newLagDiag;
        grad = newGrad;
        gradNorm = newGradNorm;
      }

      // Update the trust region size
      if(rho<shrinkRatio){
        trustRadius = shrinkRate*trustRadius; // shrink trust region

        if(printLevel>1)
          std::cout << "            Shrinking trust region because of submodel accuracy." << std::endl;

      }else if((rho>growRatio)&&(std::abs(step.norm()-trustRadius)<1e-10)) {
        trustRadius = std::min(growRate*trustRadius, maxRadius);

        if(printLevel>1)
          std::cout << "            Growing trust region." << std::endl;

      }

    }catch(LaguerreDiagram::ConstructionException& e){
      trustRadius = shrinkRate*trustRadius;
      if(printLevel>1)
        std::cout << "            Shrinking trust region because of degenerate Laguerre diagram." << std::endl;
    }
  }

  if(printLevel>0){
    std::printf("Terminating because maximum number of iterations (%d) was reached.", maxEvals);
  }
  return std::make_pair(x,fval);
}


Eigen::VectorXd SemidiscreteOT::SolveSubProblem(double obj,
                                                Eigen::Ref<const Eigen::VectorXd> const& grad,
                                                Eigen::Ref<const Eigen::SparseMatrix<double>> const& hess,
                                                double trustRadius) const
{
  const double trustTol = 1e-12;
  const unsigned int dim = grad.size();

  // Current estimate of the subproblem minimum
  Eigen::VectorXd z = Eigen::VectorXd::Zero(dim);

  // Related to the step direction
  Eigen::VectorXd r = grad;
  Eigen::VectorXd d = -r;

  // If the gradient is small enough where we're starting, then we're done
  if(r.norm()<trustTol){
    return z;
  }

  Eigen::VectorXd Bd; // the Hessian (B) applied to a vector d

  double alpha, beta, gradd, dBd, rr;

  for(int i=0; i<dim; ++i){
    Bd = hess.selfadjointView<Eigen::Lower>()*d;
    gradd = grad.dot(d);
    dBd = d.dot(Bd);
    rr = r.squaredNorm();

    // If the Hessian isn't positive definite in this direction, we can go all
    // the way to the trust region boundary
    if(dBd<=0){
      // do something

      double dz = d.dot(z);
      double dd = d.squaredNorm();
      double zz = z.squaredNorm();
      double r2 = trustRadius*trustRadius;

      double tau1 = (-dz + sqrt(dz*dz - dd*(zz-r2)))/dd;
      double tau2 = (-dz - sqrt(dz*dz - dd*(zz-r2)))/dd;

      double zBd = z.dot(Bd);
      double mval1 = tau1*gradd + tau1*zBd + tau1*tau1*dBd;
      double mval2 = tau2*gradd + tau2*zBd + tau2*tau2*dBd;

      return (mval1<mval2) ? (z+tau1*d) : (z+tau2*d);
    }

    alpha = rr / dBd;
    Eigen::VectorXd newZ = z + alpha * d;

    if(newZ.norm()>trustRadius){

      double dz = d.dot(z);
      double dd = d.squaredNorm();
      double zz = z.squaredNorm();
      double r2 = trustRadius*trustRadius;

      double tau = (-dz + sqrt(dz*dz - dd*(zz-r2)))/dd;
      return z + tau*d;
    }

    z = newZ;

    r += alpha*Bd;

    if(r.norm()<trustTol){
      return z;
    }

    beta = r.squaredNorm() / rr;
    d = (-r + beta*d).eval();
  }

  return z;
}


std::shared_ptr<LaguerreDiagram> SemidiscreteOT::BuildCentroidal(std::shared_ptr<Distribution2d> const& dist,
                                                                 Eigen::Matrix2Xd                const& initialPoints,
                                                                 Eigen::VectorXd                 const& pointProbs,
                                                                 OptionList                             opts)
{

  unsigned int maxIts =  GetOpt("Lloyd Steps", opts, 100);
  double tol = GetOpt("Lloyd Tol", opts, 1e-3);

  const unsigned int numPts = pointProbs.size();
  assert(numPts==initialPoints.cols());

  double resid = tol + 1.0;

  Eigen::MatrixXd newPts;

  Eigen::MatrixXd pts = initialPoints;
  std::shared_ptr<SemidiscreteOT> ot = std::make_shared<SemidiscreteOT>(dist, initialPoints, pointProbs);

  std::cout << "Computing constrained centroidal diagram." << std::endl;
  for(unsigned int i=0; i<maxIts; ++i){
    assert(ot);
    ot->Solve(Eigen::VectorXd::Ones(numPts), opts);

    newPts = ot->Diagram()->Centroids(dist);
    resid = (newPts - pts).cwiseAbs().maxCoeff();;
    if(resid<tol){
      std::cout << "Converged with a final step of " << resid << std::endl;
      return ot->Diagram();
    }

    pts = newPts;
    ot->SetPoints(pts);

    std::cout << "  After " << i << " iterations, change in position = " << resid << std::endl;
  }

  std::cout << "WARNING: Did not converge to constrained centroidal diagram." << std::endl;
  return ot->Diagram();
}

std::shared_ptr<LaguerreDiagram> SemidiscreteOT::BuildCentroidal(std::shared_ptr<Distribution2d> const& dist,
                                                                 Eigen::VectorXd                 const& probs,
                                                                 OptionList                             opts)
{
  BoundingBox bbox(dist->Grid()->xMin, dist->Grid()->xMax, dist->Grid()->yMin, dist->Grid()->yMax);
  Eigen::Matrix2Xd initialPts = LaguerreDiagram::LatinHypercubeSample(bbox, probs.size());
  return BuildCentroidal(dist, initialPts, probs, opts);
}

std::shared_ptr<LaguerreDiagram> SemidiscreteOT::BuildCentroidal(std::shared_ptr<Distribution2d> const& dist,
                                                                 unsigned int                           numPts,
                                                                 OptionList                             opts)
{
  Eigen::VectorXd probs = (1.0/numPts)*Eigen::VectorXd::Ones(numPts);
  return BuildCentroidal(dist, probs, opts);
}
