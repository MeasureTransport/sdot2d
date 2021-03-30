#include "SDOT/SemiDiscreteOT.h"


#include "SDOT/Assert.h"
#include "SDOT/Distances/Distances.h"

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

#include <CGAL/Kernel/global_functions.h>

#include <algorithm>
#include <chrono>

using namespace sdot;

template<typename ConjugateFunctionType>
SemidiscreteOT<ConjugateFunctionType>::SemidiscreteOT(std::shared_ptr<Distribution2d> const& distIn,
                               Eigen::Matrix2Xd                const& discrPtsIn,
                               Eigen::VectorXd                 const& discrProbsIn,
                               double                                 unbalancedPenaltyIn) : dist(distIn),
                                                                                      grid(distIn->Grid()),
                                                                                      discrPts(discrPtsIn),
                                                                                      discrProbs(discrProbsIn),
                                                                                      unbalancedPenalty(unbalancedPenaltyIn)
{
  if(discrPtsIn.cols()!=discrProbsIn.size())
  SDOT_ASSERT(discrPtsIn.cols()==discrProbsIn.size());

  SDOT_ASSERT(unbalancedPenalty>0);

  CheckNormalization();

  // Check to make sure all the points are inside the grid domain
  for(unsigned int i=0; i<discrPts.cols(); ++i){
    SDOT_ASSERT(discrPts(0,i)>=grid->xMin);
    SDOT_ASSERT(discrPts(0,i)<=grid->xMax);
    SDOT_ASSERT(discrPts(1,i)>=grid->yMin);
    SDOT_ASSERT(discrPts(1,i)<=grid->yMax);
  }
}


template<>
void SemidiscreteOT<sdot::distances::Wasserstein2>::CheckNormalization()
{
  double discrSum = discrProbs.sum();
  double contSum= dist->TotalMass();
  SDOT_ASSERT(std::abs(discrSum-contSum)<1e-10);
}

template<typename ConjugateFunctionType>
void SemidiscreteOT<ConjugateFunctionType>::SetPoints(Eigen::Matrix2Xd const& newPts){
  SDOT_ASSERT(newPts.cols()==discrProbs.size());

  // Check to make sure all the points are inside the grid domain
  for(unsigned int i=0; i<newPts.cols(); ++i){
    SDOT_ASSERT(newPts(0,i)>=grid->xMin);
    SDOT_ASSERT(newPts(0,i)<=grid->xMax);
    SDOT_ASSERT(newPts(1,i)>=grid->yMin);
    SDOT_ASSERT(newPts(1,i)<=grid->yMax);
  }

  discrPts = newPts;
}

template<typename ConjugateFunctionType>
std::tuple<double,Eigen::VectorXd, Eigen::SparseMatrix<double>> SemidiscreteOT<ConjugateFunctionType>::Objective(Eigen::VectorXd const& prices) const
{
    // Notes:
    //   - The cost c(x,y) is the squared distance between x and y
    //   - See (17) of https://arxiv.org/pdf/1710.02634.pdf

    const int numCells = discrPts.cols();
    SDOT_ASSERT(numCells==prices.size());

    // Construct the Laguerre diagram
    LaguerreDiagram lagDiag(grid->xMin, grid->xMax, grid->yMin, grid->yMax, discrPts, prices);

    double obj;
    Eigen::VectorXd grad;
    Eigen::SparseMatrix<double> hess;

    std::tie(obj,grad) = ComputeGradient(prices, lagDiag);
    hess = ComputeHessian(prices, lagDiag);

    return std::make_tuple(obj,grad,hess);
}

template<typename ConjugateFunctionType>
Eigen::Matrix2Xd SemidiscreteOT<ConjugateFunctionType>::PointGradient() const
{
  return PointGradient(optPrices, *lagDiag);
}

template<typename ConjugateFunctionType>
Eigen::Matrix2Xd SemidiscreteOT<ConjugateFunctionType>::PointGradient(Eigen::VectorXd const& prices,
                                                                      LaguerreDiagram const& lagDiag) const
{
  int numCells =  lagDiag.NumCells();
  Eigen::Matrix2Xd grad(2,numCells);

  // Loop over all ofthe cells
  for(int i=0; i<numCells; ++i){

    auto triFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2, Eigen::Vector2d const& pt3)
      {
        return ConjugateFunctionType::TriangularIntegralPointGrad(prices(i), discrPts.col(i), pt1,pt2,pt3,unbalancedPenalty);
      };
    auto rectFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
      {
        return ConjugateFunctionType::RectangularIntegralPointGrad(prices(i), discrPts.col(i), pt1,pt2,unbalancedPenalty);
      };

    grad.col(i) = lagDiag.IntegrateOverCell(i, triFunc, rectFunc, dist);
  }

  return grad;
}

template<typename ConjugateFunctionType>
Eigen::Matrix2Xd SemidiscreteOT<ConjugateFunctionType>::LloydPointHessian() const
{
  return LloydPointHessian(optPrices, *lagDiag);
}

template<typename ConjugateFunctionType>
Eigen::Matrix2Xd SemidiscreteOT<ConjugateFunctionType>::LloydPointHessian(Eigen::VectorXd const& prices,
                                                                         LaguerreDiagram const& lagDiag) const
{
  const unsigned int numCells = discrPts.cols();
  Eigen::Matrix2Xd hessVals(2,numCells);
  Eigen::Matrix2d intVal(2,2);

  // For unbalanced transport, the Diagonal of the hessian has an additional term
  for(unsigned int cellInd=0; cellInd<numCells; ++cellInd) {

    auto triFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2, Eigen::Vector2d const& pt3)
      {
        return ConjugateFunctionType::TriangularIntegralPointHessDiag(prices(cellInd), discrPts.col(cellInd), pt1,pt2,pt3, unbalancedPenalty);
      };
    auto rectFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
      {
        return ConjugateFunctionType::RectangularIntegralPointHessDiag(prices(cellInd), discrPts.col(cellInd), pt1,pt2, unbalancedPenalty);
      };

    intVal = -1.0*lagDiag.IntegrateOverCell(cellInd, triFunc, rectFunc, dist);
    hessVals(0,cellInd)   =  intVal(0,0);
    hessVals(1,cellInd)   =  intVal(1,1);
  }

  return hessVals;
}


template<typename ConjugateFunctionType>
Eigen::SparseMatrix<double> SemidiscreteOT<ConjugateFunctionType>::PointHessian(Eigen::VectorXd const& prices,
                                                                                LaguerreDiagram const& lagDiag) const
{
  const unsigned int numCells = discrPts.cols();
  typedef Eigen::Triplet<double> T;
  std::vector<T> hessVals;

  Eigen::VectorXd diagVals = Eigen::VectorXd::Zero(2*numCells);
  Eigen::Matrix2d intVal;

  // For unbalanced transport, the Diagonal of the hessian has an additional term
  for(unsigned int cellInd=0; cellInd<numCells; ++cellInd) {

    auto triFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2, Eigen::Vector2d const& pt3)
      {
        return ConjugateFunctionType::TriangularIntegralPointHessDiag(prices(cellInd), discrPts.col(cellInd), pt1,pt2,pt3, unbalancedPenalty);
      };
    auto rectFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
      {
        return ConjugateFunctionType::RectangularIntegralPointHessDiag(prices(cellInd), discrPts.col(cellInd), pt1,pt2, unbalancedPenalty);
      };

    intVal = -1.0*lagDiag.IntegrateOverCell(cellInd, triFunc, rectFunc, dist);
    hessVals.push_back(T(2*cellInd,2*cellInd,intVal(0,0)));
    hessVals.push_back(T(2*cellInd,2*cellInd+1,intVal(0,1)));
    hessVals.push_back(T(2*cellInd+1,2*cellInd,intVal(1,0)));
    hessVals.push_back(T(2*cellInd+1,2*cellInd+1,intVal(1,1)));
  }


  ///////
  // Off-diagonal parts

  unsigned int cellInd2;
  LaguerreDiagram::Point_2 srcPt, tgtPt;

  for(unsigned int cellInd1=0; cellInd1<numCells; ++cellInd1){

    for(auto edgeTuple : lagDiag.InternalEdges(cellInd1)){
      std::tie(cellInd2, srcPt, tgtPt) = edgeTuple;

      // Compute the integral of the target density along the edge
      auto func = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
      {
        return ConjugateFunctionType::LineIntegralPointHess(prices(cellInd1), discrPts.col(cellInd1), discrPts.col(cellInd2), pt1, pt2, unbalancedPenalty);
      };
      intVal = LineIntegral(func, srcPt, tgtPt)/(discrPts.col(cellInd1)-discrPts.col(cellInd2)).norm();

      hessVals.push_back(T(2*cellInd1,2*cellInd1,-intVal(0,0)));
      hessVals.push_back(T(2*cellInd1+1,2*cellInd1,-intVal(0,1)));
      hessVals.push_back(T(2*cellInd1,2*cellInd1+1,-intVal(1,0)));
      hessVals.push_back(T(2*cellInd1+1,2*cellInd1+1,-intVal(1,1)));

      hessVals.push_back(T(2*cellInd1,2*cellInd2,intVal(0,0)));
      hessVals.push_back(T(2*cellInd1+1,2*cellInd2,intVal(0,1)));
      hessVals.push_back(T(2*cellInd1,2*cellInd2+1,intVal(1,0)));
      hessVals.push_back(T(2*cellInd1+1,2*cellInd2+1,intVal(1,1)));

    }
  }

  Eigen::SparseMatrix<double> hess(2*numCells,2*numCells);
  hess.setFromTriplets(hessVals.begin(), hessVals.end());

  return hess;

}

template<typename ConjugateFunctionType>
Eigen::SparseMatrix<double> SemidiscreteOT<ConjugateFunctionType>::PointHessian() const
{
  return PointHessian(optPrices, *lagDiag);
}

template<typename ConjugateFunctionType>
std::pair<double,Eigen::VectorXd> SemidiscreteOT<ConjugateFunctionType>::ComputeGradient(Eigen::VectorXd const& prices,
                                                                                         LaguerreDiagram const& lagDiag) const
{
  const int numCells = prices.size();

  // Holds the part of the objective for each cell in the Laguerre diagram
  Eigen::VectorXd objParts = Eigen::VectorXd::Zero(numCells);
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(numCells);

  Eigen::VectorXd probs(numCells);

  //Eigen::MatrixXd cellAreas = Eigen::MatrixXd::Zero(grid->Nx, grid->Ny);
#if defined(_OPENMP)
  #pragma omp parallel for
#endif
  for(int cellInd=0; cellInd<numCells; ++cellInd){

    // auto area_integrand = std::make_shared<ConstantIntegrand>();
    //
    // auto trans_integrand = std::make_shared<TransportIntegrand>(discrPts.col(cellInd));

    // auto triArea = [](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2, Eigen::Vector2d const& pt3)
    //   {
    //     return 0.5*std::abs((pt2[0]*pt1[1]-pt1[0]*pt2[1])+(pt3[0]*pt2[1]-pt2[0]*pt3[1])+(pt1[0]*pt3[1]-pt3[0]*pt1[1]));
    //   };
    //
    // auto rectArea = [](Eigen::Vector2d const& bottomLeft, Eigen::Vector2d const& topRight)
    //   {
    //     return std::abs((topRight[0]-bottomLeft[0])*(topRight[1]-bottomLeft[1]));
    //   };
    //
    //
    // double weightedArea  = lagDiag.IntegrateOverCell(cellInd, triArea, rectArea, dist);

    objParts(cellInd) = -ConjugateFunctionType::Evaluate(-prices(cellInd), unbalancedPenalty)*discrProbs(cellInd);// - prices(cellInd)*weightedArea;

    gradient(cellInd) = ConjugateFunctionType::Derivative(-prices(cellInd), unbalancedPenalty)*discrProbs(cellInd);

    // if(weightedArea>0){

      auto triFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2, Eigen::Vector2d const& pt3)
        {
          return ConjugateFunctionType::TriangularIntegral(prices(cellInd), discrPts.col(cellInd), pt1,pt2,pt3,unbalancedPenalty);
        };
      auto rectFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
        {
          return ConjugateFunctionType::RectangularIntegral(prices(cellInd), discrPts.col(cellInd), pt1,pt2,unbalancedPenalty);
        };

      objParts(cellInd) += -lagDiag.IntegrateOverCell(cellInd, triFunc, rectFunc, dist);

      auto triFuncDeriv = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2, Eigen::Vector2d const& pt3)
        {
          return ConjugateFunctionType::TriangularIntegralDeriv(prices(cellInd), discrPts.col(cellInd), pt1,pt2,pt3,unbalancedPenalty);
        };
      auto rectFuncDeriv = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
        {
          return ConjugateFunctionType::RectangularIntegralDeriv(prices(cellInd), discrPts.col(cellInd), pt1,pt2,unbalancedPenalty);
        };

      gradient(cellInd) += -lagDiag.IntegrateOverCell(cellInd, triFuncDeriv, rectFuncDeriv, dist); //-weightedArea;
    }

  //   probs(cellInd) = weightedArea;
  // }
  //
  // double totalProb = probs.sum();
  // if(std::abs(totalProb-1.0)>1e-10){
  //
  //   std::cout << "Warning:  Total probability has an error of " << totalProb-1.0 << std::endl;
  //   // std::cout << "Prices = " << prices.transpose() << std::endl;
  //   // for(unsigned int cellInd=0; cellInd<lagDiag.NumCells(); ++cellInd){
  //   //   std::cout << "Cell " << cellInd << " has points " << std::endl;
  //   //   Eigen::MatrixXd pts = lagDiag.GetCellVertices(cellInd);
  //   //   for(unsigned int ptInd=0; ptInd<pts.cols(); ++ptInd){
  //   //      std::cout << "[" << pts(0,ptInd) << "," << pts(1,ptInd) << "], ";
  //   //   }
  //   //   std::cout << std::endl;
  //   // }
  //   throw std::runtime_error("Error in total probability.");
  //
  //
  //
  //  //SDOT_ASSERT(std::abs(weightedArea-1.0)<1e-10);
  // }

  return std::make_pair(objParts.sum(), gradient);
}


template<typename ConjugateFunctionType>
Eigen::SparseMatrix<double> SemidiscreteOT<ConjugateFunctionType>::ComputeHessian(Eigen::VectorXd const& prices,
                                                                                  LaguerreDiagram const& lagDiag) const
{
  const unsigned int numCells = discrPts.cols();
  typedef Eigen::Triplet<double> T;

  /* The diagonal entries of the Hessian are the negative sum of the off diagonals
     See equation 27 of https://arxiv.org/pdf/1710.02634.pdf
     This vector is used to keep track of this sum for each cell.
  */
  Eigen::VectorXd diagVals = Eigen::VectorXd::Zero(numCells);

  // For unbalanced transport, the Diagonal of the hessian has an additional term
  for(unsigned int cellInd=0; cellInd<numCells; ++cellInd) {

    diagVals(cellInd) -= ConjugateFunctionType::Derivative2(-prices(cellInd), unbalancedPenalty)*discrProbs(cellInd);

    auto triFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2, Eigen::Vector2d const& pt3)
      {
        return ConjugateFunctionType::TriangularIntegralDeriv2(prices(cellInd), discrPts.col(cellInd), pt1,pt2,pt3, unbalancedPenalty);
      };
    auto rectFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
      {
        return ConjugateFunctionType::RectangularIntegralDeriv2(prices(cellInd), discrPts.col(cellInd), pt1,pt2, unbalancedPenalty);
      };

    diagVals(cellInd) -= lagDiag.IntegrateOverCell(cellInd, triFunc, rectFunc, dist);
  }

  // Hold the i,j,val triplets defining the sparse Hessian
  std::vector<T> hessVals;

  double intVal;
  unsigned int cellInd2;
  LaguerreDiagram::Point_2 srcPt, tgtPt;

  for(unsigned int cellInd1=0; cellInd1<numCells; ++cellInd1){

    for(auto edgeTuple : lagDiag.InternalEdges(cellInd1)){
      std::tie(cellInd2, srcPt, tgtPt) = edgeTuple;

      // Compute the integral of the target density along the edge
      auto func = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
      {
        return ConjugateFunctionType::LineIntegralDeriv(prices(cellInd1), discrPts.col(cellInd1), pt1, pt2, unbalancedPenalty);
      };
      intVal = 0.5*LineIntegral(func, srcPt, tgtPt)/(discrPts.col(cellInd1)-discrPts.col(cellInd2)).norm();

      diagVals(cellInd1) -= intVal;
      hessVals.push_back(T(cellInd1,cellInd2,intVal));
    }
  }

  for(int i=0; i<numCells; ++i){
    hessVals.push_back(T(i,i, diagVals(i)));
  }

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


// template<typename ConjugateFunctionType>
// Eigen::VectorXd SemidiscreteOT<ConjugateFunctionType>::GetValidPrices(Eigen::VectorXd const& x0)
// {
//
// }

template<typename ConjugateFunctionType>
std::pair<Eigen::VectorXd, double> SemidiscreteOT<ConjugateFunctionType>::Solve(Eigen::VectorXd                  const& prices0,
                                                         OptionList                              options)
{
  SDOT_ASSERT(prices0.size()==discrPts.cols());
  const unsigned int dim = prices0.size();

  unsigned int printLevel = GetOpt("Print Level", options, 3);

  // Trust region approach with a double dogleg step
  double trustRadius = GetOpt("Trust Radius", options, 1.0);
  const unsigned int maxEvals = GetOpt("Max Steps", options, 100.0);

  const double xtol_abs = GetOpt("XTol Abs", options, 1e-13*std::sqrt(double(dim)));
  const double gtol_abs = GetOpt("GTol Abs", options, 2e-4*std::sqrt(double(dim)));
  const double ftol_abs = GetOpt("FTol Abs", options, 1e-11);

  const double acceptRatio = GetOpt("Accept Ratio", options, 0.1);//0.1;
  const double shrinkRatio = GetOpt("Shrink Ratio", options, 0.25);//0.1;
  const double growRatio = GetOpt("Grow Ratio", options, 0.75);
  const double growRate = GetOpt("Grow Rate", options, 2.0);
  const double shrinkRate = GetOpt("Shrink Rate", options, 0.25);
  const double maxRadius = GetOpt("Max Radius", options, 10);

  SDOT_ASSERT(shrinkRatio>=acceptRatio);

  double fval, newF, gradNorm, newGradNorm;
  Eigen::VectorXd grad, newGrad;
  Eigen::SparseMatrix<double> hess;

  Eigen::VectorXd x = prices0;
  Eigen::VectorXd newX(x);
  Eigen::VectorXd step = Eigen::VectorXd::Zero(dim);

  std::shared_ptr<LaguerreDiagram> newLagDiag;

  // Compute an initial gradient and Hessian
  lagDiag  = std::make_shared<LaguerreDiagram>(grid->xMin, grid->xMax, grid->yMin, grid->yMax, discrPts, x);
  SDOT_ASSERT(lagDiag!=nullptr);

  std::tie(fval, grad) = ComputeGradient(x, *lagDiag);
  hess = ComputeHessian(x, *lagDiag);

  fval *= -1.0;
  grad *= -1.0;
  hess *= -1.0;
  gradNorm = grad.norm();

  if(printLevel>0){
    std::cout << "Using NewtonTrust optimizer..." << std::endl;
    std::cout << "  Iteration, TrustRadius, Dual Objective,        ||g||,   ||g||/sqrt(dim)" << std::endl;
  }

  // count the non-empty Cells
  int numEmpty = 0;
  for(int i=0; i<x.size(); ++i){
    if(lagDiag->GetCellVertices(i).cols()<3)
      ++numEmpty;
  }


  for(int it=0; it<maxEvals; ++it) {

    if(printLevel>0){
      std::printf("  %9d, %11.2e, %14.3e,        %5.3e,   %15.3e\n", it, trustRadius, -1.0*fval, gradNorm, gradNorm/std::sqrt(double(dim)));
    }

    if((gradNorm < gtol_abs)&&(numEmpty==0)){
      if(printLevel>0){
        std::printf("Terminating because gradient norm (%4.2e) is smaller than gtol_abs (%4.2e).\n", gradNorm, gtol_abs);
      }
      optPrices = x;
      return std::make_pair(x,fval);
    }

    step = SolveSubProblem(fval, grad,  hess, trustRadius);

    newX = x+step;

    // Try constructing the new Laguerre diagram.
    newLagDiag  = std::make_shared<LaguerreDiagram>(grid->xMin, grid->xMax, grid->yMin, grid->yMax, discrPts, newX);

    std::tie(newF, newGrad) = ComputeGradient(newX, *newLagDiag);
    newF *= -1.0;
    newGrad *= -1.0;

    //newGradNorm = newGrad.norm();

    // Use the quadratic submodel to predict the change in the objective
    double trueDelta = newF-fval;
    double modDelta = grad.dot(step) + 0.5*step.dot(hess.selfadjointView<Eigen::Lower>()*step);

    double rho = trueDelta/modDelta;
    //std::cout << "Model, Truth = " << modDelta << ", " << trueDelta << std::endl;
    //std::cout << "          step.dot(grad) = " << step.dot(grad) << std::endl;
    // std::cout << "          delta f = " << trueDelta << std::endl;
    // std::cout << "          modDelta = " << modDelta << std::endl;
    // std::cout << "          New prices = " << newX.transpose() << std::endl;
    //std::cout << "          rho = " << rho << std::endl;

    double stepNorm = step.norm();
    if((stepNorm < xtol_abs)&&(numEmpty==0)){
      if(printLevel>0){
        std::printf("Terminating because stepsize (%4.2e) is smaller than xtol_abs (%4.2e).\n", stepNorm, xtol_abs);
      }
      optPrices = newX;
      return std::make_pair(newX,newF);
    }

    // Update the position.  If the model is really bad, we'll just stay put
    if(rho>acceptRatio){

      if(((std::abs(fval-newF)<ftol_abs))&&(numEmpty==0)){
        if(printLevel>0){
          std::printf("Terminating because change in objective (%4.2e) is smaller than ftol_abs (%4.2e).\n", fval-newF, ftol_abs);
        }
        optPrices = newX;
        return std::make_pair(newX,newF);
      }

      x = newX;
      fval = newF;
      lagDiag = newLagDiag;
      grad = newGrad;
      gradNorm = newGradNorm;

      // Recompute the Hessian at the new point
      hess = ComputeHessian(x, *lagDiag);
      hess *= -1.0;

      // Update the number of empty cells in the diagram
      numEmpty = 0;
      for(int i=0; i<x.size(); ++i){
        if(lagDiag->GetCellVertices(i).cols()<3)
          ++numEmpty;
      }
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
  }

  if(printLevel>0){
    std::printf("Terminating because maximum number of iterations (%d) was reached.", maxEvals);
  }

  optPrices = x;
  return std::make_pair(x,fval);
}

template<typename ConjugateFunctionType>
Eigen::VectorXd SemidiscreteOT<ConjugateFunctionType>::SolveSubProblem(double obj,
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

template<typename ConjugateFunctionType>
std::shared_ptr<LaguerreDiagram> SemidiscreteOT<ConjugateFunctionType>::BuildCentroidal(std::shared_ptr<Distribution2d> const& dist,
                                                                 Eigen::Matrix2Xd                const& initialPoints,
                                                                 Eigen::VectorXd                 const& pointProbs,
                                                                 OptionList                             opts)
{

  unsigned int maxIts =  GetOpt("Lloyd Steps", opts, 100);
  double xtol = GetOpt("Lloyd Tol", opts, 1e-8);
  double gtol = xtol;

  auto optIt = opts.find("Print Level");
  if(optIt == opts.end())
    opts["Print Level"] = 0;

  const unsigned int numPts = pointProbs.size();
  SDOT_ASSERT(numPts==initialPoints.cols());

  double resid = xtol + 1.0;

  Eigen::MatrixXd newPts;
  Eigen::VectorXd prices = Eigen::VectorXd::Ones(numPts);
  double dualObj;
  Eigen::MatrixXd pts = initialPoints;
  std::shared_ptr<SemidiscreteOT> ot = std::make_shared<SemidiscreteOT>(dist, initialPoints, pointProbs, GetOpt("Penalty", opts, 1.0));

  std::cout << "Computing constrained centroidal diagram..." << std::endl;
  std::cout << "  Iteration,  ||g||,   max(dx)" << std::endl;

  for(unsigned int i=0; i<maxIts; ++i){
    SDOT_ASSERT(ot);

    std::tie(prices, dualObj) = ot->Solve(prices, opts);

    Eigen::MatrixXd pointGrad = ot->PointGradient();
    Eigen::Map<Eigen::VectorXd> pointGradVec(pointGrad.data(),2*pointGrad.cols());

    double stepSize = 1.0;
    Eigen::Matrix2Xd dir = -pointGrad.array()*(ot->LloydPointHessian()+1e-12*Eigen::MatrixXd::Ones(2,numPts)).array().inverse();

    newPts = pts + stepSize*dir;

    // Backtrack until all the points are  in the domain
    while((newPts.row(0).minCoeff()<ot->grid->xMin)||(newPts.row(0).maxCoeff()>ot->grid->xMax)||(newPts.row(1).minCoeff()<ot->grid->yMin)||(newPts.row(1).maxCoeff()>ot->grid->yMax)){
      stepSize *= 0.75;
      newPts = pts + stepSize*dir;
    }

    resid = (newPts - pts).cwiseAbs().maxCoeff();;
    std::printf("  %9d, %5.3e, %5.3e\n", i, pointGrad.norm(), resid);

    if(resid<xtol){
      std::cout << "Converged due to small stepsize." << std::endl;
      return ot->Diagram();
    }

    if((pointGradVec.norm()/(2*numPts))<gtol){
      std::cout << "Converged due to small gradient." << std::endl;
      return ot->Diagram();
    }

    pts = newPts;
    ot->SetPoints(pts);
  }

  std::cout << "WARNING: Did not converge to constrained centroidal diagram." << std::endl;
  return ot->Diagram();
}

template<typename ConjugateFunctionType>
std::shared_ptr<LaguerreDiagram> SemidiscreteOT<ConjugateFunctionType>::BuildCentroidal(std::shared_ptr<Distribution2d> const& dist,
                                                                 Eigen::VectorXd                 const& probs,
                                                                 OptionList                             opts)
{
  BoundingBox bbox(dist->Grid()->xMin, dist->Grid()->xMax, dist->Grid()->yMin, dist->Grid()->yMax);
  Eigen::Matrix2Xd initialPts = LaguerreDiagram::LatinHypercubeSample(bbox, probs.size());
  return BuildCentroidal(dist, initialPts, probs, opts);
}

template<typename ConjugateFunctionType>
std::shared_ptr<LaguerreDiagram> SemidiscreteOT<ConjugateFunctionType>::BuildCentroidal(std::shared_ptr<Distribution2d> const& dist,
                                                                 unsigned int                           numPts,
                                                                 OptionList                             opts)
{
  Eigen::VectorXd probs = (1.0/numPts)*Eigen::VectorXd::Ones(numPts);
  return BuildCentroidal(dist, probs, opts);
}

template<typename ConjugateFunctionType>
Eigen::Matrix2Xd SemidiscreteOT<ConjugateFunctionType>::MarginalCentroids(Eigen::VectorXd const& prices,
                                                                          LaguerreDiagram const& lagDiag) const
{
  int numCells =  lagDiag.NumCells();
  Eigen::Matrix2Xd centroids(2,numCells);

  // Loop over all ofthe cells
  for(int i=0; i<numCells; ++i){

    auto triFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2, Eigen::Vector2d const& pt3)
      {
        return ConjugateFunctionType::TriangularIntegralMarginalCentroid(prices(i), discrPts.col(i), pt1,pt2,pt3,unbalancedPenalty);
      };
    auto rectFunc = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
      {
        return ConjugateFunctionType::RectangularIntegralMarginalCentroid(prices(i), discrPts.col(i), pt1,pt2,unbalancedPenalty);
      };


    centroids.col(i) = lagDiag.IntegrateOverCell(i, triFunc, rectFunc, dist);


    auto triFunc2 = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2, Eigen::Vector2d const& pt3)
      {
        return ConjugateFunctionType::TriangularIntegralMarginalMass(prices(i), discrPts.col(i), pt1,pt2,pt3,unbalancedPenalty);
      };
    auto rectFunc2 = [&](Eigen::Vector2d const& pt1, Eigen::Vector2d const& pt2)
      {
        return ConjugateFunctionType::RectangularIntegralMarginalMass(prices(i), discrPts.col(i), pt1,pt2,unbalancedPenalty);
      };

    centroids.col(i) /= lagDiag.IntegrateOverCell(i, triFunc2, rectFunc2, dist);
  }

  return centroids;
}

template<>
Eigen::Matrix2Xd SemidiscreteOT<sdot::distances::Wasserstein2>::MarginalCentroids(Eigen::VectorXd const& prices,
                                                                                  LaguerreDiagram const& lagDiag) const
{
  return lagDiag.Centroids(dist);
}




namespace sdot{
  template class SemidiscreteOT<sdot::distances::Wasserstein2>;
  template class SemidiscreteOT<sdot::distances::QuadraticRegularization>;
  template class SemidiscreteOT<sdot::distances::GHK>;
}
