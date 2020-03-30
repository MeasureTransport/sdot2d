#include "SDOT/SemidiscreteOT.h"

#include <algorithm>

using namespace sdot;


SemidiscreteOT::SemidiscreteOT(std::shared_ptr<Distribution2d> const& distIn,
                               Eigen::Matrix2Xd                const& discrPtsIn,
                               Eigen::VectorXd                 const& discrProbsIn) : dist(distIn),
                                                                                      grid(distIn->Grid()),
                                                                                      discrPts(discrPtsIn),
                                                                                      discrProbs(discrProbsIn)
{
  // domain.resize(2,4);
  // domain << grid->xMin, grid->xMax, grid->xMax, grid->xMin,
  //           grid->yMin, grid->yMin, grid->yMax, grid->yMax;

  assert(discrPtsIn.cols()==discrProbsIn.size());
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


std::pair<double,Eigen::VectorXd> SemidiscreteOT::ComputeGradient(Eigen::VectorXd const& prices,
                                                                  LaguerreDiagram const& lagDiag) const
{
 const int numCells = prices.size();

 // Holds the part of the objective for each cell in the Laguerre diagram
 Eigen::VectorXd objParts(numCells);
 Eigen::VectorXd gradient(numCells);

 //double totalArea = 0.0;

 for(int cellInd=0; cellInd<numCells; ++cellInd){

   //double lagCellArea = 0.0;

   objParts(cellInd) = prices(cellInd)*discrProbs(cellInd);
   gradient(cellInd) = discrProbs(cellInd);

   // Loop over the grid cells in this Laguerre cell
   auto lagCell = lagDiag.GetCell(cellInd);

   PolygonRasterizeIter gridIter(grid,lagCell);

   unsigned int xInd, yInd;
   double x1, x2, x3, y1, y2, y3, gridCellDens;

   while(gridIter.IsValid()){

     xInd = gridIter.Indices().first;
     yInd = gridIter.Indices().second;

     // The probability in this grid cell
     gridCellDens = dist->Density(xInd,yInd);


     //double interArea = 0.0;
     if(gridIter.IsBoundary()){


       // Break the intersection polygon into triangles and add contributions from each triangle
       std::shared_ptr<PolygonRasterizeIter::Polygon_2> overlapPoly = gridIter.OverlapPoly();

       auto beginVert = overlapPoly->vertices_begin();
       auto vert1 = beginVert;
       vert1++;
       auto vert2 = vert1;
       vert2++;

       x1 = CGAL::to_double( beginVert->x() );
       y1 = CGAL::to_double( beginVert->y() );

       //interArea = CGAL::to_double( overlapPoly->area() );

       for(; vert2!=overlapPoly->vertices_end(); vert2++, vert1++)
       {
         x2 = CGAL::to_double( vert1->x() );
         y2 = CGAL::to_double( vert1->y() );
         x3 = CGAL::to_double( vert2->x() );
         y3 = CGAL::to_double( vert2->y() );


         objParts(cellInd) += gridCellDens * TriangleIntegral(x1,                  y1,
                                                              x2,                  y2,
                                                              x3,                  y3,
                                                              discrPts(0,cellInd), discrPts(1,cellInd));

         // gridCellDense * triArea, where triArea=0.5*std::abs((x2*y1-x1*y2)+(x3*y2-x2*y3)+(x1*y3-x3*y1))
         double triArea = std::abs((x2*y1-x1*y2)+(x3*y2-x2*y3)+(x1*y3-x3*y1));
         //interArea += triArea;
         objParts(cellInd) -= gridCellDens*triArea*prices(cellInd);
         gradient(cellInd) -= gridCellDens*triArea;
       }

     }else{
       //interArea += grid->dx*grid->dy;

       objParts(cellInd) -= gridCellDens*grid->dx*grid->dy;
       objParts(cellInd) += gridCellDens*SquareIntegral(gridIter.LeftX(),    gridIter.RightX(),
                                                         gridIter.BottomY(),  gridIter.TopY(),
                                                         discrPts(0,cellInd), discrPts(1,cellInd));

       // gridCellDens * rectArea, where rectArea = grid->dx*grid->dy
       gradient(cellInd) -= gridCellDens*grid->dx*grid->dy;
     }

     //lagCellArea += interArea;

     gridIter.Increment();
   }

   //totalArea += lagCellArea;
 }

 return std::make_pair(objParts.sum(), gradient);
}

double SemidiscreteOT::LineIntegral(LaguerreDiagram::Point_2 const& srcPt,
                                    LaguerreDiagram::Point_2 const& tgtPt) const
{
  const double compTol = 1e-14;

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

    while(nextY<maxY){
      val += (nextY-currY);//*dist->Density(xInd,yInd);

      yInd++;
      currY = nextY;
      nextY = currY+grid->dy;
    }

    val += (maxY-currY);//*dist->Density(xInd,yInd);

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

    while(nextX<maxX){
      val += (nextX-currX);//*dist->Density(xInd,yInd);

      xInd++;
      currX = nextX;
      nextX = currX+grid->dx;
    }

    val += (maxX-currX);//*dist->Density(xInd,yInd);

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

    nextt_vert = std::min(1.0, ( (xInd+1)*grid->dx + grid->xMin - xs) / dx);

    if(posSlope){
      nextt_horz = std::min(1.0, ( (yInd+1)*grid->dy + grid->yMin - ys) / dy);
    }else{
      nextt_horz = std::min(1.0, ( (yInd-1)*grid->dy + grid->yMin - ys) / dy);
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
          nextt_horz = std::min(1.0, ((yInd-1)*grid->dy + grid->yMin-ys)/dy);
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
  hessVals.push_back(T(0.0,0.0,1e8));

  double intVal;
  unsigned int cellInd2;
  LaguerreDiagram::Point_2 srcPt, tgtPt;

  for(unsigned int cellInd1=0; cellInd1<numCells; ++cellInd1){
    for(auto edgeTuple : lagDiag.InternalEdges(cellInd1)){
      std::tie(cellInd2, srcPt, tgtPt) = edgeTuple;

      // Compute the integral of the target density along the edge
      intVal = LineIntegral(srcPt,tgtPt)/(discrPts.col(cellInd1)-discrPts.col(cellInd2)).norm();

      diagVals(cellInd1) -= intVal;
      hessVals.push_back(T(cellInd1,cellInd2,intVal));
    }
  }

  //
  // // Get a pointer to the underlying CGAL PowerDiagram so we can access the half edges
  // LaguerreDiagram::PowerDiagram* baseDiag = lagDiag.BaseDiagram();
  //
  // // Loop over all the half edges in the power diagram
  // for(auto edge = baseDiag->halfedges_begin(); edge!=baseDiag->halfedges_end(); ++edge){
  //
  //
  // }

  for(int i=0; i<numCells; ++i)
    hessVals.push_back(T(i,i,diagVals(i)));

  Eigen::SparseMatrix<double> hess(numCells,numCells);
  hess.setFromTriplets(hessVals.begin(), hessVals.end());

  return hess;
}

double SemidiscreteOT::SquareIntegral(double xmin, double xmax,
                                      double ymin, double ymax,
                                      double px,   double py)
{
 double rectInt = (0.5/3.0)*(ymax-ymin)*(std::pow(xmax-px,3.0)-std::pow(xmin-px,3.0))
                + (0.5/3.0)*(xmax-xmin)*(std::pow(ymax-py,3.0)-std::pow(ymin-py,3.0));

 return rectInt;
}

double SemidiscreteOT::TriangleIntegral(double x1, double y1,
                                       double x2, double y2,
                                       double x3, double y3,
                                       double px, double py)
{
 double triInt = (1.0/2.0)*std::pow(px, 2) - 1.0/3.0*px*x1 - 1.0/3.0*px*x2 - 1.0/3.0*px*x3 + (1.0/2.0)*std::pow(py, 2) - 1.0/3.0*py*y1 - 1.0/3.0*py*y2 - 1.0/3.0*py*y3 + (1.0/12.0)*std::pow(x1, 2) + (1.0/12.0)*x1*x2 + (1.0/12.0)*x1*x3 + (1.0/12.0)*std::pow(x2, 2) + (1.0/12.0)*x2*x3 + (1.0/12.0)*std::pow(x3, 2) + (1.0/12.0)*std::pow(y1, 2) + (1.0/12.0)*y1*y2 + (1.0/12.0)*y1*y3 + (1.0/12.0)*std::pow(y2, 2) + (1.0/12.0)*y2*y3 + (1.0/12.0)*std::pow(y3, 2);
 triInt *= 0.5*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));

 return triInt;
}
