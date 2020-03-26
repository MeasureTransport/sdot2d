#ifndef SEMIDISCRETEOT_H_
#define SEMIDISCRETEOT_H_

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

    double Objective(Eigen::VectorXd const& prices);

    // Eigen::VectorXd Gradient(Eigen::VectorXd const& prices);
    //
    // Eigen::MatrixXd Hessian(Eigen::VectorXd const& prices);

  private:

    std::shared_ptr<Distribution2d> dist;
    std::shared_ptr<RegularGrid> grid;

    Eigen::Matrix2Xd discrPts;
    Eigen::VectorXd  discrProbs;

    Eigen::Matrix2Xd domain; // <- Points defining a polygon surrounding the domain of interest


    static double SquareIntegral(double xmin, double xmax,
                                 double ymin, double ymax,
                                 double px,   double py);

    static double TriangleIntegral(double x1, double y1,
                                   double x2, double y2,
                                   double x3, double y3,
                                   double px, double py);

  };




  SemidiscreteOT::SemidiscreteOT(std::shared_ptr<Distribution2d> const& distIn,
                                 Eigen::Matrix2Xd                const& discrPtsIn,
                                 Eigen::VectorXd                 const& discrProbsIn) : dist(distIn),
                                                                                        grid(distIn->Grid()),
                                                                                        discrPts(discrPtsIn),
                                                                                        discrProbs(discrProbsIn)
  {
    domain.resize(2,4);
    domain << grid->xMin, grid->xMax, grid->xMax, grid->xMin,
              grid->yMin, grid->yMin, grid->yMax, grid->yMax;

    assert(discrPtsIn.cols()==discrProbsIn.size());
  }



  double SemidiscreteOT::Objective(Eigen::VectorXd const& prices)
  {
      // Notes:
      //   - The cost c(x,y) is the squared distance between x and y
      //   - See (17) of https://arxiv.org/pdf/1710.02634.pdf

      const int numCells = discrPts.cols();
      assert(numCells==prices.size());

      // Construct the Laguerre diagram
      LaguerreDiagram lagDiag(discrPts, prices, domain);

      // Holds the part of the objective for each cell in the Laguerre diagram
      Eigen::VectorXd cellParts(numCells);
      Eigen::VectorXd cellAreas = Eigen::VectorXd::Zero(numCells);

      for(int cellInd=0; cellInd<numCells; ++cellInd){

        cellParts(cellInd) = prices(cellInd)*discrProbs(cellInd);

        // Loop over the grid cells in this Laguerre cell
        auto lagCell = lagDiag.GetCell(cellInd);

        PolygonRasterizeIter gridIter(grid,lagCell);

        unsigned int xInd, yInd;

        while(gridIter.IsValid()){

          xInd = gridIter.Indices().first;
          yInd = gridIter.Indices().second;

          // The probability in this grid cell
          double cellProb = dist->Probability(xInd,yInd);

          cellParts(cellInd) -= cellProb*prices(cellInd);

          if(gridIter.IsBoundary()){

            // Break the intersection polygon into triangles and add contributions from each triangle
            std::shared_ptr<PolygonRasterizeIter::Polygon_2> overlapPoly = gridIter.OverlapPoly();
            cellAreas(cellInd) += CGAL::to_double( overlapPoly->area() );

            auto beginVert = overlapPoly->vertices_begin();
            auto vert1 = beginVert;
            vert1++;
            auto vert2 = vert1;
            vert2++;

            const double x1 = CGAL::to_double( beginVert->x() );
            const double y1 = CGAL::to_double( beginVert->y() );

            for(; vert2!=overlapPoly->vertices_end(); vert2++, vert1++)
            {
              cellParts(cellInd) += cellProb * TriangleIntegral(x1,                            y1,
                                                                CGAL::to_double( vert1->x() ), CGAL::to_double( vert1->y() ),
                                                                CGAL::to_double( vert2->x() ), CGAL::to_double( vert2->y() ),
                                                                discrPts(0,cellInd),           discrPts(1,cellInd));
            }

          }else{
            cellAreas(cellInd) += grid->dx*grid->dy;

            cellParts(cellInd) += cellProb*SquareIntegral(gridIter.LeftX(),    gridIter.RightX(),
                                                          gridIter.BottomY(),    gridIter.TopY(),
                                                          discrPts(0,cellInd), discrPts(1,cellInd));
          }

          gridIter.Increment();
        }
      }

      std::cout << "Cell Parts = " << cellParts.transpose() << std::endl;
      std::cout << "Cell Areas = " << cellAreas.transpose() << std::endl;
      std::cout << "Total area = " << cellAreas.sum() << std::endl;
      std::cout << "Area error = " << cellAreas.sum() - (grid->xMax - grid->xMin)*(grid->yMax - grid->yMin) << std::endl;
      return cellParts.sum();
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

}



#endif // #ifndef SEMIDISCRETEOT_H_
