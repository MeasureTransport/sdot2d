
#include <Eigen/Core>

#include <vector>
#include <memory>

#include <CGAL/Boolean_set_operations_2.h>

// standard includes
#include <iostream>

#include "SDOT/PolygonRasterize.h"
#include "SDOT/RegularGrid.h"
#include "SDOT/LaguerreDiagram.h"

using namespace sdot;

int main(int argc, char* argv[])
{
  int numPts = 3;
  Eigen::VectorXd costs(numPts);
  costs << 1.0, 1.00751, 1.00751;

  Eigen::Matrix2Xd pts(2,numPts);
  pts << 0.1, 0.2, 0.1,
         0.1, 0.1, 0.2;

  Eigen::Matrix2Xd domain(2,4);
  domain << 0.0, 0.2, 0.2, 0.0,
            0.0, 0.0, 0.2, 0.2;
  // domain << 0.0, 1.0, 1.0, 0.0,
  //           0.0, 0.0, 1.0, 1.0;

  LaguerreDiagram diag(domain(0,0), domain(0,1), domain(1,0), domain(1,2), pts, costs);

  auto grid = std::make_shared<RegularGrid>(domain(0,0),domain(1,0), domain(0,2), domain(1,2), 50, 50);

  double area = 0.0;
  Eigen::VectorXd localAreas = Eigen::VectorXd::Zero(3);
  Eigen::MatrixXd cellAreas = Eigen::MatrixXd::Zero(grid->NumCells(0), grid->NumCells(1));

  for(int polyInd=0; polyInd<3; ++polyInd){
    std::cout << "\n\n==================================\n";
    std::cout << "Polygon " << polyInd << std::endl;
    std::shared_ptr<PolygonRasterizeIter::Polygon_2> poly = diag.GetCell(polyInd);

    std::cout << "  " << *poly << std::endl;
    PolygonRasterizeIter gridIter(grid,poly);
    //
    // unsigned int oldYInd = gridIter.Indices().second;
    // std::cout << "yind = " << oldYInd << std::endl;
    // std::cout << "yval = " << grid->yMin + oldYInd*grid->dy << std::endl;
    // std::cout << "    ";

    while(gridIter.IsValid()){


      // if(gridIter.Indices().second != oldYInd){
      //   oldYInd = gridIter.Indices().second;
      //   std::cout << "\nyind = " << oldYInd << std::endl;
      //   std::cout << "yval = " << grid->yMin + oldYInd*grid->dy << std::endl;
      //   std::cout << "    ";
      // }

      double cellArea;


      if(gridIter.IsBoundary()){
        cellArea = CGAL::to_double( gridIter.OverlapPoly()->area() );
        //std::cout << "  " << gridIter.Indices().first << ", " << gridIter.Indices().second << " -> Overlap poly = " << *gridIter.OverlapPoly() << std::endl;
        // unsigned int indX = gridIter.Indices().first;
        // unsigned int indY = gridIter.Indices().second;
        //
        // PolygonRasterizeIter::Polygon_2 tempPoly;
        // tempPoly.push_back(PolygonRasterizeIter::Point_2(grid->xMin + indX*grid->dx, grid->yMin+indY*grid->dy));
        // tempPoly.push_back(PolygonRasterizeIter::Point_2(grid->xMin + (indX+1)*grid->dx, grid->yMin+indY*grid->dy));
        // tempPoly.push_back(PolygonRasterizeIter::Point_2(grid->xMin + (indX+1)*grid->dx, grid->yMin+(indY+1)*grid->dy));
        // tempPoly.push_back(PolygonRasterizeIter::Point_2(grid->xMin + indX*grid->dx, grid->yMin+(indY+1)*grid->dy));
        //
        // // Compute the intersection of P and Q.
        // std::list<PolygonRasterizeIter::Polygon_with_holes_2> interList;
        // CGAL::intersection(*poly, tempPoly, std::back_inserter(interList));
        //
        // double trueArea = CGAL::to_double( interList.begin()->outer_boundary().area() );
        //
        // //std::cout << *gridIter.OverlapPoly() << std::endl;
        // double error = cellArea-trueArea;
        // if(std::abs(error)>std::numeric_limits<double>::epsilon()){
        //   // std::cout << "  Cell area = " << cellArea << " with error " << error << std::endl;
        // }
      }else{
        cellArea = grid->dx*grid->dy;
      }

      cellAreas(gridIter.Indices().first, gridIter.Indices().second) += cellArea;
      localAreas(polyInd) += std::abs(cellArea);
      area += std::abs(cellArea);

      gridIter.Increment();
    }
    // std::cout << "\n\n";
  }

  Eigen::VectorXd trueAreas(3);
  trueAreas << 0.15*0.15, (domain(0,2)-0.15)*0.15 + 0.5*(domain(1,2)-0.15)*(domain(0,2)-0.15), 0.15*(domain(1,2)-0.15) + 0.5*(domain(1,2)-0.15)*(domain(0,2)-0.15);

  // List all cells with errors
  for(unsigned int yInd=0; yInd<grid->NumCells(1); ++yInd){
    for(unsigned int xInd=0; xInd<grid->NumCells(0); ++xInd){
      if(std::abs(cellAreas(xInd,yInd) - grid->dx*grid->dy)>1e-15){
        std::cout << "Cell " << xInd << "," << yInd << " error = " << cellAreas(xInd,yInd) - grid->dx*grid->dy << std::endl;
      }
    }
  }

  std::cout << "Local areas = " << localAreas.transpose() << std::endl;

  //std::cout << "dx*dy = " << grid->dx * grid->dy << std::endl;
  std::cout << "True total area = " << grid->dx*grid->dy*grid->NumCells() << std::endl;
  std::cout << "Total Area = " << area << std::endl;
  std::cout << "Error = " << area - grid->dx*grid->dy*grid->NumCells() << std::endl;
  return 0;
}
