
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
  int numPts = 4;
  Eigen::VectorXd costs = Eigen::VectorXd::Ones(numPts);

  Eigen::Matrix2Xd pts(2,numPts);
  pts <<  0.680375,  0.566198,  0.823295,  0.329554,//  0.444451, 0.0452059,  0.270431,  0.904459,  0.271423,  0.716795,
          0.211234,   0.59688,  0.604897,  0.536459;//,   0.10794,  0.257742, 0.0268018,   0.83239,  0.434594,  0.213938;


  std::cout << "Points = \n";
  std::cout << "[[" << pts(0,0) << "," << pts(1,0) << "]";
  for(int i=1; i<numPts; ++i){
    std::cout << ", [" << pts(0,i) << "," << pts(1,i) << "]";
  }
  std::cout << "]" << std::endl;

  Eigen::Matrix2Xd domain(2,4);
  domain << 0.0, 1.0, 1.0, 0.0,
            0.0, 0.0, 1.0, 1.0;
  // domain << 0.0, 1.0, 1.0, 0.0,
  //           0.0, 0.0, 1.0, 1.0;

  LaguerreDiagram diag(domain(0,0), domain(0,1), domain(1,0), domain(1,2), pts, costs);

  auto grid = std::make_shared<RegularGrid>(domain(0,0),domain(1,0), domain(0,2), domain(1,2), 10, 10);

  double area = 0.0;
  Eigen::VectorXd localAreas = Eigen::VectorXd::Zero(numPts);
  Eigen::MatrixXd cellAreas = Eigen::MatrixXd::Zero(grid->NumCells(0), grid->NumCells(1));

  for(int polyInd=0; polyInd<numPts; ++polyInd){
    std::cout << "\n\n==================================\n";
    std::cout << "Polygon " << polyInd << std::endl;
    std::shared_ptr<PolygonRasterizeIter::Polygon_2> poly = diag.GetCell(polyInd);

    auto vertIt = poly->vertices_begin();
    std::cout << "[[" << vertIt->x() << "," << vertIt->y() << "]";
    vertIt++;
    for(;  vertIt != poly->vertices_end(); ++vertIt){
      std::cout << ", [" << vertIt->x() << "," << vertIt->y() << "]";
    }
    std::cout << "]" << std::endl;

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

      std::cout << "    area(" << gridIter.Indices().first << "," << gridIter.Indices().second << ")= " << cellArea << std::endl;
      cellAreas(gridIter.Indices().first, gridIter.Indices().second) += cellArea;
      localAreas(polyInd) += std::abs(cellArea);
      area += std::abs(cellArea);

      gridIter.Increment();
    }

    std::cout << "  Polygon area = " << localAreas(polyInd) << std::endl;
    std::cout << "  Polygon area error = " << localAreas(polyInd) - CGAL::to_double( poly->area() ) << std::endl;
    // std::cout << "\n\n";
  }

  // Eigen::VectorXd trueAreas(3);
  // trueAreas << 0.15*0.15, (domain(0,2)-0.15)*0.15 + 0.5*(domain(1,2)-0.15)*(domain(0,2)-0.15), 0.15*(domain(1,2)-0.15) + 0.5*(domain(1,2)-0.15)*(domain(0,2)-0.15);

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
