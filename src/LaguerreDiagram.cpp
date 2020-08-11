#include "SDOT/LaguerreDiagram.h"

#include "SDOT/PolygonRasterize.h"

#include  "SDOT/Integrands/ConstantIntegrand.h"

// standard includes
#include <iostream>
#include <random>
#include <chrono>

// Additional CGAL includes
#include <CGAL/squared_distance_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Boolean_set_operations_2.h>

using namespace sdot;

LaguerreDiagram::LaguerreDiagram(double xBndLeftIn,   double xBndRightIn,
                                 double yBndBottomIn, double yBndTopIn,
                                 Eigen::Matrix2Xd const& pts,
                                 Eigen::VectorXd  const& costs) : bbox(xBndLeftIn, xBndRightIn, yBndBottomIn, yBndTopIn),
                                                                  infDist(2.0*std::sqrt(std::pow(xBndRightIn-xBndLeftIn,2.0) + std::pow(yBndTopIn-yBndBottomIn,2.0)))
{
  numPts = pts.cols();
  assert(costs.size()==numPts);

  CreateUnboundedDiagram(pts,costs);
  CreateBoundedCells(pts);
}


LaguerreDiagram::LaguerreDiagram(BoundingBox const& bbox,
                                 Eigen::Matrix2Xd  const& pts,
                                 Eigen::VectorXd  const& costs) : LaguerreDiagram(bbox.xMin, bbox.xMax, bbox.yMin, bbox.yMax, pts, costs)
{
}


Eigen::Matrix2Xd LaguerreDiagram::GetCellVertices(int ind) const
{
  auto& poly = laguerreCells.at(ind);
  Eigen::Matrix2Xd points(2,poly->size());

  // Copy the vertices from the polygon to the Eigen::Matrix
  unsigned int col = 0;
  for(auto vertIter = poly->vertices_begin(); vertIter!=poly->vertices_end(); ++vertIter, ++col){
    points(0,col) = CGAL::to_double( vertIter->x() );
    points(1,col) = CGAL::to_double( vertIter->y() );
  }

  return points;
}

Eigen::VectorXd LaguerreDiagram::Areas(std::shared_ptr<Distribution2d> const& dist) const
{
  Eigen::VectorXd areas(NumCells());
  for(int i=0; i<areas.size();++i)
    areas(i) = CellArea(i, dist);
  return areas;
}


double LaguerreDiagram::CellArea(unsigned int cellInd,
                                 std::shared_ptr<Distribution2d> const& dist) const
{
  auto integrand = std::make_shared<ConstantIntegrand>(1.0);
  return IntegrateOverCell(cellInd, integrand, dist);
  //
  // double x1, x2, y1, y2, x3, y3;
  //
  // auto beginVert = laguerreCells.at(cellInd)->vertices_begin();
  // auto vert1 = beginVert;
  // vert1++;
  // auto vert2 = vert1;
  // vert2++;
  //
  // x1 = CGAL::to_double( beginVert->x() );
  // y1 = CGAL::to_double( beginVert->y() );
  //
  // //interArea = gridCellDens*CGAL::to_double( overlapPoly->area() );
  //
  // double polyArea = 0.0;
  // double triArea;
  //
  // for(; vert2!=laguerreCells.at(cellInd)->vertices_end(); vert2++, vert1++)
  // {
  //   x2 = CGAL::to_double( vert1->x() );
  //   y2 = CGAL::to_double( vert1->y() );
  //   x3 = CGAL::to_double( vert2->x() );
  //   y3 = CGAL::to_double( vert2->y() );
  //
  //   triArea = 0.5*std::abs((x2*y1-x1*y2)+(x3*y2-x2*y3)+(x1*y3-x3*y1));
  //   polyArea += triArea;
  // }
  //
  // return polyArea;
}

double LaguerreDiagram::IntegrateOverCell(unsigned int                     cellInd,
                         std::shared_ptr<Integrand> const& integrand) const
{
  Eigen::Vector2d pt1(2), pt2(2), pt3(2);

  auto beginVert = laguerreCells.at(cellInd)->vertices_begin();
  auto vert1 = beginVert;
  vert1++;
  auto vert2 = vert1;
  vert2++;

  pt1(0) = CGAL::to_double( beginVert->x() );
  pt1(0) = CGAL::to_double( beginVert->y() );

  //interArea = gridCellDens*CGAL::to_double( overlapPoly->area() );

  double result = 0.0;
  double triArea;

  for(; vert2!=laguerreCells.at(cellInd)->vertices_end(); vert2++, vert1++)
  {
    pt2(0) = CGAL::to_double( vert1->x() );
    pt2(1) = CGAL::to_double( vert1->y() );
    pt3(0) = CGAL::to_double( vert2->x() );
    pt3(1) = CGAL::to_double( vert2->y() );

    result += integrand->TriangularIntegral(pt1,pt2,pt3);
  }

  return result;
}

double LaguerreDiagram::IntegrateOverCell(unsigned int                           cellInd,
                                          std::shared_ptr<Integrand>       const& integrand,
                                          std::shared_ptr<Distribution2d> const& dist) const
{

  if(dist==nullptr)
    return IntegrateOverCell(cellInd, integrand);

  auto& grid = dist->Grid();

  double result = 0.0;

  // Loop over the grid cells in this Laguerre cell
  if(laguerreCells.at(cellInd)->size()==0)
    return 0.0;

  PolygonRasterizeIter gridIter(dist->Grid(), laguerreCells.at(cellInd));

  unsigned int xInd, yInd;
  Eigen::Vector2d pt1(2), pt2(2), pt3(2);
  double gridCellDens;

  while(gridIter.IsValid()){

    xInd = gridIter.Indices().first;
    yInd = gridIter.Indices().second;

    // The probability in this grid cell
    gridCellDens = dist->Density(xInd,yInd);

    if(gridCellDens>0.0){

      // The area of the intersection of this grid cell and the Laguerre cell
      double interResult = 0.0;

      if(gridIter.IsBoundary()){

        // Break the intersection polygon into triangles and add contributions from each triangle
        std::shared_ptr<PolygonRasterizeIter::Polygon_2> overlapPoly = gridIter.OverlapPoly();
        assert(overlapPoly);
        assert(overlapPoly->size()>2); // <- Makes sure there is at least 3 nodes in the polygon

        auto beginVert = overlapPoly->vertices_begin();
        auto vert1 = beginVert;
        vert1++;
        auto vert2 = vert1;
        vert2++;

        pt1[0] = CGAL::to_double( beginVert->x() );
        pt1[1] = CGAL::to_double( beginVert->y() );

        for(; vert2!=overlapPoly->vertices_end(); vert2++, vert1++)
        {
          pt2[0] = CGAL::to_double( vert1->x() );
          pt2[1] = CGAL::to_double( vert1->y() );
          pt3[0] = CGAL::to_double( vert2->x() );
          pt3[1] = CGAL::to_double( vert2->y() );

          interResult += integrand->TriangularIntegral(pt1,pt2,pt3);
        }

      }else{
        pt1[0] = grid->xMin + grid->dx*xInd;
        pt1[1] = grid->yMin + grid->dy*yInd;
        pt2[0] = grid->xMin + grid->dx*(xInd+1);
        pt2[1] = grid->xMin + grid->dy*(yInd+1);

        interResult += integrand->RectangularIntegral(pt1,pt2);

      }

      result += interResult*gridCellDens;
    }
    gridIter.Increment();
  }

  return result;
}


Eigen::Matrix2Xd LaguerreDiagram::Centroids(std::shared_ptr<Distribution2d> const& dist) const
{
  Eigen::Matrix2Xd centroids(2,NumCells());

  for(int i=0; i<NumCells(); ++i){
    centroids.col(i) = CellCentroid(i, dist);
  }

  return centroids;
}


Eigen::Matrix2Xd LaguerreDiagram::SeedPts() const
{
  Eigen::Matrix2Xd pts(2,NumCells());

  for(auto faceIter = unboundedDiagram.faces_begin(); faceIter!=unboundedDiagram.faces_end(); ++faceIter){
    auto pt = faceIter->dual();// ->info()) = *faceIter;

    unsigned int ind = pt->info();
    pts(0,ind) = CGAL::to_double(pt->point().x());
    pts(1,ind) = CGAL::to_double(pt->point().y());
  }

  return pts;
}

Eigen::Vector2d LaguerreDiagram::CellCentroid(unsigned int cellInd) const
{
  double x1, x2, y1, y2, x3, y3;

  auto beginVert = laguerreCells.at(cellInd)->vertices_begin();
  auto vert1 = beginVert;
  vert1++;
  auto vert2 = vert1;
  vert2++;

  x1 = CGAL::to_double( beginVert->x() );
  y1 = CGAL::to_double( beginVert->y() );

  //interArea = gridCellDens*CGAL::to_double( overlapPoly->area() );
  double polyArea = 0.0;
  double triArea;
  Eigen::Vector2d centroid = Eigen::Vector2d::Zero(2);

  for(; vert2!=laguerreCells.at(cellInd)->vertices_end(); vert2++, vert1++)
  {
    x2 = CGAL::to_double( vert1->x() );
    y2 = CGAL::to_double( vert1->y() );
    x3 = CGAL::to_double( vert2->x() );
    y3 = CGAL::to_double( vert2->y() );

    triArea = 0.5*std::abs((x2*y1-x1*y2)+(x3*y2-x2*y3)+(x1*y3-x3*y1));
    polyArea += triArea;

    centroid(0) += triArea*(x1+x2+x3)/3.0;
    centroid(1) += triArea*(y1+y2+y3)/3.0;
  }

  return centroid / polyArea;
}


std::shared_ptr<LaguerreDiagram> LaguerreDiagram::BuildCentroidal(BoundingBox       const& bbox,
                                                                  Eigen::Matrix2Xd  const& initialPts,
                                                                  Eigen::VectorXd   const& prices,
                                                                  unsigned int maxIts,
                                                                  double tol,
                                                                  std::shared_ptr<Distribution2d> const& dist)
{
  assert(prices.size()==initialPts.cols());

  double resid;

  Eigen::Matrix2Xd currPts;
  Eigen::Matrix2Xd newPts = initialPts;

  std::shared_ptr<LaguerreDiagram> diagram;

  std::cout << "Constructing Centroidal Power diagram with " << prices.size() << " points." << std::endl;

  for(unsigned int i=0; i<maxIts; ++i){
    std::swap(currPts,newPts);

    diagram = std::make_shared<LaguerreDiagram>(bbox, currPts, prices);
    newPts = diagram->Centroids(dist);

    resid = (newPts-currPts).cwiseAbs().maxCoeff();

    std::cout << "  After " << i << " iterations, residual = " << resid << std::endl;

    // Check to see if we've converged
    if(resid<tol){
      std::cout << "Converged after " << i << " iterations." << std::endl;
      return diagram;
    }
  }

  std::cout << "WARNING: Residual after " << maxIts << " is still " << resid << " which is larger than the specified tolerance of " << tol << "." << std::endl;
  return diagram;
}

Eigen::Matrix2Xd LaguerreDiagram::LatinHypercubeSample(BoundingBox const& bbox,
                                                       unsigned int       numPts)
{
  double dx = (bbox.xMax - bbox.xMin)/numPts;
  double dy = (bbox.yMax - bbox.yMin)/numPts;

  // Latin hypercube sampling
  std::vector<double> xInds(numPts);
  std::vector<double> yInds(numPts);
  for(unsigned int i=0; i<numPts; ++i){
    xInds[i] = i;
    yInds[i] = i;
  }

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine rg(seed);
  std::uniform_real_distribution<double> u(0.0,1.0);

  std::shuffle(xInds.begin(), xInds.end(), rg);
  std::shuffle(yInds.begin(), yInds.end(), rg);

  // Generate the points
  Eigen::Matrix2Xd initialPoints(2,numPts);
  for(unsigned int i=0; i<numPts; ++i){
    initialPoints(0,i) = bbox.xMin + dx*(u(rg)+xInds[i]);
    initialPoints(1,i) = bbox.yMin + dy*(u(rg)+yInds[i]);
  }

  return initialPoints;
}

std::shared_ptr<LaguerreDiagram> LaguerreDiagram::BuildCentroidal(BoundingBox const& bbox,
                                                                  unsigned int       numPts,
                                                                  unsigned int maxIts,
                                                                  double tol,
                                                                  std::shared_ptr<Distribution2d> const& dist)
{
  Eigen::Matrix2Xd initialPoints = LatinHypercubeSample(bbox,numPts);
  return BuildCentroidal(bbox, initialPoints, Eigen::VectorXd::Ones(numPts), maxIts, tol, dist);
}



Eigen::Vector2d LaguerreDiagram::CellCentroid(unsigned int cellInd, std::shared_ptr<Distribution2d> const& dist) const
{
  if(dist==nullptr)
    return CellCentroid(cellInd);

  auto& grid = dist->Grid();

  Eigen::Vector2d centroid = Eigen::Vector2d::Zero(2);

  double polyArea = 0.0;

  // Loop over the grid cells in this Laguerre cell
  PolygonRasterizeIter gridIter(dist->Grid(), laguerreCells.at(cellInd));

  unsigned int xInd, yInd;
  double x1, x2, x3, y1, y2, y3, gridCellDens;

  while(gridIter.IsValid()){

    xInd = gridIter.Indices().first;
    yInd = gridIter.Indices().second;

    // The probability in this grid cell
    gridCellDens = dist->Density(xInd,yInd);

    // The area of the intersection of this grid cell and the Laguerre cell
    //double interArea = 0.0;

    if(gridIter.IsBoundary()){

      // Break the intersection polygon into triangles and add contributions from each triangle
      std::shared_ptr<PolygonRasterizeIter::Polygon_2> overlapPoly = gridIter.OverlapPoly();
      assert(overlapPoly);
      assert(overlapPoly->size()>2); // <- Makes sure there is at least 3 nodes in the polygon

      auto beginVert = overlapPoly->vertices_begin();
      auto vert1 = beginVert;
      vert1++;
      auto vert2 = vert1;
      vert2++;

      x1 = CGAL::to_double( beginVert->x() );
      y1 = CGAL::to_double( beginVert->y() );

      //interArea = gridCellDens*CGAL::to_double( overlapPoly->area() );

      for(; vert2!=overlapPoly->vertices_end(); vert2++, vert1++)
      {
        x2 = CGAL::to_double( vert1->x() );
        y2 = CGAL::to_double( vert1->y() );
        x3 = CGAL::to_double( vert2->x() );
        y3 = CGAL::to_double( vert2->y() );

        double triArea = 0.5*std::abs((x2*y1-x1*y2)+(x3*y2-x2*y3)+(x1*y3-x3*y1));

        centroid(0) += gridCellDens*triArea * (x1+x2+x3)/3.0;
        centroid(1) += gridCellDens*triArea * (y1+y2+y3)/3.0;

        polyArea += gridCellDens*triArea;
      }

    }else{
      polyArea += gridCellDens*grid->dx*grid->dy;
      centroid += gridCellDens*(grid->dx*grid->dy)*grid->Center(xInd,yInd);
    }

    gridIter.Increment();
  }

  centroid /= polyArea;
  return centroid;
}


void LaguerreDiagram::CreateUnboundedDiagram(Eigen::Matrix2Xd const& pts,
                                             Eigen::VectorXd  const& costs)
{
  assert(costs.size()==numPts);

  //std::vector<Site_2> wpoints(numPts);
  for(unsigned int i=0; i<numPts; ++i){
    assert(costs(i)>=0);
    unboundedDiagram.insert( Site_2(Point_2(pts(0,i), pts(1,i)), std::sqrt(costs(i))) )->dual()->info()= i;
  }

  assert(unboundedDiagram.is_valid());
}

bool LaguerreDiagram::HasInternalSource(Ccb_halfedge_circulator halfEdge) const
{
  if(!halfEdge->has_source()){
    return false;
  }else{
    double x = CGAL::to_double(halfEdge->source()->point().x());
    double y = CGAL::to_double(halfEdge->source()->point().y());
    return (x<bbox.xMax+compTol)&&(x>bbox.xMin-compTol)&&(y<bbox.yMax+compTol)&&(y>bbox.yMin-compTol);
  }
}

bool LaguerreDiagram::HasInternalTarget(Ccb_halfedge_circulator halfEdge) const
{
  if(!halfEdge->has_target()){
    return false;
  }else{
    double x = CGAL::to_double(halfEdge->target()->point().x());
    double y = CGAL::to_double(halfEdge->target()->point().y());
    return (x<bbox.xMax+compTol)&&(x>bbox.xMin-compTol)&&(y<bbox.yMax+compTol)&&(y>bbox.yMin-compTol);
  }
}

std::pair<LaguerreDiagram::Point_2,LaguerreDiagram::Point_2> LaguerreDiagram::GetEdgeVerts(Ccb_halfedge_circulator& halfEdge)
{
  bool hasSrc = halfEdge->has_source();
  bool hasTgt = halfEdge->has_target();

  Point_2 midPt;
  Vector_2 dir;

  if((!hasSrc)||(!hasTgt)){
    // Get the node in the triangulation on the other side of this edge
    auto vt = halfEdge->twin()->face()->dual()->point();
    auto vs = halfEdge->face()->dual()->point();

    // Compute the direction of dividing line between this cell and the other cell
    dir = Vector_2(vs.y()-vt.y(), -(vs.x()-vt.x()));
    dir /= std::sqrt(CGAL::to_double(dir.squared_length()));

    // Find the weighted mid point between the triangulation sites
    double ts = CGAL::to_double( CGAL::squared_distance(vt.point(), vs.point()).exact() );
    double h = (0.5/ts)*CGAL::to_double(ts + vs.weight() - vt.weight());
    midPt = vs.point() + h*(vt.point()-vs.point());
  }

  Point_2 srcPt;
  if(hasSrc){
    srcPt = halfEdge->source()->point();
  } else {
    srcPt = midPt - infDist*dir;
  }

  Point_2 tgtPt;
  if(hasTgt){
    tgtPt = halfEdge->target()->point();
  }else{
    tgtPt = midPt + infDist*dir;
  }

  return std::make_pair(srcPt, tgtPt);
}


std::shared_ptr<LaguerreDiagram::Polygon_2> LaguerreDiagram::BoundOneCell(PowerDiagram::Face_handle const& face)
{
      // Get the point in the Regular triangulation at the center of this Laguerre cell
      auto vs = face->dual()->point();

      // Loop over the edges in the face
      Ccb_halfedge_circulator halfEdgeStart = face->ccb();
      Ccb_halfedge_circulator halfEdge = halfEdgeStart;

      // Construct a polygon from the face
      //std::vector<Point_2> polyPts;
      std::shared_ptr<Polygon_2> poly = std::make_shared<Polygon_2>();

      do {

        bool hasSrc = halfEdge->has_source();
        bool hasTgt = halfEdge->has_target();

        Point_2 srcPt, tgtPt;
        std::tie(srcPt, tgtPt) = GetEdgeVerts(halfEdge);

        // always add the source point
        poly->push_back(srcPt);

        // If the target was at infinity, then make sure to add the finite point
        if(!hasTgt)
          poly->push_back(tgtPt);

      } while ( ++halfEdge != halfEdgeStart);

      // Clip the polygon to the bounding box
      auto outPoly = bbox.ClipPolygon(poly);
      //
      // // Loop over all of the edges in the face
      // do {
      //
      //   bool hasSrc = halfEdge->has_source();
      //   bool hasTgt = halfEdge->has_target();
      //
      //   Point_2 srcPt, tgtPt;
      //   std::tie(srcPt, tgtPt) = GetEdgeVerts(halfEdge);
      //
      //   bool srcInside = bbox.IsInside(srcPt);
      //   bool tgtInside = bbox.IsInside(tgtPt);
      //
      //   if(srcInside){
      //     polyPts.push_back(srcPt);
      //
      //     // If the src is inside but the target is not, then we need to figure out where the edges cross
      //     if(!tgtInside){
      //
      //       bbox.ClipSegment(srcPt,tgtPt);
      //
      //       polyPts.push_back(tgtPt);
      //
      //       // Figure out the next edge that intersects the bounding box
      //       auto nextEdge = halfEdge;
      //       nextEdge++;
      //       std::tie(srcPt, tgtPt) = GetEdgeVerts(nextEdge);
      //       while(!bbox.ClipSegment(srcPt,tgtPt)){
      //         nextEdge++;
      //         std::tie(srcPt, tgtPt) = GetEdgeVerts(nextEdge);
      //       }
      //
      //       // Add any necessary corners and the point where the polygon reenters
      //       bbox.AddCorners(srcPt, polyPts);
      //       polyPts.push_back(srcPt);
      //     }
      //
      //   // source and target are outside
      //   }else{
      //
      //     if(bbox.ClipSegment(srcPt,tgtPt)){
      //
      //       polyPts.push_back(srcPt);
      //       // If the target is not inside, we might also have to add some corners
      //       if(!tgtInside){
      //         polyPts.push_back(tgtPt);
      //
      //         auto nextEdge = halfEdge;
      //         nextEdge++;
      //         std::tie(srcPt, tgtPt) = GetEdgeVerts(nextEdge);
      //         while(!bbox.ClipSegment(srcPt,tgtPt)){
      //           nextEdge++;
      //           std::tie(srcPt, tgtPt) = GetEdgeVerts(nextEdge);
      //         }
      //         bbox.AddCorners(srcPt, polyPts);
      //       }
      //     }
      //   }
      //
      // } while ( ++halfEdge != halfEdgeStart);
      //
      // auto outPoly = std::make_shared<Polygon_2>(polyPts.begin(), polyPts.end());//temp.begin()->outer_boundary());
      assert(outPoly);

      if(outPoly->size()<3){
        std::stringstream msg;
        msg << "Could not construct Laguerre diagram. Found empty Laguerre cell.\n";
        msg << "  Laguerre cell vertices:\n";

        for(auto it = outPoly->vertices_begin(); it!=outPoly->vertices_end(); ++it)
          msg << "    (" <<it->x() << "," << it->y() << ")\n";
        msg << "\n";
        throw LaguerreDiagram::ConstructionException(msg.str());
      }

      // Remove near duplicate vertices
      auto currVert = outPoly->vertices_begin();
      auto nextVert = currVert;
      nextVert++;
      while(nextVert != outPoly->vertices_end()){

        if( CGAL::squared_distance(*nextVert, *currVert)< compTol){
          nextVert = outPoly->erase(nextVert);
        }else{
          currVert = nextVert;
          nextVert++;
        }
      }

      if(CGAL::squared_distance(*outPoly->vertices_begin(), *currVert)< compTol){
        outPoly->erase(currVert);
      }

      // Sanity check to make sure there are no duplicate vertices or edge crossings
      if(!outPoly->is_simple()){
        std::cout << "Created a polygon that is not simple:\n";
        for(auto it = outPoly->vertices_begin(); it!=outPoly->vertices_end(); ++it){
          std::cout << "(" << CGAL::to_double(it->x()) << "," << CGAL::to_double(it->y()) << ") -> ";
        }
        std::cout << std::endl;
        assert(outPoly->is_simple());
      }

      return outPoly;
}


void LaguerreDiagram::CreateBoundedCells(Eigen::Matrix2Xd const& pts)
{
  internalEdges.clear();
  internalEdges.resize(numPts);

  // Figure out which faces correspond to the original points
  std::vector<Face_handle> faces(numPts);
  unsigned int faceCount = 0;
  for(auto faceIter = unboundedDiagram.faces_begin(); faceIter!=unboundedDiagram.faces_end(); ++faceIter){
    faces.at(faceIter->dual()->info()) = *faceIter;
    faceCount++;
  }

  if(faceCount!=numPts){
    std::stringstream msg;
    msg << "Could not construct Laguerre diagram. ";
    msg << "There should be " << numPts << " faces, but only " << faceCount << " faces exist in the Laguerre diagram.";
    throw LaguerreDiagram::ConstructionException(msg.str());
  }

  // for(auto& face : faces){
  //
  //   auto ec_start = face->ccb();
  //   auto ec = ec_start;
  //   do {
  //     if ( ec->has_source() )  std::cout << ec->source()->point();
  //     else  std::cout << "inf";
  //     std::cout << " -> ";
  //     if ( ec->has_target() )  std::cout << ec->target()->point();
  //     else  std::cout << "inf";
  //     std::cout  << std::endl;
  //
  //   } while ( ++ec != ec_start );
  //
  // }
  for(auto& face : faces){
    // std::cout << "Current point = " << face->dual()->point() << std::endl;

    laguerreCells.push_back( BoundOneCell(face) );
    AddInternalEdges(face);
  }
}

void LaguerreDiagram::AddInternalEdges(PowerDiagram::Face_handle const& face)
{

  // Get the point in the Regular triangulation at the center of this Laguerre cell
  auto vs = face->dual()->point();

  // Loop over the edges in the face
  Ccb_halfedge_circulator halfEdgeStart = face->ccb();
  Ccb_halfedge_circulator halfEdge = halfEdgeStart;

  Point_2 srcPt, tgtPt;

  // Loop over all of the edges in the face
  do {

    // Get the node in the triangulation on the other side of this edge
    auto vt = halfEdge->twin()->face()->dual()->point();

    // Compute the direction of dividing line between this cell and the other cell
    Vector_2 dir(vt.y()-vs.y(), -(vt.x()-vs.x()));


    // If the edge has a source, add the source point to the polygon
    if(halfEdge->has_source()){
      srcPt = halfEdge->source()->point();

      if(halfEdge->has_target()){
        tgtPt = halfEdge->target()->point();
      }else{
        tgtPt = halfEdge->source()->point() - infDist*dir;
      }

    // This edge starts at infinity, but has a target
    }else if(halfEdge->has_target()){
      srcPt = halfEdge->target()->point() + infDist*dir;
      tgtPt = halfEdge->target()->point();

    // This edge has neither a source nor a target
    }else{

      // Find the weighted mid point between the triangulation sites
      double ts = CGAL::to_double( CGAL::squared_distance(vt.point(), vs.point()).exact() );
      double h = (0.5/ts)*CGAL::to_double(ts + vs.weight() - vt.weight());
      Point_2 midPt = vs.point() + h*(vt.point()-vs.point());

      srcPt = midPt - infDist*dir;
      tgtPt = midPt + infDist*dir;
    }

    bool hasIntersection = bbox.ClipSegment(srcPt, tgtPt);

    // If the edge does not intersect the domain at all, just skip it
    if(hasIntersection){
      // This edge has both a source and target, so it is an internal edge.  Save it!
      unsigned int currInd = face->dual()->info();
      unsigned int otherInd = halfEdge->twin()->face()->dual()->info();
      internalEdges.at(currInd).push_back(std::make_tuple(otherInd, srcPt, tgtPt));
    }

  } while ( ++halfEdge != halfEdgeStart);

}


// bool LaguerreDiagram::ClipToBoundary(Point_2& srcPt, Point_2& tgtPt) const
// {
//
//   double sx = CGAL::to_double(srcPt.x());
//   double sy = CGAL::to_double(srcPt.y());
//   double tx = CGAL::to_double(tgtPt.x());
//   double ty = CGAL::to_double(tgtPt.y());
//
//   // Check if it's vertical
//   if(std::abs(sx-tx)<compTol){
//
//     // Check to see if it crosses the domain at all
//     if((sx>xBndLeft-compTol)&&(sx<xBndRight+compTol)&&(std::max(sy,ty)>yBndBottom-compTol)&&(std::min(sy,ty)<yBndTop+compTol)){
//       srcPt = Point_2(sx, std::min(yBndTop, std::max(yBndBottom, sy)));
//       tgtPt = Point_2(tx, std::min(yBndTop, std::max(yBndBottom, ty)));
//       return true;
//     }else{
//       return false;
//     }
//
//   // Check if it's horizontal
//   }else if(std::abs(sy-ty)<compTol){
//
//     if((sy>yBndBottom-compTol)&&(sy<yBndTop-compTol)&&(std::max(sx,tx)>xBndLeft-compTol)&&(std::min(sx,tx)<xBndRight+compTol)){
//       srcPt = Point_2(std::min(xBndRight, std::max(xBndLeft, sx)), sy);
//       tgtPt = Point_2(std::min(xBndRight, std::max(xBndLeft, tx)), ty);
//       return true;
//     }else{
//       return false;
//     }
//
//   // Otherwise it has a non zero finite slope
//   }else{
//
//
//     // Parameterize the line as srcPt + t*(tgtPt-srcPt) with t in [0,1]
//     Eigen::Matrix<double,4,1> qs(4), ps(4);
//
//     ps(0) = sx-tx;
//     ps(1) = tx-sx;
//     ps(2) = sy-ty;
//     ps(3) = ty-sy;
//
//     qs(0) = (sx-xBndLeft);///(tx-sx);
//     qs(1) = (xBndRight-sx);//(tx-sx);
//     qs(2) = (sy-yBndBottom);///(ty-sy);
//     qs(3) = (yBndTop-sy);///(ty-sy);
//
//
//     double u1 = 0.0;
//     for(unsigned int i=0; i<4; ++i){
//       if(ps(i)<0.0)
//         u1 = std::max(u1, qs(i)/ps(i));
//     }
//
//     double u2 = 1.0;
//     for(unsigned int i=0; i<4; ++i){
//       if(ps(i)>0.0)
//         u2 = std::min(u2, qs(i)/ps(i));
//     }
//
//     if(u1>u2){
//       return false;
//     }else{
//       srcPt = Point_2(sx + u1*ps(1), sy + u1*ps(3));
//       tgtPt = Point_2(sx + u2*ps(1), sy + u2*ps(3));
//       return true;
//     }
//   }
//
// }
