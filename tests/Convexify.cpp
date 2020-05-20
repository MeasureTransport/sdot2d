#include <Eigen/Core>

#include <vector>
#include <memory>

#include "SDOT/Polygon.h"

using namespace sdot;

void MakePolygonConvex(){

  // First, construct a non-convex polygon
  /*
      v4     v2
      | \   /|
      |  \ / |
      |  v3  |
      |      |
      v0-----v1
 */
  Eigen::Matrix2Xd pts(2,5);
  pts << 0, 1, 1, 0.5, 0,
         0, 0, 1, 0.5, 1;

  Polygon poly(pts);

  // Test the IsConvex() function.
  std::cout << "Is non-convex base polygon being detected? ";
  if(poly.IsConvex()){
    std::cout << " NO!"<< std::endl;
  }else{
    std::cout << " YES!" << std::endl;
  }


  // Test the ability to construct a convex partition
  std::vector<std::shared_ptr<Polygon>> polys = poly.ConvexPartition();

  std::cout << "Convex Polygons:" << std::endl;
  for(int i=0; i<polys.size(); ++i){
    std::cout << "  Poly"<< i+1 << " IsConvex()? " << polys.at(i)->IsConvex() << std::endl;
  }
}


int main(){

  MakePolygonConvex();
  return 0;
}
