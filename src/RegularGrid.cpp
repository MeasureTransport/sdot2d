#include "SDOT/RegularGrid.h"

#include <iostream>
#include <cassert>
#include <cmath>

using namespace sdot;

RegularGrid::RegularGrid(double       xMin_, double       yMin_,
                         double       xMax_, double       yMax_,
                         unsigned int Nx_,   unsigned int Ny_) : xMin(xMin_), yMin(yMin_),
                                                                 xMax(xMax_), yMax(yMax_),
                                                                 Nx(Nx_),     Ny(Ny_),
                                                                 dx((xMax_-xMin_)/double(Nx_)), dy((yMax_-yMin_)/double(Ny_))
{
  assert(xMin<xMax);
  assert(yMin<yMax);
};

unsigned int RegularGrid::NumCells(int dim) const
{
  if(dim==0){
    return Nx;

  }else if(dim==1){
    return Ny;

  }else{
    return Nx*Ny;
  }
}

unsigned int RegularGrid::NumNodes(int dim) const
{
  if(dim==0){
    return Nx+1;

  }else if(dim==1){
    return Ny+1;

  }else{
    return (Nx+1)*(Ny+1);
  }
}

double RegularGrid::LeftSide(double x) const
{
  return xMin + dx*LeftNode(x);
}

/** Returns the x coordinate of the left side of the cells containig x. */
double RegularGrid::RightSide(double x) const
{
  return xMin + dx*RightNode(x);
}

unsigned int RegularGrid::RightNode(double x) const
{
  double indDouble = (x-xMin)/dx;
  double distToNode = std::abs(std::round(indDouble) - indDouble);

  if(distToNode<compTol){
    return std::round(indDouble);
  }else{
    return std::ceil((x-xMin)/dx);
  }
}

unsigned int RegularGrid::TopNode(double y) const
{
  double indDouble = (y-yMin)/dy;
  double distToNode = std::abs(std::round(indDouble) - indDouble);
  if(distToNode<compTol){
    return std::round(indDouble);
  }else{
    return std::ceil((y-xMin)/dy);
  }
}

unsigned int RegularGrid::LeftNode(double x) const
{
  double indDouble = (x-xMin)/dx;
  double distToNode = std::abs(std::round(indDouble) - indDouble);

  if(distToNode<compTol){
    return std::round(indDouble);
  }else{
    return std::floor((x-xMin)/dx);
  }
}


unsigned int RegularGrid::BottomNode(double y) const
{
  double indDouble = (y-yMin)/dy;
  double distToNode = std::abs(std::round(indDouble) - indDouble);

  if(distToNode<compTol){
    return std::round(indDouble);
  }else{
    return std::floor((y-xMin)/dy);
  }
}
