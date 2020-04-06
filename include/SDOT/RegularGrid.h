#ifndef REGULARGRID_H_
#define REGULARGRID_H_

#include "SDOT/BoundingBox.h"

namespace sdot {

/** Class for defining a regular grid in 2d. */
class RegularGrid
{
public:
  RegularGrid(double       xMin_, double          yMin_,
              double       xMax_, double          yMax_,
              unsigned int Nx_,   unsigned int    Ny_);

  RegularGrid(BoundingBox const& bbox,
              unsigned int Nx_,
              unsigned int Ny_);

  /** Returns the number of cells in each direction (dim=0,1) or the total number
      of cells in the grid (dim=-1).  Defaults to total number in grid.
  */
  unsigned int NumCells(int dim=-1) const;

  /** Returns the number of nodes in each direction (dim=0,1) or the total number
      of nodes in the grid (dim=-1).  Defaults to total number in grid.
  */
  unsigned int NumNodes(int dim=-1) const;

  /** Returns the x coordinate of the left side of the cells containig x. */
  double LeftSide(double x) const;

  /** Returns the x coordinate of the left side of the cells containig x. */
  double RightSide(double x) const;

  /** Returns the y coordinate of the bottom side of the cells containig y. */
  double BottomSide(double y) const;

  /** Returns the y coordinate of the top side of the cells containig y. */
  double TopSide(double y) const;

  /** Returns the index of the node to the left side of the cell containing x.
      If x is within compTol of a node location, that location is returned.
  */
  unsigned int LeftNode(double x) const;

  /** Returns the index of the node to the right side of the cell containing x.
      If x is within compTol of a node location, that location is returned.
  */
  unsigned int RightNode(double x) const;

  /** Returns the index of the node on the top of the the cell containing y.
      If y is within compTol of a node location, that location is returned.
  */
  unsigned int TopNode(double y) const;

  /** Returns the index of the node on the bottom of the the cell containing y.
      If y is within compTol of a node location, that location is returned.
  */
  unsigned int BottomNode(double y) const;

  /** Returns the center of a grid cell.
      @param[in] xInd The x index of the grid cell.
      @param[in] yInd The y index of the grid cell.
  */
  Eigen::Vector2d Center(unsigned int xInd, unsigned int yInd) const;

  /// Bounding box
  const double xMin, yMin, xMax, yMax;

  /// Number of cells in x and y directions
  const unsigned int Nx, Ny;

  const double dx, dy;

  const double compTol = 1e-13;
};

} // namespace sdot

#endif
