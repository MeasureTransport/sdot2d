#include "SDOT/LaguerreDiagram.h"

#include <gtest/gtest.h>

using namespace sdot;

TEST(LaguerreDiagram, Construction)
{
  // Create a Laguerre diagram with two points in a unit box
  Eigen::MatrixXd pts(2,2);
  pts << 0.25, 0.75,
         0.5, 0.5;

  Eigen::VectorXd costs = Eigen::VectorXd::Ones(2);

  LaguerreDiagram diag(0,1,0,1, pts, costs);

  Eigen::MatrixXd verts = diag.GetCellVertices(0);
  EXPECT_EQ(4, verts.cols());

  EXPECT_DOUBLE_EQ(0, verts(0,0));
  EXPECT_DOUBLE_EQ(0, verts(1,0));
  EXPECT_DOUBLE_EQ(0.5, verts(0,1));
  EXPECT_DOUBLE_EQ(0, verts(1,1));
  EXPECT_DOUBLE_EQ(0.5, verts(0,2));
  EXPECT_DOUBLE_EQ(1, verts(1,2));
  EXPECT_DOUBLE_EQ(0, verts(0,3));
  EXPECT_DOUBLE_EQ(1, verts(1,3));

  verts = diag.GetCellVertices(1);
  EXPECT_EQ(4, verts.cols());

  EXPECT_DOUBLE_EQ(0.5, verts(0,0));
  EXPECT_DOUBLE_EQ(0, verts(1,0));
  EXPECT_DOUBLE_EQ(1.0, verts(0,1));
  EXPECT_DOUBLE_EQ(0, verts(1,1));
  EXPECT_DOUBLE_EQ(1.0, verts(0,2));
  EXPECT_DOUBLE_EQ(1, verts(1,2));
  EXPECT_DOUBLE_EQ(0.5, verts(0,3));
  EXPECT_DOUBLE_EQ(1, verts(1,3));

  EXPECT_DOUBLE_EQ(0.5, diag.CellArea(0));
  EXPECT_DOUBLE_EQ(0.5, diag.CellArea(1));
}


TEST(LaguerreDiagram, GroundCost)
{
  // Create a Laguerre diagram with two points in a unit box
  Eigen::MatrixXd pts(2,2);
  pts << 0.25, 0.75,
         0.5, 0.5;

  Eigen::VectorXd costs(2);
  costs << 1.1, 1.0;

  LaguerreDiagram diag(0,1,0,1, pts, costs);

  Eigen::MatrixXd verts = diag.GetCellVertices(0);
  EXPECT_EQ(4, verts.cols());

  double midx = verts(0,1);
  EXPECT_DOUBLE_EQ(std::pow(midx-pts(0,0), 2.0) - costs(0),  std::pow(midx-pts(0,1), 2.0) - costs(1));
}
