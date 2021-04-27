#include "SDOT/DiscretizedDistribution.h"
#include "SDOT/Assert.h"

#include <iostream>
using namespace sdot;

DiscretizedDistribution::DiscretizedDistribution(std::shared_ptr<RegularGrid> const& gridIn,
                                                 Eigen::MatrixXd              const& densIn) : Distribution2d(gridIn),
                                                                                               densVals(densIn)
{
  SDOT_ASSERT(densVals.rows()==grid->NumCells(0));
  SDOT_ASSERT(densVals.cols()==grid->NumCells(1));
};

Eigen::MatrixXd DiscretizedDistribution::Normalize(std::shared_ptr<RegularGrid> const& grid,
                                                   Eigen::MatrixXd const& vals)
{
  double sum = (grid->dx*grid->dy*vals).sum();
  return vals / sum;
}

double DiscretizedDistribution::TotalMass() const
{
  return (grid->dx*grid->dy*densVals).sum();
}


double DiscretizedDistribution::Density(unsigned int xInd, unsigned int yInd) const
{
  return densVals(xInd,yInd);
}
