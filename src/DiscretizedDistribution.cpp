#include "SDOT/DiscretizedDistribution.h"

using namespace sdot;

DiscretizedDistribution::DiscretizedDistribution(std::shared_ptr<RegularGrid> const& gridIn,
                                                 Eigen::MatrixXd              const& probsIn) : Distribution2d(gridIn),
                                                                                                probs(probsIn)
{
  assert(probs.rows()==grid->NumCells(0));
  assert(probs.cols()==grid->NumCells(1));
};

double DiscretizedDistribution::Probability(unsigned int xInd, unsigned int yInd) const
{
  return probs(xInd,yInd);
}
