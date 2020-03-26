#ifndef DISCRETIZEDDISTRIBUTION_H_
#define DISCRETIZEDDISTRIBUTION_H_

#include "SDOT/Distribution2d.h"
#include <Eigen/Core>

namespace sdot{

  /**
  Defines a probability measure on a uniform cartesian grid, where the probability
  assigned to each grid cell is known or has already been computed.
  */
  class DiscretizedDistribution : public Distribution2d {
  public:

    /** Construct the discretized distribution.
        @param[in] gridIn The uniform regular grid used to define this distribution.
        @param[in] probsIn A matrix of probabilities for each grid cell.  The rows of
                   the matrix correspond to the x index.  i.e., probsIn(xInd,yInd)
                   is the probability assigned to cell xInd, yInd.
    */
    DiscretizedDistribution(std::shared_ptr<RegularGrid> const& gridIn,
                            Eigen::MatrixXd              const& probsIn);

    virtual double Probability(unsigned int xInd, unsigned int yInd) const override;

  protected:
    const Eigen::MatrixXd probs;

  }; // class DiscretizedDistribution
}

#endif // #ifndef DISCRETIZEDDISTRIBUTION_H_
