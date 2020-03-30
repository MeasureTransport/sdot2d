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
        @param[in] probsIn A matrix of the constant density for each grid cell.  The rows of
                   the matrix correspond to the x index.  This is a density not a
                   probability measure.  Thus, densIn(xInd,yInd)*dx*dy is the probability
                   assigned to cell xInd,yInd.
    */
    DiscretizedDistribution(std::shared_ptr<RegularGrid> const& gridIn,
                            Eigen::MatrixXd              const& densIn);

    virtual double Density(unsigned int xInd, unsigned int yInd) const override;

  protected:

    /** Normalizes density values so that the density integrates to 1 over the
        grid domain.
        @param[in] grid Grid defining the domain and grid sizing.
        @param[in] vals Unnormalized density values that are assumed constant
                        over each grid cell.
    */
    static Eigen::MatrixXd Normalize(std::shared_ptr<RegularGrid> const& grid,
                                     Eigen::MatrixXd              const& vals);

    const Eigen::MatrixXd densVals;

  }; // class DiscretizedDistribution
}

#endif // #ifndef DISCRETIZEDDISTRIBUTION_H_
