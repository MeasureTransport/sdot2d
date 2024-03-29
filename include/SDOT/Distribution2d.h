#ifndef DISTRIBUTION2D_H_
#define DISTRIBUTION2D_H_

#include "SDOT/RegularGrid.h"

#include <memory>

namespace sdot{


  /** Abstract base class for defining distributions on regular grids.
  */
  class Distribution2d {

  public:
    Distribution2d(std::shared_ptr<RegularGrid> const& gridIn) : grid(gridIn){};

    virtual ~Distribution2d() = default;

    /** Returns the constant density over a particular grid cell. */
    virtual double Density(unsigned int xInd, unsigned int yInd) const = 0;

    virtual std::shared_ptr<RegularGrid> const& Grid() const{return grid;};

    virtual double TotalMass() const = 0;
    
  protected:
    std::shared_ptr<RegularGrid> grid;

  }; // class Distribution2d
}

#endif // #ifndef DISTRIBUTION2D_H_
