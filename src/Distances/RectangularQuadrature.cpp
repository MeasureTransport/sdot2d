#include "SDOT/Distances/RectangularQuadrature.h"
#include "SDOT/Distances/LineQuadrature.h"

#include <Eigen/Core>
#include "SDOT/Assert.h"

using namespace sdot::distances;


std::pair<Eigen::Matrix2Xd, Eigen::VectorXd> RectangularQuadrature::Get(unsigned int degree)
{
  // Get the 1d scalar points and weights
  Eigen::VectorXd pts1d, wts1d;
  std::tie(pts1d, wts1d) = LineQuadrature::Get(degree);

  const int N1 = pts1d.size();
  int ind;

  Eigen::Matrix2Xd pts(2,N1*N1);
  Eigen::VectorXd wts(N1*N1);
  for(int xind=0; xind<N1; ++xind){
    for(int yind=0; yind<N1; ++yind){
      ind = yind + xind*N1;
      pts(0,ind) = pts1d(xind);
      pts(1,ind) = pts1d(yind);
      wts(ind) = wts1d(xind)*wts1d(yind);
    }
  }

  return std::make_pair(pts,wts);
}
