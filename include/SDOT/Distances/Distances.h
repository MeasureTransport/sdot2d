#ifndef DISTANCES_DISTANCES_H
#define DISTANCES_DISTANCES_H

#include "SDOT/Distances/QuadratureDistance.h"
#include "SDOT/Distances/QuadraticRegularizationFunctions.h"
#include "SDOT/Distances/GHKFunctions.h"
#include "SDOT/Distances/Wasserstein2.h"

namespace sdot {
namespace distances{

/** @defgroup Distances
    @details In [(Bourne et al., 2018)](https://arxiv.org/pdf/1808.01962.pdf), a family of unbalanced SDOT approaches
    is described.  The standard Kantorovich optimal transport problem, which
    requires that the marginals of the coupling are the source and target measure,
    is relaxed.  The hard marginal constraints are translated into penalties that
    encourage the marginals to match the source and target measures but do not
    require an exact match.   [Bourne et al., 2018] considers penalties that are
    based on discrepancies of the form
    \f[
    \mathcal{F}(\rho | \mu) = \left\{ \begin{array}{ll} \int_\Omega F\left(\frac{d\rho}{d\mu}\right) & \rho \ll \mu \\ \infty & \text{otherwise} \end{array}\right.,
    \f]
    where \f$F\f$ is some convex lower semi-continuous function and \f$\rho \ll \mu\f$
    indicates that \f$\rho\f$ is absolutely continuous with respect to \f$\mu\f$.
    (i.e., the support of \f$\rho\f$ is contained in the support of \f$\mu\f$.)
    See Definition 2.2 of [(Bourne et al., 2018)](https://arxiv.org/pdf/1808.01962.pdf)
    for details.   The function
    \f$F\f$ completely defines the discrepancy.  The dual formulation of the OT
    transport problem however, depends on the conjugate function
    \f[
    F^\ast(z) = \sup_{s\geq 0}\left(z\cdot s - F(s) \right).
    \f]
    To solve the unbalance SDOT problem, we will need to evaluate \f$F^\ast\f$,
    as well as its first derivative. We will also need to compute integrals involving
    integrands of the form \f$F^|ast(w_i- c(x,x_i))\f$.  This group contains classes
    that support evaluating and integrating these conjugate functions.

    Each "distance" class defines the conjugate function for a specific unbalanced
    OT "distance".  It does this through static functions for evaluating the
    function and computing integrals over triangular and rectangular domains.
    The use of static functions is used, as opposed to inheritance and virtual
    functions, to make it easier to leverage OpenMP and Kokkos.

    In all ddistances, the cost is assumed to be
    \f[
    c(x,y) = \frac{1}{2}\|x-y\|^2.
    \f]

    Each class in this group defines the following static functions, which are
    called by the templated SemiDiscreteOT class.
    - `Evaluate`
      Evaluates the conjugate function \f$F^\ast(z)\f$ at a point \f$z\f$.
    - `Derivative`
      Evaluates the derivative of the conjugate function, i.e., \f$(F^\ast)^\prime(z) = \frac{\partial}{\partial z}F^\ast \f$.
    - `TriangularIntegral`
      Returns the integral of \f$F^\ast(w_i-c(x,x_i))\f$ wrt \f$x\f$ over a triangular region.  The cost \f$w_i\in \mathbb{R}\f$ and point \f$x_i\in\mathbb{R}^2\f$ are fixed.
    - `RectangularIntegral`
      Returns the integral of \f$F^\ast(w_i-c(x,x_i))\f$ wrt \f$x\f$ over a rectangular region.  Like the `TriangularIntegral` function, the cost \f$w_i\in \mathbb{R}\f$ and point \f$x_i\in\mathbb{R}^2\f$ are fixed.
    - `TriangularIntegralDeriv`
      Returns the integral of \f$(F^\ast)^\prime(w_i-c(x,x_i))\f$ wrt \f$x\f$ over a triangular region.
    - `RectangularIntegralDeriv`
      Returns the integral of \f$(F^\ast)^\prime(w_i-c(x,x_i))\f$ wrt \f$x\f$ over a line segment.
    - `LineIntegralDeriv`
      Returns the integral of \f$(F^\ast)^\prime(w_i-c(x,x_i))\f$ wrt \f$x\f$ over a rectangular region
    - `TriangularIntegralDeriv2`
      Returns the integral of \f$(F^\ast)^{\prime\prime}(w_i-c(x,x_i))\f$ wrt \f$x\f$ over a triangular region.
    - `RectangularIntegralDeriv2`
      Returns the integral of \f$(F^\ast)^{\prime\prime}(w_i-c(x,x_i))\f$ wrt \f$x\f$ over a rectangular region.

*/

  typedef QuadratureDistance<QuadraticRegularizationFunctions> QuadraticRegularization;
  typedef QuadratureDistance<GHKFunctions> GHK;

} // namespace distances
} // namespace sdot

#endif
