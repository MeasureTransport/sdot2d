#ifndef CONJUGATEFUNCTION_H
#define CONJUGATEFUNCTION_H

namespace sdot {

/** @defgroup ConjugateFunctions
    @brief The Fenchel-Legendre conjugate functions for unbalanced SDOT.
    @details In [Bourne et al., 2018], a family of unbalanced SDOT approaches
    is described.  The standard Kantorovich optimal transport problem, which
    requires that the marginals of the coupling are the source and target measure,
    is relaxed.  The hard marginal constraints are translated into penalties that
    encourage the marginals to match the source and target measures but do not
    require an exact match.   [Bourne et al., 2018] considers penalties that are
    based on discrepancies of the form
    \f[
    \mathcal{F}(\rho | \mu) = \left\{ \begin{array}{ll} \int_\Omega F\left(\frac{d\rho}{d\mu}\right)\end{array}\right.,
    \f]
    where \f$F\f$ is some convex lower semi-continuous function.   The function
    \f$F\f$ completely defines the discrepancy.  The dual formulation of the OT
    transport problem however, depends on the conjugate
    \f[
    F^\ast(z) = \sup_{s\geq 0}\left(z\dot s - F(s) \right).
    \]
    To solve the unbalance SDOT problem, we will need to evaluate \f$F^\ast\f$ as
    well as its first derivative.   We will also need to compute integral involving
    integrands of the form \f$F^|ast(w- c(x,x_i))\f$.  This group contains classes
    that support evaluating and integrating these conjugate functions.

    Each conjugate function class defines 8 static functions for evaluating the
    function and computing integrals over triangular and rectangular domains.
    The use of static functions is used, as opposed to inheritance and virtual
    functions, to make it easier to leverage OpenMP and Kokkos.
*/

} // namespace sdot

#endif
