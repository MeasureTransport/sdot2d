#ifndef DISTANCES_GHKFUNCTIONS_H
#define DISTANCES_GHKFUNCTIONS_H


namespace sdot {
namespace distances{

  /** @class GHKFunctions
      @ingroup Distances
      @brief Defines conjugate function for KL divergence marginal penalty.
      @details  Using the terminology of [Bourne et al.], the Gaussian Hellinger-Kantorovich
      (or GHK) distance arises when the KL divergence is used as a penalty function in the
      unbalanced OT setting.  This class provides implementations of the conjugate function
      \f$F^\ast(z)\f$ corresponding to the GHK distance.  It can be used with
      the QuadratureDistance template class to provide all the integrals necesssary
      for unbalanced OT with the quadratic penalty.
  */
  class GHKFunctions{

  public:
    /** Evaluates the conjugate function \f$F^\ast(z)\f$ at a point \f$z\f$. */
    static double Evaluate(double z, double penaltyCoeff=1);

    /** Evaluates the derivative of the conjugate function \f$\frac{\partial}{\partial z}F^\ast(z)\f$ at a point \f$z\f$. */
    static double Derivative(double z, double penaltyCoeff=1);

    /** Evaluates the second derivative of the conjugate function \f$\frac{\partial}{\partial z}F^\ast(z)\f$ at a point \f$z\f$. */
    static double Derivative2(double z, double penaltyCoeff=1);
  };

} // namespace distances
} // namespace sdot

#endif
