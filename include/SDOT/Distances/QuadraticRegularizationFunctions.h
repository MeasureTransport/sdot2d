#ifndef DISTANCES_QUADRATICREGULARIZATIONFUNCTIONS_H
#define DISTANCES_QUADRATICREGULARIZATIONFUNCTIONS_H


namespace sdot {
namespace distances{

  /** @class QuadraticRegularizationFunctions
      @ingroup Distances
      @details The Quadratic Regularization described in Example 2.8 of
      [Bourne et al., 2018](https://arxiv.org/pdf/1808.01962.pdf) uses the function
      \f[
      F(s) = (s-1)^2 ,
      \f]
      which results in a conjugate function \f$F^\ast(z)\f$ of the form
      \f[
      F^\ast(z) = \left\{ \begin{array}{ll} \frac{z^2}{4} + z & \text{if} z\geq -2\\ -1 & \text{otherwise}\endd{array}\right.
      \f]
      This class provides implementations of \f$F^\ast(z)\f$ that can be used with
      the QuadratureDistance template class to provide all the integrals necesssary
      for unbalanced OT with the quadratic penalty.
  */
  class QuadraticRegularizationFunctions{

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
