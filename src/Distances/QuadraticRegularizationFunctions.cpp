#include "SDOT/Distances/QuadraticRegularizationFunctions.h"

using namespace sdot;
using namespace sdot::distances;


double QuadraticRegularizationFunctions::Evaluate(double z, double penaltyCoeff){
  if(z > -2*penaltyCoeff){
    return (0.25/penaltyCoeff)*z*z+ z;
  }else{
    return -penaltyCoeff;
  }
};

double QuadraticRegularizationFunctions::Derivative(double z, double penaltyCoeff){
  if(z > -2*penaltyCoeff){
    return (0.5/penaltyCoeff)*z + 1.0;
  }else{
    return 0;
  }
};

double QuadraticRegularizationFunctions::Derivative2(double z, double penaltyCoeff){
  if(z > -2*penaltyCoeff){
    return 0.5/penaltyCoeff;
  }else{
    return 0;
  }
};
