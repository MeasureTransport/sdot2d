#include "SDOT/Distances/GHKFunctions.h"

#include <math.h>

using namespace sdot;
using namespace sdot::distances;


double GHKFunctions::Evaluate(double z, double penaltyCoeff){
  return penaltyCoeff*(exp(z/penaltyCoeff)-1.0);
};

double GHKFunctions::Derivative(double z, double penaltyCoeff){
  return exp(z/penaltyCoeff);
};

double GHKFunctions::Derivative2(double z, double penaltyCoeff){
  return exp(z/penaltyCoeff)/penaltyCoeff;
}
