#ifndef OPTIONUTILITIES_H
#define OPTIONUTILITIES_H

#include <unordered_map>
#include <string>
#include <limits>

namespace sdot{

  typedef std::unordered_map<std::string, double> OptionList;

  /**
  This function is primarily used inside the SemiDiscreteOT class to extract
  optimization options from an unordered_map mapping strings (option names) to
  doubles.   If the option name (key) does not exist in the unordered_map, the
  default value is returned.  If the default is not specified, a NaN is returned.

  @param[in]  key The name of the option to extract.
  @param[in] opts The unordered_map containing all of the possible options.
  @param[in] default Optional specification of default option.

  @return The value of the option.
  */
  double GetOpt(std::string const& key,
                OptionList  const& opts,
                double             def = std::numeric_limits<double>::quiet_NaN());

} // namespace sdot

#endif
