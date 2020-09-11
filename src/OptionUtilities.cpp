#include "SDOT/OptionUtilities.h"

using namespace sdot;


double sdot::GetOpt(std::string const& key,
                    OptionList  const& opts,
                    double             def)
{
  auto it = opts.find(key);
  if(it == opts.end()){
    return def;
  }else{
    return it->second;
  }
}
