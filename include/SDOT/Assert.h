#ifndef SDOT_ASSERT_H
#define SDOT_ASSERT_H

#include <cstdlib>
#include <cstdio>



#ifndef SKIP_SDOT_ASSERT

  #define SDOT_STR(x) #x
  #define SDOT_ASSERT(x) if (!(x)) { printf("\nAssertion failed: (%s)\n  function: \"%s\"\n  file: \"%s\"\n  line: %d\n\n", SDOT_STR(x), __PRETTY_FUNCTION__, __FILE__, __LINE__); abort(); }

#else

  #define SDOT_ASSERT(cond) \
      do { (void)sizeof(cond); } while(0)

#endif // #ifndef SKIP_SDOT_ASSERT

#endif // #ifndef SDOT_ASSERT_H
