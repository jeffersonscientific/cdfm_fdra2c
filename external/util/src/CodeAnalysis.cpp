#include "util/include/CodeAnalysis.hpp"

namespace util {
  namespace Ca {
    
    util::Timer* GetTimer() {
      static Timer t;
      return &t;
    }

  }
}
