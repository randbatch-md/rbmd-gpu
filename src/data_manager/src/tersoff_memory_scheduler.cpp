#include "./include/scheduler/tersoff_memory_scheduler.h"

bool TersoffMemoryScheduler::asyncMemoryH2D() {
  if (false == MemoryScheduler::asyncMemoryH2D()) {
    // log
    return false;
  }

  return true;
}

bool TersoffMemoryScheduler::asyncMemoryD2H() { return true; }
