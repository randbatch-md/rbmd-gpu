#include <rocm-core/rocm_version.h>
#include <cstdio>
#ifndef ROCM_VERSION_PATCH
#define ROCM_VERSION_PATCH 0
#endif
#define STRINGIFYHELPER(x) #x
#define STRINGIFY(x) STRINGIFYHELPER(x)
int main() {
  printf("%d.%d.%s", ROCM_VERSION_MAJOR, ROCM_VERSION_MINOR, STRINGIFY(ROCM_VERSION_PATCH));
  return 0;
}
