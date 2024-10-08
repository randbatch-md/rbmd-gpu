#include "application/md_application.h"
#include <hip/hip_runtime.h>
int main(int argc, char* argv[]) {
  CHECK_RUNTIME(hipSetDevice(0));
  std::shared_ptr<Application> app = std::make_shared<MDApplication>(argc, argv);
  app->Run();

  return 0;
}
