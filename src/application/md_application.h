#pragma once
#include "application.h"
#include "../common/thermo_tui_view.hpp"

class MDApplication : public Application {
 public:
  MDApplication(int argc, char* argv[]);
  ~MDApplication() = default;

 public:
  int Execute() override;
  void AddSimulate();

 private:
  int ReadMDData();
  std::shared_ptr<CommandLine> _cmd; //
  std::shared_ptr<ThermoStatsTUIView> _thermo_view;
};
