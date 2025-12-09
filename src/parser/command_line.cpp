#include "command_line.h"
#include "output/include/Logger.hpp"
#include <iostream>

CommandLine::CommandLine(int argc, char* argv[])
  : _opts("rbmd")
{
  //
  _opts.add_options()
      ("j,json", "JSON config file", cxxopts::value<std::string>())
      ("t,tui", "Enable terminal UI", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
      ("h,help", "Print help")
      ("v,version", "Print version");

  try {
    //_opts.add_options()("j", "json file", cxxopts::value<std::string>());
    _co = _opts.parse(argc, argv);

    for (int i = 1; i < argc; ++i) {
      if (strcmp(argv[i], "-j") == 0 || strcmp(argv[i], "--json") == 0) {
        if (i + 1 >= argc || argv[i+1][0] == '-') {
          Logger::Instance().error("\033[31m Missing value for -j parameter {}\033[0m", GetHelpText() );
        }
      }
    }

    ParseCommand();

  } catch (const std::exception& e) {
  }
}

void CommandLine::ParseCommand() {
  if (_co.count("h")) {
    std::cout << GetHelpText();
    exit(0); //
  }

  if (_co.count("v")) {
    std::cout << GetVersionText();
    exit(0);
  }
  //
  if (!_co.count("j")) {
    Logger::Instance().error("\033[31m Missing required parameter -j: {}\033[0m", GetHelpText() );
  }
}
