#include "command_line.h"

#include <iostream>


CommandLine::CommandLine(int argc, char* argv[])
  : _opts("rbmd")
{
  //
  _opts.add_options()
      ("j,json", "JSON config file", cxxopts::value<std::string>())
      ("h,help", "Print help")
      ("v,version", "Print version");

  try {
    //_opts.add_options()("j", "json file", cxxopts::value<std::string>());
    _co = _opts.parse(argc, argv);
    ParseCommand();

    } catch (const std::exception&) {
    std::cerr << "Error parsing options: "  << "\n"
              << GetHelpText() << std::endl;
    throw std::runtime_error("Invalid command line");
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
  if (!_co.count("j") && !_co.count("h") && !_co.count("v")) {
    std::cerr << "Error: Missing required parameter -j\n"
              << GetHelpText() << std::endl;
    throw std::runtime_error("Missing required parameter");
  }
}
