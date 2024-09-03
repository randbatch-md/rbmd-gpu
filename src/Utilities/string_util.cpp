#include "string_util.h"

namespace rbmd {
bool IsLegalLine(const std::string& line) {
  return (!line.empty() && line != "\n" && line != "\r" && line != "\n\r" &&
          line != "\r\n");
}
}  // namespace rbmd
