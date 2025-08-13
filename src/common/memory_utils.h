#pragma once
#include <cstddef>
#include <string>

// 在您的日志记录系统中包含适当的头文件，这里以 spdlog 为例
#include <spdlog/spdlog.h>

#if defined(__linux__)
#include <fstream>
#include <string>
#include <unistd.h>
#elif defined(_WIN32)
#include <windows.h>
#endif

namespace rbmd {
namespace utils {

// 获取当前系统可用物理内存（以字节为单位）
inline size_t get_available_memory() {
#if defined(__linux__)
  std::ifstream meminfo("/proc/meminfo");
  std::string line;
  while (std::getline(meminfo, line)) {
    if (line.rfind("MemAvailable:", 0) == 0) {
      try {
        size_t mem_kb = std::stoull(line.substr(line.find_first_of("0123456789")));
        return mem_kb * 1024;
      } catch (...) {
        break;
      }
    }
  }
  long pages = sysconf(_SC_AVPHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  if (pages != -1 && page_size != -1) {
    return static_cast<size_t>(pages) * static_cast<size_t>(page_size);
  }
#elif defined(_WIN32)
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if (GlobalMemoryStatusEx(&status)) {
    return status.ullAvailPhys;
  }
#endif
  spdlog::warn("Cannot determine available system memory. Using a conservative default of 2GB.");
  return 2ULL * 1024 * 1024 * 1024; // 返回一个保守的默认值，2GB
}

} // namespace utils
} // namespace rbmd