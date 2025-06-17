#include "../include/mmap_reader.h"
#include <system_error>  // 用于错误处理

MmapReader::MmapReader(const std::string& filePath)
    : BaseReader(filePath),
      _file_size(0),
      _mapped_memory(nullptr),
      _line_start(nullptr),
      _locate(0) {}


int MmapReader::Execute() {
  try {
    // 使用 mio 进行内存映射（只读模式）
    mio::mmap_source mmap(_file_path);

    // 获取文件大小
    _file_size = mmap.size();

    // 获取映射的内存地址
    _mapped_memory = const_cast<char*>(mmap.data());  // 由于 mmap_source 是 const，需去除 const（谨慎操作）
    _line_start = _mapped_memory;

    // 将 mmap 对象移动到成员变量，确保生命周期延续
    _mmap = std::move(mmap);

    return 0;
  } catch (const std::system_error& e) {
    // 捕获并处理错误（如文件不存在、权限问题等）
    // log: e.what()
    return -1;
  }
}