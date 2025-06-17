#pragma once
#include <memory>
#include "common/mio.hpp" // 新增：引入 mio 头文件

#include "base_reader.h"
#include "common/types.h"

class MmapReader : public BaseReader {
 public:
  MmapReader(const std::string& filePath);
  virtual ~MmapReader() = default;  // 修改：不再需要手动释放资源

  int Execute() override;

 protected:
  off_t _file_size;
  char* _mapped_memory;    // 保持兼容性（外部可能依赖此指针）
  char* _line_start;
  rbmd::Id _locate;

 private:
  mio::mmap_source _mmap;  // 新增：mio 的内存映射对象（RAII 管理）
};