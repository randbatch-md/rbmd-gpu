#pragma once
#include <memory>

#include "base_reader.h"
#include "common/types.h"

class MmapReader : public BaseReader {
 public:
  MmapReader(const std::string& filePath);
  virtual ~MmapReader();

  int Execute() override;

 protected:
  off_t _file_size;
  char* _mapped_memory;
  char* _line_start;
  rbmd::Id _locate;
};