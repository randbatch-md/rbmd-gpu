#pragma once

#include <memory>

#include "linked_cell.h"

class LinkedCellLocator {
public:
  LinkedCellLocator(const LinkedCellLocator &) = delete;

  LinkedCellLocator &operator=(const LinkedCellLocator &) = delete;
private:
  LinkedCellLocator() {
  }

  ~LinkedCellLocator() = default;


  std::shared_ptr<LinkedCell> _linked_cell = nullptr;
public:
  static LinkedCellLocator &GetInstance();

  std::shared_ptr<LinkedCell> GetLinkedCell();

};

