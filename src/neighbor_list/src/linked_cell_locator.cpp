#include "linked_cell_locator.h"
LinkedCellLocator& LinkedCellLocator::GetInstance() {
  static LinkedCellLocator instance;
  return instance;
}
std::shared_ptr<LinkedCell> LinkedCellLocator::GetLinkedCell() {
  if (nullptr == _linked_cell) {
    this->_linked_cell = std::make_shared<LinkedCell>();
  }
  return this->_linked_cell;
}