#include "data_manager.h"
#include "linked_cell/linked_cell_locator.h"
#include "model/md_data.h"
LinkedCellLocator& LinkedCellLocator::GetInstance() {
  static LinkedCellLocator instance;
  return instance;
}
std::shared_ptr<LinkedCell> LinkedCellLocator::GetLinkedCell() {
  if (nullptr == _linked_cell) {
    this->_linked_cell = std::make_shared<LinkedCell>();
    _linked_cell->Rebuild(DataManager::getInstance().getMDData()->_h_box.get());
    _linked_cell->InitializeCells();
  }
  return this->_linked_cell;
}