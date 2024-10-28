#include "linked_cell/linked_cell_locator.h"

#include "data_manager.h"
#include "model/md_data.h"
LinkedCellLocator& LinkedCellLocator::GetInstance() {
  static LinkedCellLocator instance;
  return instance;
}
std::shared_ptr<LinkedCell> LinkedCellLocator::GetLinkedCell() {
  if (nullptr == _linked_cell) {
    this->_linked_cell = std::make_shared<LinkedCell>();
    auto& box = DataManager::getInstance().getMDData()->_h_box;
    auto short_edge =
        MIN(MIN(box->_length[0], box->_length[1]), box->_length[2]);
    if (this->_linked_cell->_cutoff >= short_edge) {
      std::cout << "\033[31mError: cutoff must be less than the shortest side "
                   "of the box.\033[0m"
                << std::endl;
      exit(0);
    }
    // RBL
    if (DataManager::getInstance().getConfigData()->Get<std::string>(
            "type", "hyper_parameters", "neighbor") == "RBL") {
      const auto _r_core =
          DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
              "r_core", "hyper_parameters", "neighbor");
      if (_r_core >= _linked_cell->_cutoff) {
        std::cout << "\033[31mError r_core must be less than  the "
                     "cutoff.\033[0m"
                  << std::endl;
        exit(0);
      }
      _linked_cell->_cell_count_within_cutoff = static_cast<rbmd::Id>(
          std::ceil(static_cast<double>(_linked_cell->_cutoff / _r_core)));
    }
    _linked_cell->Build(DataManager::getInstance().getMDData()->_h_box.get());
    _linked_cell->InitializeCells();
  }
  return this->_linked_cell;
}