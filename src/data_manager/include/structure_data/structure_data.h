#pragma once
#include "../common/object.h"
#include "../common/types.h"
#include "../common/rbmd_define.h"

// 内存泄漏
class StructureData : public Object {
public:
  virtual bool checkStructure() const = 0;

  ~StructureData() override;

  /// position on host
  rbmd::Real* _h_px = nullptr;
  rbmd::Real* _h_py = nullptr;
  rbmd::Real* _h_pz = nullptr;

  /// atoms is on host
  rbmd::Id* _h_atoms_id = nullptr;

  /// atoms type on host
  rbmd::Id* _h_atoms_type = nullptr;

  /// molecular id on host
  rbmd::Id* _h_molecular_id = nullptr;

  // atoms flag on host
  rbmd::Id* _h_flagX = nullptr;
  rbmd::Id* _h_flagY = nullptr;
  rbmd::Id* _h_flagZ = nullptr;

  /// velocity on host
  rbmd::Real* _h_vx = nullptr;
  rbmd::Real* _h_vy = nullptr;
  rbmd::Real* _h_vz = nullptr;

  /// force on host
  rbmd::Real* _h_fx = nullptr;
  rbmd::Real* _h_fy = nullptr;
  rbmd::Real* _h_fz = nullptr;

  //
  rbmd::Real* _h_evdwl = nullptr;
};

inline StructureData::~StructureData() {
  FREE(_h_px);
  FREE(_h_py);
  FREE(_h_pz);

  FREE(_h_atoms_id);
  FREE(_h_atoms_type);
  FREE(_h_molecular_id);

  FREE(_h_flagX);
  FREE(_h_flagY);
  FREE(_h_flagZ);

  FREE(_h_vx);
  FREE(_h_vy);
  FREE(_h_vz);

  FREE(_h_fx);
  FREE(_h_fy);
  FREE(_h_fz);

  FREE(_h_evdwl);
}