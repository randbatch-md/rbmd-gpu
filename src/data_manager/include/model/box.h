#pragma once
#include "../../../common/rbmd_define.h"
#include "../../../common/types.h"

class Box {
 public:
  // 定义枚举类型 BoxType
  enum class BoxType {
    ORTHOGONAL,  // 正交
    TRICLINIC    // 三斜
  };

  BoxType _type = BoxType::ORTHOGONAL;  // 默认正交盒子
  Box() = default;

  /**
   * 读取配置文件时初始化盒子
   * @param box_type  盒子类型，枚举变量
   * @param coord_min 盒子右下角坐标（xyz最小值）
   * @param coord_max 盒子左上角坐标（xyz最达值）
   * @param pbc 是否使用周期性边界条件 bool数组，对应xyz维度
   */
  void Setup(BoxType box_type, const rbmd::Real coord_min[3],
            const rbmd::Real coord_max[3], const bool pbc[3]);

  // TODO 边长角度  不经常访问的就直接做成对齐数组了    xyz 半的xyz
  /// 在x维度使用周期性边界条件
  bool _pbc_x = true;
  /// 在y维度使用周期性边界条件
  bool _pbc_y = true;
  /// 在z维度使用周期性边界条件
  bool _pbc_z = true;

  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 6)) _length[6]{};
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 6)) _length_inv[6]{};

  /// （local）盒子的左下角坐标（x,y,z）
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) _coord_min[3]{};
  /// （local）盒子的右上角坐标（x,y,z）
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) _coord_max[3]{};

  /// 以单元格为单位的box的宽度（x,y,z）
  rbmd::Id ALIGN(ALIGN_SIZE(rbmd::Id, 3)) _box_width_as_cell_units[3]{};
};

static __host__ __device__ __forceinline__ void MinImageDistance( Box box,
                                                          rbmd::Real& dx,
                                                          rbmd::Real& dy,
                                                          rbmd::Real& dz) {
  if (box._type == Box::BoxType::ORTHOGONAL) {
    if (box._pbc_x) {
      if (ABS(dx) > box._length[0] * 0.5) {
        dx -= (dx > 0 ? box._length[0] : -box._length[0]);
      }
    }
    if (box._pbc_y) {
      if (ABS(dy) > box._length[1] * 0.5) {
        dy -= (dy > 0 ? box._length[1] : -box._length[1]);
      }
    }
    if (box._pbc_z) {
      if (ABS(dz) > box._length[2] * 0.5) {
        dz -= (dz > 0 ? box._length[2] : -box._length[2]);
      }
    }
  }
  // TODO else: tri
}

static __host__ __device__ __forceinline__ void MinImageDistance_fix(
     Box box, rbmd::Real& dx, rbmd::Real& dy, rbmd::Real& dz) {
  if (box._type == Box::BoxType::ORTHOGONAL) {
    if (box._pbc_x) {
      dx -= box._length[0] * RINT(dx * box._length_inv[0]);
    }
    if (box._pbc_y) {
      dy -= box._length[1] * RINT(dy * box._length_inv[1]);
    }
    if (box._pbc_z) {
      dz -= box._length[2] * RINT(dz * box._length_inv[2]);
    }
  }
  // TODO else: tri
}

static __host__ __device__ __forceinline__ Real3 MinImageDistanceVec(
  const rbmd::Real& px_1,
  const rbmd::Real& py_1,
  const rbmd::Real& pz_1,
  const rbmd::Real& px_2,
  const rbmd::Real& py_2,
  const rbmd::Real& pz_2,
   Box box)
{
  Real3 vec;
  vec.x = px_1 - px_2;
  vec.y = py_1 - py_2;
  vec.z = pz_1 - pz_2;

  if (box._type == Box::BoxType::ORTHOGONAL) {
     // X
      if (box._pbc_x)
      {
        if (ABS(vec.x) > box._length[0] * 0.5)
        {
          vec.x -= (vec.x > 0 ? box._length[0] : -(box._length[0]));
        }
      }
      // Y
      if (box._pbc_y)
      {
        if (ABS(vec.y) > box._length[1] * 0.5)
        {
          vec.y -= (vec.y > 0 ? box._length[1] : -(box._length[1]));
        }
      }

      // Z
      if (box._pbc_z)
      {
        if (ABS(vec.z) > box._length[2] * 0.5)
        {
          vec.z -= (vec.z > 0 ? box._length[2] : -(box._length[2]));
        }
      }
  }
  return vec;
}

__host__ __device__ __forceinline__ void ApplyPBC(
     Box box, rbmd::Real& px, rbmd::Real& py, rbmd::Real& pz,
    rbmd::Id& flag_px_tid, rbmd::Id& flag_py_tid, rbmd::Id& flag_pz_tid) {
  if (box._type == Box::BoxType::ORTHOGONAL) {
    // x
    if (box._pbc_x) {
      if (px > box._coord_max[0]) {
        flag_px_tid += 1;
        px -= box._length[0];

      } else if (px < box._coord_min[0]) {
        flag_px_tid -= 1;
        px += box._length[0];
      }
    }

    // y
    if (box._pbc_y) {
      if (py > box._coord_max[1]) {
        flag_py_tid += 1;
        py -= box._length[1];
      } else if (py < box._coord_min[1]) {
        flag_py_tid -= 1;
        py += box._length[1];
      }
    }

    // z
    if (box._pbc_z) {
      if (pz > box._coord_max[2]) {
        flag_pz_tid += 1;
        pz -= box._length[2];
      } else if (pz < box._coord_min[2]) {
        flag_pz_tid -= 1;
        pz += box._length[2];
      }
    }
  }
}

__host__ __device__ __forceinline__ void ApplyPBC_unflag(
     Box box, rbmd::Real& px, rbmd::Real& py, rbmd::Real& pz) {
  if (box._type == Box::BoxType::ORTHOGONAL) {
    // x
    if (box._pbc_x) {
      if (px > box._coord_max[0]) {
        px -= box._length[0];

      } else if (px < box._coord_min[0]) {
        px += box._length[0];
      }
    }

    // y
    if (box._pbc_y) {
      if (py > box._coord_max[1]) {
        py -= box._length[1];
      } else if (py < box._coord_min[1]) {
        py += box._length[1];
      }
    }

    // z
    if (box._pbc_z) {
      if (pz > box._coord_max[2]) {
        pz -= box._length[2];
      } else if (pz < box._coord_min[2]) {
        pz += box._length[2];
      }
    }
  }
}

__host__ __device__ __forceinline__ rbmd::Real CalculateVolume( Box box) {
  if (box._type == Box::BoxType::ORTHOGONAL) {
    return box._length[0] * box._length[1] * box._length[2];
  } else {
    return 0;
    exit(0);  // todo
  }
}

