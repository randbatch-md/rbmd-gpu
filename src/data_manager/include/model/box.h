#pragma once
#include "../../../common/rbmd_define.h"
#include "../../../common/types.h"

class Box {
 public:
  //
  enum class BoxType {
    ORTHOGONAL,  //
    TRICLINIC    //
  };

  BoxType _type = BoxType::ORTHOGONAL;  //default
  Box() = default;

  /**
   * initialize the box
   * @param box_type
   * @param coord_min
   * @param coord_max
   * @param pbc
   */
  void Setup(BoxType box_type, const rbmd::Real coord_min[3],
            const rbmd::Real coord_max[3], const bool pbc[3]);

  // TODO
  ///
  bool _pbc_x = true;
  ///
  bool _pbc_y = true;
  ///
  bool _pbc_z = true;

  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 6)) _length[6]{};
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 6)) _length_inv[6]{};
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) _median_point[3]{};

  /// （local）The coordinates of the lower left corner of the box
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) _coord_min[3]{};
  /// （local）The coordinates of the upper right corner of the box
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) _coord_max[3]{};

  /// individual cells
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

static __host__ __device__ __forceinline__ void MinImageDistance_while( Box box,
                                                          rbmd::Real& dx,
                                                          rbmd::Real& dy,
                                                          rbmd::Real& dz) {
  if (box._type == Box::BoxType::ORTHOGONAL) {
    if (box._pbc_x) {
      while (ABS(dx) > box._length[0] * 0.5) {
        dx -= (dx > 0 ? box._length[0] : -box._length[0]);
      }
    }
    if (box._pbc_y) {
      while (ABS(dy) > box._length[1] * 0.5) {
        dy -= (dy > 0 ? box._length[1] : -box._length[1]);
      }
    }
    if (box._pbc_z) {
      while (ABS(dz) > box._length[2] * 0.5) {
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

__host__ __device__ __forceinline__ void ApplyPBC_Robust(
    const Box& box,
    rbmd::Real& px, rbmd::Real& py, rbmd::Real& pz,
    rbmd::Id& flag_px, rbmd::Id& flag_py, rbmd::Id& flag_pz)
{
  if (box._type == Box::BoxType::ORTHOGONAL) {
    // X
    if (box._pbc_x) {
      // The floor function directly calculates how many "box units"
      // the atomic coordinates are away from the lower boundary of the box.
      rbmd::Real wrapped_px = px - box._coord_min[0];
      rbmd::Id crossings = FLOOR(wrapped_px * box._length_inv[0]);

      px -= crossings * box._length[0];
      flag_px -= crossings;
    }

    // Y
    if (box._pbc_y) {
      rbmd::Real wrapped_py = py - box._coord_min[1];
      rbmd::Id crossings = FLOOR(wrapped_py * box._length_inv[1]);

      py -= crossings * box._length[1];
      flag_py -= crossings;
    }

    // Z
    if (box._pbc_z) {
      rbmd::Real wrapped_pz = pz - box._coord_min[2];
      rbmd::Id crossings = FLOOR(wrapped_pz * box._length_inv[2]);

      pz -= crossings * box._length[2];
      flag_pz -= crossings;
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

