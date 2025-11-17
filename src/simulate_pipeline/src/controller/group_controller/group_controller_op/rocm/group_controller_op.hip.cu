#include "rbmd_define.h"
#include "group_controller_op.h"

namespace op {
#define THREADS_PER_BLOCK 256

__device__ inline void unwarp(const rbmd::Real* pos, const int* image,
  Box box, rbmd::Real* unwrap_pos)
{

  unwrap_pos[0] = pos[0] + image[0] * box._length[0];
  unwrap_pos[1] = pos[1] + image[1] * box._length[1];
  unwrap_pos[2] = pos[2] + image[2] * box._length[2];
}

// 新增的VCM计算内核 (并行归约)
__global__ void compute_vcm_kernel(const int num_atoms, const int* atoms_type,
    const  rbmd::Real* mass,const  rbmd::Real* vx, const  rbmd::Real* vy,
    const  rbmd::Real* vz,rbmd::Real* vcm_contrib)
{
    extern __shared__  rbmd::Real sdata[];
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // sdata[0..3] -> mass, mom_x, mom_y, mom_z
    sdata[tid] = 0.0;
    sdata[tid + blockDim.x] = 0.0;
    sdata[tid + 2 * blockDim.x] = 0.0;
    sdata[tid + 3 * blockDim.x] = 0.0;

    if (i < num_atoms) {
        rbmd::Real m = mass[atoms_type[i]];
        sdata[tid] = m;
        sdata[tid + blockDim.x] = m * vx[i];
        sdata[tid + 2 * blockDim.x] = m * vy[i];
        sdata[tid + 3 * blockDim.x] = m * vz[i];
    }
    __syncthreads();

    // 并行归约
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
            sdata[tid + blockDim.x] += sdata[tid + blockDim.x + s];
            sdata[tid + 2 * blockDim.x] += sdata[tid + 2 * blockDim.x + s];
            sdata[tid + 3 * blockDim.x] += sdata[tid + 3 * blockDim.x + s];
        }
        __syncthreads();
    }

    // 线程0将块内结果原子加到全局变量
    if (tid == 0) {
        atomicAdd(&vcm_contrib[0], sdata[0]);
        atomicAdd(&vcm_contrib[1], sdata[blockDim.x]);
        atomicAdd(&vcm_contrib[2], sdata[2 * blockDim.x]);
        atomicAdd(&vcm_contrib[3], sdata[3 * blockDim.x]);
    }
}

__global__ void compute_vcm_kernel1(
    const int num_atoms,const int* atoms_type,const rbmd::Real* mass,
    const rbmd::Real* vx,const rbmd::Real* vy,const rbmd::Real* vz,
    rbmd::Real* vcm_contrib)
{
  // 每个线程块需要 4 个归约存储（质量 + 动量 x/y/z）
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_mass;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_mv_x;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_mv_y;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_mv_z;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  // 初始化局部变量
  rbmd::Real local_mass = 0.0;
  rbmd::Real local_mv_x = 0.0;
  rbmd::Real local_mv_y = 0.0;
  rbmd::Real local_mv_z = 0.0;

  // 计算当前线程处理的原子（如果有效）
  if (tid < num_atoms) {
    rbmd::Real m = mass[atoms_type[tid]];  // 原子质量
    local_mass = m;
    local_mv_x = m * vx[tid];  // 动量 x
    local_mv_y = m * vy[tid];  // 动量 y
    local_mv_z = m * vz[tid];  // 动量 z
  }

  // 块内归约（分别计算总质量、总动量 x/y/z）
  rbmd::Real block_mass = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
   (temp_storage_mass).Sum(local_mass);
  rbmd::Real block_mv_x = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
    (temp_storage_mv_x).Sum(local_mv_x);
  rbmd::Real block_mv_y = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
   (temp_storage_mv_y).Sum(local_mv_y);
  rbmd::Real block_mv_z = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
    (temp_storage_mv_z).Sum(local_mv_z);

  // 线程 0 将块内结果原子性地累加到全局变量
  if (threadIdx.x == 0) {
    atomicAdd(&vcm_contrib[0], block_mass);  // 总质量
    atomicAdd(&vcm_contrib[1], block_mv_x);    // 总动量 x
    atomicAdd(&vcm_contrib[2], block_mv_y);    // 总动量 y
    atomicAdd(&vcm_contrib[3], block_mv_z);    // 总动量 z
  }
}


__global__ void unwarp_position(const rbmd::Id num_atoms, Box  box,
   const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
   const rbmd::Id* flag_px, const rbmd::Id* flag_py, const rbmd::Id* flag_pz,
   rbmd::Real* unwarp_px,rbmd::Real* unwarp_py,rbmd::Real* unwarp_pz)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < num_atoms)
  {
    unwarp_px[i] = px[i] + flag_px[i] * box._length[0];
    unwarp_py[i] = py[i] + flag_py[i] * box._length[1];
    unwarp_pz[i] = pz[i] + flag_pz[i] * box._length[2];
  }
}

/**
 * @brief GPU内核：计算质心位置(XCM)所需的质量和质量矩
 */
__global__ void compute_xcm_kernel(const rbmd::Id num_atoms,
  const rbmd::Real* mass,const rbmd::Id* atoms_type,
  const rbmd::Real* unwrap_px, const rbmd::Real* unwrap_py,
  const rbmd::Real* unwrap_pz,rbmd::Real* xcm_contrib)
{
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_mass;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_mp_x;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_mp_y;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_mp_z;

  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  rbmd::Real local_mass = 0.0;
  rbmd::Real local_mp_x = 0.0;
  rbmd::Real local_mp_y = 0.0;
  rbmd::Real local_mp_z = 0.0;

  if (tid < num_atoms)
  {
    local_mass  = mass[atoms_type[tid]];
    local_mp_x  = local_mass * unwrap_px[tid];  //m * p_x
    local_mp_y  = local_mass * unwrap_py[tid]; //m * p_y
    local_mp_z  = local_mass * unwrap_pz[tid]; //m * p_z
  }

  rbmd::Real block_mass = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
   (temp_storage_mass).Sum(local_mass);
  rbmd::Real block_mp_x = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
    (temp_storage_mp_x).Sum(local_mp_x);
  rbmd::Real block_mp_y = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
   (temp_storage_mp_y).Sum(local_mp_y);
  rbmd::Real block_mp_z = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
    (temp_storage_mp_z).Sum(local_mp_z);

  if (threadIdx.x == 0) {
    atomicAdd(&xcm_contrib[0], block_mass);
    atomicAdd(&xcm_contrib[1], block_mp_x);
    atomicAdd(&xcm_contrib[2], block_mp_y);
    atomicAdd(&xcm_contrib[3], block_mp_z);
  }
}


/**
 * @brief GPU内核：高效计算角动量
 */
__global__ void compute_angmom_kernel(const rbmd::Id num_atoms,Real3 cm,
                                              const int* d_atoms_type, const rbmd::Real* d_mass,
                                              const rbmd::Real* d_x, const rbmd::Real* d_y, const rbmd::Real* d_z,
                                              const rbmd::Real* d_vx, const rbmd::Real* d_vy, const rbmd::Real* d_vz,
                                              const int* d_image_x, const int* d_image_y, const int* d_image_z,
                                              Box  box,rbmd::Real* angmom_contrib) {
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_angmom_x;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_angmom_y;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_angmom_z;

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  rbmd::Real local_angmom_x = 0.0;
  rbmd::Real local_angmom_y = 0.0;
  rbmd::Real local_angmom_z = 0.0;

  if (i < num_atoms)
  {
    rbmd::Real pos_wrap[3] = {d_x[i], d_y[i], d_z[i]};
    int image[3] = {d_image_x[i], d_image_y[i], d_image_z[i]};
    rbmd::Real pos_unwrap[3];

    unwarp(pos_wrap, image, box, pos_unwrap);
     // printf("pos_unwrap: %f %f %f\n",pos_unwrap[0],pos_unwrap[1],pos_unwrap[2]);

    rbmd::Real dx = pos_unwrap[0] - cm.x;
    rbmd::Real dy = pos_unwrap[1] - cm.y;
    rbmd::Real dz = pos_unwrap[2] - cm.z;
    rbmd::Real mass_i = d_mass[d_atoms_type[i]];
    // L = m * (r x v)
    local_angmom_x = mass_i * (dy * d_vz[i] - dz * d_vy[i]);  // Lx
    local_angmom_y = mass_i * (dz * d_vx[i] - dx * d_vz[i]); // Ly
    local_angmom_z = mass_i * (dx * d_vy[i] - dy * d_vx[i]); // Lz

  }

  rbmd::Real block_angmom_x = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
    (temp_storage_angmom_x).Sum(local_angmom_x);
  rbmd::Real block_angmom_y = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
   (temp_storage_angmom_y).Sum(local_angmom_y);
  rbmd::Real block_angmom_z = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
    (temp_storage_angmom_z).Sum(local_angmom_z);

  if (threadIdx.x == 0) {
    atomicAdd(&angmom_contrib[0], block_angmom_x);
    atomicAdd(&angmom_contrib[1], block_angmom_y);
    atomicAdd(&angmom_contrib[2], block_angmom_z);
  }
}

/**
 * @brief GPU内核：高效计算转动惯量张量
 */
__global__ void compute_inertia_kernel(const rbmd::Id num_atoms,Real3 cm,
                                              const int* d_atoms_type, const rbmd::Real* d_mass,
                                              const rbmd::Real* d_x, const rbmd::Real* d_y, const rbmd::Real* d_z,
                                              const int* d_image_x, const int* d_image_y, const int* d_image_z,
                                              Box  box,rbmd::Real* inertia_contrib) {

  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_inertia_xx;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_inertia_yy;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_inertia_zz;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_inertia_xy;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_inertia_xz;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_inertia_yz;

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  rbmd::Real local_inertia_xx = 0.0;
  rbmd::Real local_inertia_yy = 0.0;
  rbmd::Real local_inertia_zz = 0.0;
  rbmd::Real local_inertia_xy = 0.0;
  rbmd::Real local_inertia_xz = 0.0;
  rbmd::Real local_inertia_yz = 0.0;

  if (i < num_atoms)
  {
    rbmd::Real pos_wrap[3] = {d_x[i], d_y[i], d_z[i]};
    int image[3] = {d_image_x[i], d_image_y[i], d_image_z[i]};
    rbmd::Real pos_unwrap[3];
    unwarp(pos_wrap, image, box, pos_unwrap);

    rbmd::Real dx = pos_unwrap[0] - cm.x;
    rbmd::Real dy = pos_unwrap[1] - cm.y;
    rbmd::Real dz = pos_unwrap[2] - cm.z;
    rbmd::Real mass_i = d_mass[d_atoms_type[i]];

    // Ixx = m*(y^2+z^2), Ixy = -m*x*y, etc.
    local_inertia_xx = mass_i * (dy * dy + dz * dz);  // Ixx
    local_inertia_yy = mass_i * (dx * dx + dz * dz);  // Iyy
    local_inertia_zz = mass_i * (dx * dx + dy * dy);  // Izz
    local_inertia_xy = -mass_i * dx * dy;             // Ixy
    local_inertia_xz = -mass_i * dx * dz;             // Ixz
    local_inertia_yz = -mass_i * dy * dz;             // Iyz
  }

  rbmd::Real block_inertia_xx = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
  (temp_storage_inertia_xx).Sum(local_inertia_xx);
  rbmd::Real block_inertia_yy = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
   (temp_storage_inertia_yy).Sum(local_inertia_yy);
  rbmd::Real block_inertia_zz = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
    (temp_storage_inertia_zz).Sum(local_inertia_zz);

  rbmd::Real block_inertia_xy = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
  (temp_storage_inertia_xy).Sum(local_inertia_xy);
  rbmd::Real block_inertia_xz = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
   (temp_storage_inertia_xz).Sum(local_inertia_xz);
  rbmd::Real block_inertia_yz = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
    (temp_storage_inertia_yz).Sum(local_inertia_yz);

  if (threadIdx.x == 0) {
    atomicAdd(&inertia_contrib[0], block_inertia_xx);
    atomicAdd(&inertia_contrib[1], block_inertia_yy);
    atomicAdd(&inertia_contrib[2], block_inertia_zz);
    atomicAdd(&inertia_contrib[3], block_inertia_xy);
    atomicAdd(&inertia_contrib[4], block_inertia_xz);
    atomicAdd(&inertia_contrib[5], block_inertia_yz);
  }

}



// =============================================================================
//  (Functor Implementations)
// =============================================================================

void ComputeVCMOp<device::DEVICE_GPU>::operator()(const int num_atoms, const int* atoms_type,
    const  rbmd::Real* mass,const  rbmd::Real* vx, const  rbmd::Real* vy,
    const  rbmd::Real* vz,rbmd::Real* vcm_contrib) {
  // ... 计算blocks和threads ...
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

  // 需要为每个归约变量分配共享内存
//   CHECK_KERNEL(compute_vcm_kernel<<<blocks, threads, 4 * threads * sizeof(rbmd::Real)>>>
// (num_atoms,atoms_type,mass,vx,vy,vz,vcm_contrib));
  CHECK_KERNEL(compute_vcm_kernel1<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>
  (num_atoms,atoms_type,mass,vx,vy,vz,vcm_contrib));
}

void UnwarpPositionOp<device::DEVICE_GPU>::operator()(
   const rbmd::Id num_atoms, Box  box,
   const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
   const rbmd::Id* flag_px, const rbmd::Id* flag_py, const rbmd::Id* flag_pz,
   rbmd::Real* unwarp_px,rbmd::Real* unwarp_py,rbmd::Real* unwarp_pz) {

  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(unwarp_position<<<blocks_per_grid, BLOCK_SIZE,0,0>>>
    (num_atoms,box, px,py,pz,flag_px,flag_py,flag_pz,
      unwarp_px, unwarp_py, unwarp_pz));
}


void ComputeXCMOp<device::DEVICE_GPU>::operator()(
  const rbmd::Id num_atoms,
  const rbmd::Real* mass,const rbmd::Id* atoms_type,
  const rbmd::Real* unwrap_px, const rbmd::Real* unwrap_py,
  const rbmd::Real* unwrap_pz,rbmd::Real* xcm_contrib) {

  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(compute_xcm_kernel<<<blocks_per_grid, BLOCK_SIZE,0,0>>>(num_atoms,
    mass,atoms_type,unwrap_px, unwrap_py, unwrap_pz,xcm_contrib));
}

void ComputeAngMomOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, Real3 cm,
    const rbmd::Id* d_atoms_type, const rbmd::Real* d_mass,
    const rbmd::Real* d_x, const rbmd::Real* d_y, const rbmd::Real* d_z,
    const rbmd::Real* d_vx, const rbmd::Real* d_vy, const rbmd::Real* d_vz,
    const int* d_image_x, const int* d_image_y, const int* d_image_z,
     Box box,rbmd::Real* d_angmom_contrib) {

  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(compute_angmom_kernel<<<blocks_per_grid, BLOCK_SIZE,0,0>>>(num_atoms, cm,
                                                     d_atoms_type, d_mass,
                                                     d_x, d_y, d_z, d_vx, d_vy, d_vz,
                                                     d_image_x, d_image_y, d_image_z,
                                                     box, d_angmom_contrib));

}

void ComputeInertiaOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, Real3 cm,
    const rbmd::Id* d_atoms_type, const rbmd::Real* d_mass,
    const rbmd::Real* d_x, const rbmd::Real* d_y, const rbmd::Real* d_z,
    const int* d_image_x, const int* d_image_y, const int* d_image_z,
     Box box,rbmd::Real* d_inertia_contrib) {

  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(compute_inertia_kernel<<<blocks_per_grid, BLOCK_SIZE,0,0>>>(num_atoms, cm,
                                                     d_atoms_type, d_mass,
                                                     d_x, d_y, d_z,d_image_x, d_image_y, d_image_z,
                                                     box, d_inertia_contrib));
}

}