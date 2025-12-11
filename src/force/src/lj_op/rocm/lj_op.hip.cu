#include "../common/rbmd_define.h"
#include "lj_op.h"
#include "model/box.h"

namespace op {

  //---------device---------//
  __device__ rbmd::Id index2D(rbmd::Id row, rbmd::Id col, rbmd::Id num_columns) {
    return row * num_columns + col;
  }

  __device__ rbmd::Id index3D(rbmd::Id depth, rbmd::Id row, rbmd::Id col,
                              rbmd::Id num_rows, rbmd::Id num_columns) {
    return depth * (num_rows * num_columns) + row * num_columns + col;
  }

  // lj126
  inline __device__ void lj126(rbmd::Real cut_off, rbmd::Real px12, rbmd::Real py12,
                        rbmd::Real pz12, rbmd::Real eps_ij, rbmd::Real sigma_ij,
                        rbmd::Real& force_lj, rbmd::Real& energy_lj) {
    const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
    const rbmd::Real cut_off_2 = cut_off * cut_off;

    //if (dis_2 < cut_off_2 && dis_2 > EPSILON){
    if (dis_2 < cut_off_2) {
      rbmd::Real sigmaij_6 = POW(sigma_ij, 6.0);
      rbmd::Real dis_6 = POW(dis_2, 3.0);
      rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

      force_lj = -24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;//+
      energy_lj =
          0.5 * (4 * eps_ij * (sigmaij_6 / dis_6 - 1) * sigmaij_dis_6);
    } else {
      force_lj = 0.0;
      energy_lj = 0.0;
    }
  }

  // lj126_rs
  inline __device__ void lj126_rs(rbmd::Real rs, rbmd::Real px12, rbmd::Real py12,
                           rbmd::Real pz12, rbmd::Real eps_ij,
                           rbmd::Real sigma_ij, rbmd::Real& fs_ij) {
    const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
    const rbmd::Real rs_2 = rs * rs;

    if (dis_2 < rs_2 && dis_2 > EPSILON) {
      rbmd::Real sigmaij_6 = POW(sigma_ij, 6.0);
      rbmd::Real dis_6 = POW(dis_2, 3.0);
      rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;
      fs_ij = -24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;
    } else
      fs_ij = 0.0;
  }

  // lj126_rcs
  inline __device__ void lj126_rcs(rbmd::Real rc, rbmd::Real rs, rbmd::Id pice_num,
                            rbmd::Real px12, rbmd::Real py12, rbmd::Real pz12,
                            rbmd::Real eps_ij, rbmd::Real sigma_ij,
                            rbmd::Real& fcs_ij) {
    const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
    const rbmd::Real rc_2 = rc * rc;
    const rbmd::Real rs_2 = rs * rs;

    if (dis_2 < rc_2 && dis_2 > rs_2) {
      rbmd::Real sigmaij_6 = POW(sigma_ij, 6.0);
      rbmd::Real dis_6 = POW(dis_2, 3.0);
      rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

      fcs_ij = pice_num *
               (-24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2);
    } else
      fcs_ij = 0.0;
  }


  inline __device__ void ComputeVirial(rbmd::Real px12, rbmd::Real py12,
                        rbmd::Real pz12,rbmd::Real force,
                        rbmd::Real& local_virial_xx,rbmd::Real& local_virial_yy,
                        rbmd::Real& local_virial_zz,rbmd::Real& local_virial_xy,
                        rbmd::Real& local_virial_xz,rbmd::Real& local_virial_yz)
  {
    local_virial_xx = -0.5* px12 *px12 * force;
    local_virial_yy = -0.5* py12 *py12 * force;
    local_virial_zz = -0.5* pz12 *pz12 * force;
    local_virial_xy = -0.5* px12 *py12 * force;
    local_virial_xz = -0.5* px12 *pz12 * force;
    local_virial_yz = -0.5* py12 *pz12 * force;
  }

  //------global---------//
  // verlet-list: LJForce
  __global__ void ComputeLJForce(
       Box box, const rbmd::Real cut_off, const rbmd::Id num_atoms,
      const rbmd::Id* atoms_type, const rbmd::Real* sigma, const rbmd::Real* eps,
      const rbmd::Id* start_id, const rbmd::Id* end_id,
      const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
      const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
      rbmd::Real* flat_virial,rbmd::Real* total_evdwl) {
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage;

    rbmd::Real sum_fx = 0;
    rbmd::Real sum_fy = 0;
    rbmd::Real sum_fz = 0;

    //virial init
    rbmd::Real sum_virial[6];
    for (int i = 0; i < 6; ++i)
    {
      sum_virial[i] = 0.0;
    }
    rbmd::Real sum_elj = 0;

    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_atoms)
    {
      rbmd::Id typei = atoms_type[tid1];
      rbmd::Real eps_i = eps[typei];
      rbmd::Real sigma_i = sigma[typei];
      rbmd::Real x1 = px[tid1];
      rbmd::Real y1 = py[tid1];
      rbmd::Real z1 = pz[tid1];

      for (int j = start_id[tid1]; j < end_id[tid1]; ++j) {
        rbmd::Id tid2 = id_verletlist[j];
        rbmd::Id typej = atoms_type[tid2];
        rbmd::Real eps_j = eps[typej];
        rbmd::Real sigma_j = sigma[typej];
        // mix
        rbmd::Real eps_ij = SQRT(eps_i * eps_j);
        rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;
        rbmd::Real x2 = px[tid2];
        rbmd::Real y2 = py[tid2];
        rbmd::Real z2 = pz[tid2];
        rbmd::Real px12 = x2 - x1;
        rbmd::Real py12 = y2 - y1;
        rbmd::Real pz12 = z2 - z1;
        MinImageDistance(box, px12, py12, pz12);

        rbmd::Real force_lj;
        rbmd::Real energy_lj;
        rbmd::Real local_virial_xx,local_virial_yy,local_virial_zz,
        local_virial_xy,local_virial_xz,local_virial_yz;

        lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij, force_lj, energy_lj);

        ComputeVirial(px12, py12, pz12,force_lj,local_virial_xx,local_virial_yy,
    local_virial_zz,local_virial_xy,local_virial_xz,local_virial_yz);
        sum_fx += force_lj * px12;
        sum_fy += force_lj * py12;
        sum_fz += force_lj * pz12;
        sum_elj += energy_lj;

        //
        sum_virial[0] +=local_virial_xx;
        sum_virial[1] +=local_virial_yy;
        sum_virial[2] +=local_virial_zz;
        sum_virial[3] +=local_virial_xy;
        sum_virial[4] +=local_virial_xz;
        sum_virial[5] +=local_virial_yz;
      }
    }
      fx[tid1] = sum_fx;
      fy[tid1] = sum_fy;
      fz[tid1] = sum_fz;
      //
      for(int i =0;i<6;++i) {
        flat_virial[  i * num_atoms + tid1 ] = sum_virial[i];
      }

    rbmd::Real block_sum =
        BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage).Sum(sum_elj);
    if (threadIdx.x == 0) {
      atomicAdd(total_evdwl, block_sum);
    }
  }

  // RBL: LJForce
  __global__ void ComputeLJRBLForce(
      const Box  box, const rbmd::Real rs, const rbmd::Real rc,
      const rbmd::Id num_atoms, const rbmd::Id neighbor_sample_num,
      const rbmd::Id pice_num, const rbmd::Id* __restrict__ atoms_type,
      const rbmd::Real* __restrict__ sigma, const rbmd::Real* __restrict__ eps, const rbmd::Id* __restrict__ start_id,
      const rbmd::Id* __restrict__ end_id, const rbmd::Id* __restrict__ id_verletlist,
      const rbmd::Id* __restrict__ id_random_neighbor, const rbmd::Id* __restrict__ random_neighbor_num,
      const rbmd::Real* __restrict__ px, const rbmd::Real* __restrict__ py, const rbmd::Real* __restrict__ pz,
      rbmd::Real* __restrict__ fx, rbmd::Real* __restrict__ fy, rbmd::Real* __restrict__ fz) {
    rbmd::Real sum_fx = 0;
    rbmd::Real sum_fy = 0;
    rbmd::Real sum_fz = 0;

    rbmd::Real sum_fsx = 0;
    rbmd::Real sum_fsy = 0;
    rbmd::Real sum_fsz = 0;

    rbmd::Real sum_fcsx = 0;
    rbmd::Real sum_fcsy = 0;
    rbmd::Real sum_fcsz = 0;

    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_atoms) {
      __shared__ rbmd::Real shared_px[BLOCK_SIZE];
      __shared__ rbmd::Real shared_py[BLOCK_SIZE];
      __shared__ rbmd::Real shared_pz[BLOCK_SIZE];
      if (threadIdx.x < BLOCK_SIZE) {
        shared_px[threadIdx.x] = __ldg(&px[tid1]);
        shared_py[threadIdx.x] = __ldg(&py[tid1]);
        shared_pz[threadIdx.x] = __ldg(&pz[tid1]);
      }
      __syncthreads();
      rbmd::Id typei = __ldg(&atoms_type[tid1]);
      rbmd::Real eps_i = __ldg(&eps[typei]);
      rbmd::Real sigma_i = __ldg(&sigma[typei]);
      rbmd::Real x1 = shared_px[threadIdx.x];
      rbmd::Real y1 = shared_py[threadIdx.x];
      rbmd::Real z1 = shared_pz[threadIdx.x];

      rbmd::Real fs_ij, fcs_ij;
      // rs
      for (rbmd::Id j = __ldg(&start_id[tid1]); j < __ldg(&end_id[tid1]); ++j) {
        rbmd::Id tid2 = __ldg(&id_verletlist[j]);
        rbmd::Id typej = __ldg(&atoms_type[tid2]);
        rbmd::Real eps_j = __ldg(&eps[typej]);
        rbmd::Real sigma_j = __ldg(&sigma[typej]);

        // mix
        rbmd::Real eps_ij = SQRT(eps_i * eps_j);
        rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

        rbmd::Real x2 = __ldg(&px[tid2]);
        rbmd::Real y2 = __ldg(&py[tid2]);
        rbmd::Real z2 = __ldg(&pz[tid2]);
        rbmd::Real px12 = x2 - x1;
        rbmd::Real py12 = y2 - y1;
        rbmd::Real pz12 = z2 - z1;

        MinImageDistance(box, px12, py12, pz12);

        // compute the force_rs
        lj126_rs(rs, px12, py12, pz12, eps_ij, sigma_ij, fs_ij);

        sum_fsx += fs_ij * px12;
        sum_fsy += fs_ij * py12;
        sum_fsz += fs_ij * pz12;
      }

      // rcs
      rbmd::Id real_random_num = __ldg(&random_neighbor_num[tid1]);
      for (rbmd::Id jj = 0; jj < real_random_num; ++jj) {
        rbmd::Id tid2 =
            __ldg(&id_random_neighbor[tid1 * neighbor_sample_num + jj]);
        rbmd::Id typej = __ldg(&atoms_type[tid2]);
        rbmd::Real eps_j = __ldg(&eps[typej]);
        rbmd::Real sigma_j = __ldg(&sigma[typej]);

        // mix
        rbmd::Real eps_ij = SQRT(eps_i * eps_j);
        rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

        rbmd::Real x2 = __ldg(&px[tid2]);
        rbmd::Real y2 = __ldg(&py[tid2]);
        rbmd::Real z2 = __ldg(&pz[tid2]);
        rbmd::Real px12 = x2 - x1;
        rbmd::Real py12 = y2 - y1;
        rbmd::Real pz12 = z2 - z1;
        MinImageDistance(box, px12, py12, pz12);

        // compute the force_rcs
        lj126_rcs(rc, rs, pice_num, px12, py12, pz12, eps_ij, sigma_ij, fcs_ij);

        sum_fcsx += fcs_ij * px12;
        sum_fcsy += fcs_ij * py12;
        sum_fcsz += fcs_ij * pz12;
      }

      // total force = fs + fcs
      sum_fx = sum_fsx + sum_fcsx;
      sum_fy = sum_fsy + sum_fcsy;
      sum_fz = sum_fsz + sum_fcsz;

      fx[tid1] = sum_fx;
      fy[tid1] = sum_fy;
      fz[tid1] = sum_fz;
    }
  }

  // verlet-list: LJEnergy
  __global__ void ComputeLJEnergy(
       Box box, const rbmd::Real cut_off, const rbmd::Id num_atoms,
      const rbmd::Id* atoms_type, const rbmd::Real* sigma, const rbmd::Real* eps,
      const rbmd::Id* start_id, const rbmd::Id* end_id,
      const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
      const rbmd::Real* pz, rbmd::Real* flat_virial,rbmd::Real* total_evdwl) {
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage;
    rbmd::Real sum_elj = 0;

    //virial init
    rbmd::Real sum_virial[6];
    for (int i = 0; i < 6; ++i)
    {
      sum_virial[i] = 0.0;
    }

    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_atoms) {
      rbmd::Id typei = atoms_type[tid1];
      rbmd::Real eps_i = eps[typei];
      rbmd::Real sigma_i = sigma[typei];
      rbmd::Real x1 = px[tid1];
      rbmd::Real y1 = py[tid1];
      rbmd::Real z1 = pz[tid1];
      for (int j = start_id[tid1]; j < end_id[tid1]; ++j) {
        rbmd::Id tid2 = id_verletlist[j];
        rbmd::Id typej = atoms_type[tid2];
        rbmd::Real eps_j = eps[typej];
        rbmd::Real sigma_j = sigma[typej];
        // mix
        rbmd::Real eps_ij = SQRT(eps_i * eps_j);
        rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;
        rbmd::Real x2 = px[tid2];
        rbmd::Real y2 = py[tid2];
        rbmd::Real z2 = pz[tid2];
        rbmd::Real px12 = x2 - x1;
        rbmd::Real py12 = y2 - y1;
        rbmd::Real pz12 = z2 - z1;
        MinImageDistance(box, px12, py12, pz12);
        rbmd::Real force_lj;
        rbmd::Real energy_ij;
        lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij, force_lj, energy_ij);
        sum_elj += energy_ij;

        rbmd::Real local_virial_xx,local_virial_yy,local_virial_zz,
    local_virial_xy,local_virial_xz,local_virial_yz;
        ComputeVirial(px12, py12, pz12,force_lj,local_virial_xx,local_virial_yy,
      local_virial_zz,local_virial_xy,local_virial_xz,local_virial_yz);

        sum_virial[0] +=local_virial_xx;
        sum_virial[1] +=local_virial_yy;
        sum_virial[2] +=local_virial_zz;
        sum_virial[3] +=local_virial_xy;
        sum_virial[4] +=local_virial_xz;
        sum_virial[5] +=local_virial_yz;
      }
      //
      for(int i =0;i<6;++i) {
        flat_virial[ i * num_atoms +  tid1] = sum_virial[i];
      }
    }

    rbmd::Real block_sum =
        BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage).Sum(sum_elj);
    if (threadIdx.x == 0) {
      atomicAdd(total_evdwl, block_sum);
    }
  }

  /////////////////////
  // verlet-list: LJForce
  void LJForceOp<device::DEVICE_GPU>::operator()(
       Box box, const rbmd::Real cut_off, const rbmd::Id num_atoms,
      const rbmd::Id* atoms_type, const rbmd::Real* sigma, const rbmd::Real* eps,
      const rbmd::Id* start_id, const rbmd::Id* end_id,
      const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
      const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
      rbmd::Real* flat_virial ,rbmd::Real* total_evdwl) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeLJForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, cut_off, num_atoms, atoms_type, sigma, eps, start_id, end_id,
        id_verletlist, px, py, pz, fx, fy, fz, flat_virial,total_evdwl));
  }

  // RBL:  LJLForce
  void LJRBLForceOp<device::DEVICE_GPU>::operator()(
      const  Box box, const rbmd::Real rs, const rbmd::Real rc,
      const rbmd::Id num_atoms, const rbmd::Id neighbor_sample_num,
      const rbmd::Id pice_num, const rbmd::Id* atoms_type,
      const rbmd::Real* sigma, const rbmd::Real* eps, const rbmd::Id* start_id,
      const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
      const rbmd::Id* id_random_neighbor, const rbmd::Id* random_neighbor_num,
      const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
      rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeLJRBLForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, rs, rc, num_atoms, neighbor_sample_num, pice_num, atoms_type, sigma,
        eps, start_id, end_id, id_verletlist, id_random_neighbor,
        random_neighbor_num, px, py, pz, fx, fy, fz));
  }

  // verlet-list: LJEnergy
  void LJEnergyOp<device::DEVICE_GPU>::operator()(
       Box box, const rbmd::Real cut_off, const rbmd::Id num_atoms,
      const rbmd::Id* atoms_type, const rbmd::Real* sigma, const rbmd::Real* eps,
      const rbmd::Id* start_id, const rbmd::Id* end_id,
      const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
      const rbmd::Real* pz, rbmd::Real* flat_virial,rbmd::Real* total_evdwl) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeLJEnergy<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, cut_off, num_atoms, atoms_type, sigma, eps, start_id, end_id,
        id_verletlist, px, py, pz, flat_virial,total_evdwl));
  }
}

