//#include <hip/hip_runtime.h>

#include "../common/rbmd_define.h"
#include "cvff_op.h"
#include "model/box.h"
#include "../lj_cut_coul_kspace_op/rocm/lj_cut_coul_kspace_op.hip.cu"


namespace op{
const rbmd::Real SMALL = 0.001;

//verlet-list : SpecialLJCutCoul
__global__ void ComputeSpecialLJCutCoulForce(
     Box box, ERFTable* erf_table, const rbmd::Real cut_off,
    const rbmd::Id num_atoms, const rbmd::Real alpha, const rbmd::Real qqr2e,
    const rbmd::Id* atoms_type, const rbmd::Id* atoms_id,
    const rbmd::Real* sigma, const rbmd::Real* eps, const rbmd::Id* start_id,
    const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
    const rbmd::Id* special_ids, const rbmd::Real* special_weights,
    const rbmd::Id* special_offset, const rbmd::Id* special_count,
    const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial,rbmd::Real* total_evdwl, rbmd::Real* total_ecoul) {
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
      temp_storage_elj;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
      temp_storage_ecoul;

  rbmd::Real sum_fx = 0;
  rbmd::Real sum_fy = 0;
  rbmd::Real sum_fz = 0;
  rbmd::Real sum_elj = 0;
  rbmd::Real sum_ecoul = 0;
  //virial init
  rbmd::Real sum_virial[6];
  for (int i = 0; i < 6; ++i)
  {
    sum_virial[i] = 0.0;
  }

  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_atoms) {
    rbmd::Id atom_id1 = atoms_id[tid1];
    rbmd::Id num_components = special_offset[atom_id1];
    rbmd::Id typei = atoms_type[tid1];
    rbmd::Real eps_i = eps[typei];
    rbmd::Real sigma_i = sigma[typei];
    rbmd::Real charge_i = charge[tid1];
    rbmd::Real x1 = px[tid1];
    rbmd::Real y1 = py[tid1];
    rbmd::Real z1 = pz[tid1];

    for (int j = start_id[tid1]; j < end_id[tid1]; ++j) {
      rbmd::Id tid2 = id_verletlist[j];
      rbmd::Id atom_id2 = atoms_id[tid2];
      rbmd::Id typej = atoms_type[tid2];
      rbmd::Real eps_j = eps[typej];
      rbmd::Real sigma_j = sigma[typej];
      rbmd::Real charge_j = charge[tid2];
      // mix
      rbmd::Real eps_ij = SQRT(eps_i * eps_j);
      rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;
      rbmd::Real x2 = px[tid2];
      rbmd::Real y2 = py[tid2];
      rbmd::Real z2 = pz[tid2];
      rbmd::Real x12 = x2 - x1;
      rbmd::Real y12 = y2 - y1;
      rbmd::Real z12 = z2 - z1;
      // rbmd::Real x12 = x1 - x2;
      // rbmd::Real y12 = y1 - y2;
      // rbmd::Real z12 = z1 - z2;
      MinImageDistance(box, x12, y12, z12);
      // erf value
      rbmd::Real dis = SQRT(x12 * x12 + y12 * y12 + z12 * z12);
      rbmd::Id index_table_pij = Extract(dis);
      rbmd::Real table_pij = TableGnearValue(erf_table,dis, index_table_pij);

      rbmd::Real force_lj, force_coul, force_pair,force_coul_factor;
      rbmd::Real energy_lj, energy_coul,energy_coul_factor;
      // lj cut
      lj126(cut_off, x12, y12, z12, eps_ij, sigma_ij, force_lj, energy_lj);

      // special lj  weight
      rbmd::Real weight = 1.0;
      for (rbmd::Id k = 0; k < special_count[atom_id1]; ++k) {
        rbmd::Id special_id = special_ids[num_components + k];
        if (special_id == atom_id2) {
          weight = special_weights[num_components + k];
          // printf("weight %f\n", weight);
        }
      }

      // Coul cut
      CoulCutForce_fix(cut_off, alpha, qqr2e, charge_i, charge_j,
      x12, y12, z12, force_coul_factor,energy_coul_factor,
      force_coul, energy_coul);
      force_coul = force_coul-(1-weight)*force_coul_factor;
      energy_coul = energy_coul-(1-weight)*energy_coul_factor;

      // sum force of special_lj_cut  and coul_cut
      force_pair = weight * force_lj + force_coul;
      sum_fx += force_pair * x12;
      sum_fy += force_pair * y12;
      sum_fz += force_pair * z12;
      // sum energy  of special_lj_cut  and coul_cut
      sum_elj += weight * energy_lj;
      sum_ecoul += energy_coul;

      //rbmd::Real local_virial[6];
      rbmd::Real local_virial_xx,local_virial_yy,local_virial_zz,
        local_virial_xy,local_virial_xz,local_virial_yz;
      ComputeVirial_fix(x12, y12, z12,force_pair,
        local_virial_xx,local_virial_yy,local_virial_zz,
        local_virial_xy,local_virial_xz,local_virial_yz);

      //
      sum_virial[0] +=local_virial_xx;
      sum_virial[1] +=local_virial_yy;
      sum_virial[2] +=local_virial_zz;
      sum_virial[3] +=local_virial_xy;
      sum_virial[4] +=local_virial_xz;
      sum_virial[5] +=local_virial_yz;
    }

    fx[tid1] = sum_fx;
    fy[tid1] = sum_fy;
    fz[tid1] = sum_fz;
    //
    for(int i =0;i<6;++i) {
      flat_virial[  i * num_atoms + tid1] = sum_virial[i];
    }
  }

  rbmd::Real block_sum_elj =
      BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage_elj)
          .Sum(sum_elj);
  rbmd::Real block_sum_ecoul =
      BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage_ecoul)
          .Sum(sum_ecoul);

  if (threadIdx.x == 0) {
    atomicAdd(total_evdwl, block_sum_elj);
    atomicAdd(total_ecoul, block_sum_ecoul);
  }
}

//RBL : SpecialLJCutCoul
__global__ void ComputeSpecialLJCutCoulRBLForce(
     Box box, ERFTable* erf_table, const rbmd::Real rs, const rbmd::Real rc,
    const rbmd::Id num_atoms, const rbmd::Id neighbor_sample_num,
    const rbmd::Id pice_num, const rbmd::Real alpha, const rbmd::Real qqr2e,
    const rbmd::Id* atoms_type, const rbmd::Id* atoms_id,
    const rbmd::Real* sigma, const rbmd::Real* eps, const rbmd::Id* start_id,
    const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
    const rbmd::Id* id_random_neighbor, const rbmd::Id* random_neighbor_num,
    const rbmd::Id* special_ids, const rbmd::Real* special_weights,
    const rbmd::Id* special_offset, const rbmd::Id* special_count,
    const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz) {
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
    rbmd::Id atom_id1 = atoms_id[tid1];
    rbmd::Id num_components = special_offset[atom_id1];
    rbmd::Id typei = atoms_type[tid1];
    rbmd::Real eps_i = eps[typei];
    rbmd::Real sigma_i = sigma[typei];
    rbmd::Real charge_i = charge[tid1];
    rbmd::Real x1 = px[tid1];
    rbmd::Real y1 = py[tid1];
    rbmd::Real z1 = pz[tid1];
    // rs
    for (rbmd::Id j = start_id[tid1]; j < end_id[tid1]; ++j) {
      rbmd::Id tid2 = id_verletlist[j];
      rbmd::Id atom_id2 = atoms_id[tid2];
      rbmd::Id typej = atoms_type[tid2];
      rbmd::Real eps_j = eps[typej];
      rbmd::Real sigma_j = sigma[typej];
      rbmd::Real charge_j = charge[tid2];
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

      // erf
      rbmd::Real dis = SQRT(px12 * px12 + py12 * py12 + pz12 * pz12);
      rbmd::Id index_table_pij = Extract(dis);
      rbmd::Real table_pij = TableGnearValue(erf_table,dis, index_table_pij);

      // compute the force_rs
      rbmd::Real force_lj_rs, force_coul_rs,force_coul_factor_rs;
      rbmd::Real fs_ij;
      lj126_rs(rs, px12, py12, pz12, eps_ij, sigma_ij, force_lj_rs);

      rbmd::Real weight = 1.0;
      for (rbmd::Id k = 0; k < special_count[atom_id1]; ++k) {
        rbmd::Id special_id = special_ids[num_components + k];
        if (special_id == atom_id2) {
          weight = special_weights[num_components + k];
          // printf("weight %f\n", weight);
        }
      }
      force_lj_rs = weight * force_lj_rs;

      CoulCutForce_rs_fix(rs,alpha,qqr2e,charge_i,charge_j,
        px12, py12, pz12, force_coul_factor_rs,force_coul_rs);
      force_coul_rs = force_coul_rs - (1-weight)*force_coul_factor_rs;

      fs_ij = force_lj_rs + force_coul_rs;
      sum_fsx += fs_ij * px12;
      sum_fsy += fs_ij * py12;
      sum_fsz += fs_ij * pz12;
    }

    // rcs
    rbmd::Id real_random_num = random_neighbor_num[tid1];
    for (rbmd::Id jj = 0; jj < real_random_num; ++jj) {
      rbmd::Id tid2 = id_random_neighbor[tid1 * neighbor_sample_num + jj];
      rbmd::Id atom_id2 = atoms_id[tid2];
      rbmd::Id typej = atoms_type[tid2];
      rbmd::Real eps_j = eps[typej];
      rbmd::Real sigma_j = sigma[typej];
      rbmd::Real charge_j = charge[tid2];

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

      // erf
      rbmd::Real dis = SQRT(px12 * px12 + py12 * py12 + pz12 * pz12);
      rbmd::Id index_table_pij = Extract(dis);
      rbmd::Real table_pij = TableGnearValue(erf_table,dis, index_table_pij);

      // compute the force_rcs
      rbmd::Real force_lj_rcs, force_coul_rcs,force_coul_factor_rcs;
      rbmd::Real fcs_ij;
      lj126_rcs(rc, rs, pice_num, px12, py12, pz12, eps_ij, sigma_ij,
                force_lj_rcs);

      rbmd::Real weight = 1.0;
      for (rbmd::Id k = 0; k < special_count[atom_id1]; ++k) {
        rbmd::Id special_id = special_ids[num_components + k];
        if (special_id == atom_id2) {
          weight = special_weights[num_components + k];
          // printf("weight %f\n", weight);
        }
      }
      force_lj_rcs = weight * force_lj_rcs;

      CoulCutForce_rcs_fix(rc,rs, pice_num, alpha, qqr2e,
      charge_i, charge_j, px12, py12, pz12, force_coul_factor_rcs,
    force_coul_rcs);
      force_coul_rcs = force_coul_rcs - (1-weight)*force_coul_factor_rcs;

      fcs_ij = force_lj_rcs + force_coul_rcs;
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

//verlet-list : SpecialLJCutCoul Energy
__global__ void ComputeSpecialLJCutCoulEnergy(
     Box box, ERFTable* erf_table, const rbmd::Real cut_off,
    const rbmd::Id num_atoms, const rbmd::Real alpha, const rbmd::Real qqr2e,
    const rbmd::Id* atoms_type, const rbmd::Id* atoms_id,
    const rbmd::Real* sigma, const rbmd::Real* eps, const rbmd::Id* start_id,
    const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
    const rbmd::Id* special_ids, const rbmd::Real* special_weights,
    const rbmd::Id* special_offset, const rbmd::Id* special_count,
    const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz,  rbmd::Real* flat_virial,
    rbmd::Real* total_evdwl, rbmd::Real* total_ecoul) {
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
      temp_storage_elj;
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
      temp_storage_ecoul;
  rbmd::Real sum_elj = 0;
  rbmd::Real sum_ecoul = 0;
  //virial init
  rbmd::Real sum_virial[6];
  for (int i = 0; i < 6; ++i)
  {
    sum_virial[i] = 0.0;
  }

  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_atoms) {
    rbmd::Id atom_id1 = atoms_id[tid1];
    rbmd::Id num_components = special_offset[atom_id1];
    rbmd::Id typei = atoms_type[tid1];
    rbmd::Real eps_i = eps[typei];
    rbmd::Real sigma_i = sigma[typei];
    rbmd::Real charge_i = charge[tid1];
    rbmd::Real x1 = px[tid1];
    rbmd::Real y1 = py[tid1];
    rbmd::Real z1 = pz[tid1];
    for (int j = start_id[tid1]; j < end_id[tid1]; ++j) {
      rbmd::Id tid2 = id_verletlist[j];
      rbmd::Id atom_id2 = atoms_id[tid2];

      rbmd::Id typej = atoms_type[tid2];
      rbmd::Real eps_j = eps[typej];
      rbmd::Real sigma_j = sigma[typej];
      rbmd::Real charge_j = charge[tid2];

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

      // erf
      rbmd::Real dis = SQRT(px12 * px12 + py12 * py12 + pz12 * pz12);
      rbmd::Id index_table_pij = Extract(dis);
      rbmd::Real table_pij = TableGnearValue(erf_table,dis, index_table_pij);

      rbmd::Real force_lj, force_coul, force_pair,force_coul_factor;
      rbmd::Real energy_lj, energy_coul,energy_coul_factor;

      // lj cut
      lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij, force_lj, energy_lj);

      rbmd::Real weight = 1.0;
      for (rbmd::Id k = 0; k < special_count[atom_id1]; ++k) {
        rbmd::Id special_id = special_ids[num_components + k];
        if (special_id == atom_id2) {
          weight = special_weights[num_components + k];
        }
      }
      force_lj = weight * force_lj;
      energy_lj = weight * energy_lj;

      // Coul cut
      // CoulCutForce_erf(cut_off, alpha, qqr2e, table_pij, charge_i, charge_j,
      //                  px12, py12, pz12, force_coul, energy_coul);
      CoulCutForce_fix(cut_off, alpha, qqr2e, charge_i, charge_j,
        px12, py12, pz12, force_coul_factor,energy_coul_factor,
      force_coul,energy_coul);
      force_coul = force_coul-(1-weight)*force_coul_factor;
      energy_coul = energy_coul- (1-weight)*energy_coul_factor;

      //sum force of special_lj_cut  and coul
      force_pair = force_lj + force_coul;

      //sum energy  of special_lj_cut  and coul
      sum_elj += energy_lj;
      sum_ecoul += energy_coul;


      //rbmd::Real local_virial[6];
      rbmd::Real local_virial_xx,local_virial_yy,local_virial_zz,
  local_virial_xy,local_virial_xz,local_virial_yz;
      ComputeVirial_fix(px12, py12, pz12,force_pair,local_virial_xx,local_virial_yy,
    local_virial_zz,local_virial_xy,local_virial_xz,local_virial_yz);

      //
      sum_virial[0] +=local_virial_xx;
      sum_virial[1] +=local_virial_yy;
      sum_virial[2] +=local_virial_zz;
      sum_virial[3] +=local_virial_xy;
      sum_virial[4] +=local_virial_xz;
      sum_virial[5] +=local_virial_yz;
    }

    //
    for(int i =0;i<6;++i) {
      flat_virial[ tid1 * 6 + i ] = sum_virial[i];
    }
  }
  rbmd::Real block_sum_elj =
      BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage_elj)
          .Sum(sum_elj);
  rbmd::Real block_sum_ecoul =
      BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage_ecoul)
          .Sum(sum_ecoul);

  if (threadIdx.x == 0) {
    atomicAdd(total_evdwl, block_sum_elj);
    atomicAdd(total_ecoul, block_sum_ecoul);
  }
}

 //bond
__global__ void ComputeBondForce(
     Box box, const rbmd::Id num_atoms,const rbmd::Id num_bonds,
     const rbmd::Id* atom_id_to_idx,
    const rbmd::Real* bond_coeffs_k, const rbmd::Real* bond_coeffs_equilibrium,
    const rbmd::Id* bond_type, const rbmd::Id* bondlisti,
    const rbmd::Id* bondlistj, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial,rbmd::Real* energy_bond) {
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
      temp_storage;
  rbmd::Real local_energy_bond = 0;
  // initialize the virial to 0
  rbmd::Real sum_virial[6];
  for (int i = 0; i < 6; i++) {
    sum_virial[i] = 0.0;
  }

  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_bonds) {
    rbmd::Id bondi = bondlisti[tid1];
    rbmd::Id bondj = bondlistj[tid1];
    rbmd::Id bondii = atom_id_to_idx[bondi];
    rbmd::Id bondjj = atom_id_to_idx[bondj];

    rbmd::Id bondtype = bond_type[tid1];
    rbmd::Real k = bond_coeffs_k[bondtype];
    rbmd::Real equilibrium_bond = bond_coeffs_equilibrium[bondtype];

    rbmd::Real x12 = px[bondii] - px[bondjj];
    rbmd::Real y12 = py[bondii] - py[bondjj];
    rbmd::Real z12 = pz[bondii] - pz[bondjj];
    MinImageDistance(box, x12, y12, z12);
    rbmd::Real dis_12 = SQRT(x12 * x12 + y12 * y12 + z12 * z12);
    rbmd::Real dr = dis_12 - equilibrium_bond;
    rbmd::Real rk = k * dr;

    // energy
    local_energy_bond = rk * dr; //double counting??

    rbmd::Real forcebondij;
    if (dis_12 > 0.0) //??
      forcebondij = -2.0 * rk / dis_12;
    else
      forcebondij = 0.0;

    rbmd::Real fx_ij = forcebondij * x12;
    rbmd::Real fy_ij = forcebondij * y12;
    rbmd::Real fz_ij = forcebondij * z12;

    // apply force to each of 2 atoms
    atomicAdd(&fx[bondii], fx_ij);
    atomicAdd(&fy[bondii], fy_ij);
    atomicAdd(&fz[bondii], fz_ij);

    atomicAdd(&fx[bondjj], -fx_ij);
    atomicAdd(&fy[bondjj], -fy_ij);
    atomicAdd(&fz[bondjj], -fz_ij);

    rbmd::Real local_virial[6];
    local_virial[0]  = 0.5 *x12 *x12 *forcebondij; //double counting
    local_virial[1]  = 0.5 *y12 *y12 *forcebondij;
    local_virial[2]  = 0.5 *z12 *z12 *forcebondij;
    local_virial[3]  = 0.5 *x12 *y12 *forcebondij;
    local_virial[4]  = 0.5 *x12 *z12 *forcebondij;
    local_virial[5]  = 0.5 *y12 *z12 *forcebondij;

    // for (int i = 0; i < 6; i++)
    //   sum_virial[i] += local_virial[i];

    // 将每个 bond 的 virial 分量加到相应的原子
    atomicAdd(&flat_virial[0 * num_atoms + bondii], local_virial[0]);
    atomicAdd(&flat_virial[1 * num_atoms + bondii], local_virial[1]);
    atomicAdd(&flat_virial[2 * num_atoms + bondii], local_virial[2]);
    atomicAdd(&flat_virial[3 * num_atoms + bondii], local_virial[3]);
    atomicAdd(&flat_virial[4 * num_atoms + bondii], local_virial[4]);
    atomicAdd(&flat_virial[5 * num_atoms + bondii], local_virial[5]);

    atomicAdd(&flat_virial[0 * num_atoms + bondjj], local_virial[0]);
    atomicAdd(&flat_virial[1 * num_atoms + bondjj], local_virial[1]);
    atomicAdd(&flat_virial[2 * num_atoms + bondjj], local_virial[2]);
    atomicAdd(&flat_virial[3 * num_atoms + bondjj], local_virial[3]);
    atomicAdd(&flat_virial[4 * num_atoms + bondjj], local_virial[4]);
    atomicAdd(&flat_virial[5 * num_atoms + bondjj], local_virial[5]);
    //
    // for(int i =0;i<6;++i) {
    //   flat_virial[  i * num_atoms + tid1 ] = local_virial[i];
    // }
  }
  rbmd::Real block_sum =
      BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage)
          .Sum(local_energy_bond);

  if (threadIdx.x == 0) {
    atomicAdd(energy_bond, block_sum);
  }
}

//angle
  __global__ void ComputeAngleForce(
       Box box, const rbmd::Id num_atoms,const rbmd::Id num_anglels,
       const rbmd::Id* atom_id_to_idx,
      const rbmd::Real* anglel_coeffs_k,
      const rbmd::Real* anglel_coeffs_equilibrium, const rbmd::Id* anglel_type,
      const rbmd::Id* anglelisti, const rbmd::Id* anglelistj,
      const rbmd::Id* anglelistk, const rbmd::Real* px, const rbmd::Real* py,
      const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
      rbmd::Real* flat_virial,rbmd::Real* energy_angle) {
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage;
    rbmd::Real local_energy_angle = 0;

    // initialize the virial to 0
    rbmd::Real sum_virial[6];
    for (int i = 0; i < 6; i++) {
      sum_virial[i] = 0.0;
    }

    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_anglels) {
      rbmd::Id angleli = anglelisti[tid1];
      rbmd::Id anglelj = anglelistj[tid1];
      rbmd::Id anglelk = anglelistk[tid1];
      rbmd::Id anglelii = atom_id_to_idx[angleli];
      rbmd::Id angleljj = atom_id_to_idx[anglelj];
      rbmd::Id anglelkk = atom_id_to_idx[anglelk];

      rbmd::Id angleltype = anglel_type[tid1];
      rbmd::Real k = anglel_coeffs_k[angleltype];
      rbmd::Real equilibrium_angle = anglel_coeffs_equilibrium[angleltype];

      rbmd::Real x12 = px[anglelii] - px[angleljj];  // i j
      rbmd::Real y12 = py[anglelii] - py[angleljj];
      rbmd::Real z12 = pz[anglelii] - pz[angleljj];
      MinImageDistance(box, x12, y12, z12);

      rbmd::Real x23 = px[anglelkk] - px[angleljj];  // k j
      rbmd::Real y23 = py[anglelkk] - py[angleljj];
      rbmd::Real z23 = pz[anglelkk] - pz[angleljj];
      MinImageDistance(box, x23, y23, z23);

      rbmd::Real dis_12_2 = x12 * x12 + y12 * y12 + z12 * z12;
      rbmd::Real dis_12 = SQRT(dis_12_2);

      rbmd::Real dis_23_2 = x23 * x23 + y23 * y23 + z23 * z23;
      rbmd::Real dis_23 = SQRT(dis_23_2);

      rbmd::Real cosangle = x12 * x23 + y12 * y23 + z12 * z23;
      cosangle /= dis_12 * dis_23;

      if (cosangle > 1.0) cosangle = 1.0;
      if (cosangle < -1.0) cosangle = -1.0;
      rbmd::Real s = SQRT(1.0 - cosangle * cosangle);

      if (s < SMALL) s = SMALL;
      s = 1.0 / s;

      rbmd::Real dtheta = ACOS(cosangle) - (equilibrium_angle * M_PI) / 180;
      rbmd::Real tk = k * dtheta;

      // energy
      local_energy_angle = tk * dtheta;

      rbmd::Real a = -2.0 * tk * s;
      rbmd::Real a11 = a * cosangle / dis_12_2;
      rbmd::Real a12 = -a / (dis_12 * dis_23);
      rbmd::Real a22 = a * cosangle / dis_23_2;

      rbmd::Real force_anglei_x, force_anglei_y, force_anglei_z;
      rbmd::Real force_anglek_x, force_anglek_y, force_anglek_z;
      rbmd::Real force_anglej_x, force_anglej_y, force_anglej_z;

      force_anglei_x = a11 * x12 + a12 * x23;
      force_anglei_y = a11 * y12 + a12 * y23;
      force_anglei_z = a11 * z12 + a12 * z23;

      force_anglek_x = a22 * x23 + a12 * x12;
      force_anglek_y = a22 * y23 + a12 * y12;
      force_anglek_z = a22 * z23 + a12 * z12;

      force_anglej_x = -(force_anglei_x + force_anglek_x);
      force_anglej_y = -(force_anglei_y + force_anglek_y);
      force_anglej_z = -(force_anglei_z + force_anglek_z);

      // fx[anglelii] = force_anglei_x;
      // fy[anglelii] = force_anglei_y;
      // fz[anglelii] = force_anglei_z;
      //
      // fx[anglelkk] = force_anglek_x;
      // fy[anglelkk] = force_anglek_y;
      // fz[anglelkk] = force_anglek_z;
      //
      // fx[angleljj] = force_anglej_x;
      // fy[angleljj] = force_anglej_y;
      // fz[angleljj] = force_anglej_z;

      atomicAdd(&fx[anglelii], force_anglei_x);
      atomicAdd(&fy[anglelii], force_anglei_y);
      atomicAdd(&fz[anglelii], force_anglei_z);

      atomicAdd(&fx[anglelkk], force_anglek_x);
      atomicAdd(&fy[anglelkk], force_anglek_y);
      atomicAdd(&fz[anglelkk], force_anglek_z);

      atomicAdd(&fx[angleljj], force_anglej_x);
      atomicAdd(&fy[angleljj], force_anglej_y);
      atomicAdd(&fz[angleljj], force_anglej_z);

      //virial
      rbmd::Real local_virial[6];
      local_virial[0]  =  0.3333333333*(x12 * force_anglei_x + x23 * force_anglek_x);
      local_virial[1]  =  0.3333333333*(y12 * force_anglei_y + y23 * force_anglek_y);
      local_virial[2]  =  0.3333333333*(z12 * force_anglei_z + z23 * force_anglek_z);
      local_virial[3]  =  0.3333333333*(x12 * force_anglei_x + x23 * force_anglek_y);
      local_virial[4]  =  0.3333333333*(x12 * force_anglei_z + x23 * force_anglek_z);
      local_virial[5]  =  0.3333333333*(y12 * force_anglei_z + y23 * force_anglek_z);
      // for (int i = 0; i < 6; i++)
      //   sum_virial[i] += local_virial[i];

      //
      // for(int i =0;i<6;++i) {
      //
      //   flat_virial[ i * num_atoms + tid1 ] = sum_virial[i];
      // }

      // 将每个 angle 的 virial 分量加到相应的原子
      //Column-Major Order : j * M + i
      atomicAdd(&flat_virial[0 * num_atoms + anglelii], local_virial[0]);
      atomicAdd(&flat_virial[1 * num_atoms + anglelii], local_virial[1]);
      atomicAdd(&flat_virial[2 * num_atoms + anglelii], local_virial[2]);
      atomicAdd(&flat_virial[3 * num_atoms + anglelii], local_virial[3]);
      atomicAdd(&flat_virial[4 * num_atoms + anglelii], local_virial[4]);
      atomicAdd(&flat_virial[5 * num_atoms + anglelii], local_virial[5]);

      atomicAdd(&flat_virial[0 * num_atoms + angleljj], local_virial[0]);
      atomicAdd(&flat_virial[1 * num_atoms + angleljj], local_virial[1]);
      atomicAdd(&flat_virial[2 * num_atoms + angleljj], local_virial[2]);
      atomicAdd(&flat_virial[3 * num_atoms + angleljj], local_virial[3]);
      atomicAdd(&flat_virial[4 * num_atoms + angleljj], local_virial[4]);
      atomicAdd(&flat_virial[5 * num_atoms + angleljj], local_virial[5]);

      atomicAdd(&flat_virial[0 * num_atoms + anglelkk], local_virial[0]);
      atomicAdd(&flat_virial[1 * num_atoms + anglelkk], local_virial[1]);
      atomicAdd(&flat_virial[2 * num_atoms + anglelkk], local_virial[2]);
      atomicAdd(&flat_virial[3 * num_atoms + anglelkk], local_virial[3]);
      atomicAdd(&flat_virial[4 * num_atoms + anglelkk], local_virial[4]);
      atomicAdd(&flat_virial[5 * num_atoms + anglelkk], local_virial[5]);

    }

    rbmd::Real block_sum =
        BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage)
            .Sum(local_energy_angle);

    if (threadIdx.x == 0) {
      atomicAdd(energy_angle, block_sum);
    }
  }

  //Dihedral
  __global__ void ComputeDihedralForce(
       Box box, const rbmd::Id num_atoms,const rbmd::Id num_dihedrals,
       const rbmd::Id* atom_id_to_idx,const rbmd::Real* dihedral_coeffs_k,
       const rbmd::Id* dihedral_coeffs_sign,
      const rbmd::Id* dihedral_coeffs_multiplicity, const rbmd::Id* dihedral_type,
      const rbmd::Id* dihedrallisti, const rbmd::Id* dihedrallistj,
      const rbmd::Id* dihedrallistk, const rbmd::Id* dihedrallistw,
      const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
      rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,rbmd::Real* flat_virial,
      rbmd::Real* energy_dihedral) {
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage;
    rbmd::Real local_energy_dihedral = 0;

    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_dihedrals) {
      rbmd::Id dihedrali = dihedrallisti[tid1];
      rbmd::Id dihedralj = dihedrallistj[tid1];
      rbmd::Id dihedralk = dihedrallistk[tid1];
      rbmd::Id dihedralw = dihedrallistw[tid1];

      rbmd::Id dihedralii = atom_id_to_idx[dihedrali];
      rbmd::Id dihedraljj = atom_id_to_idx[dihedralj];
      rbmd::Id dihedralkk = atom_id_to_idx[dihedralk];
      rbmd::Id dihedralww = atom_id_to_idx[dihedralw];
      rbmd::Real cos_shift, sin_shift;

      rbmd::Id dihedraltype = dihedral_type[tid1];
      if (dihedral_coeffs_sign[dihedraltype] == 1) {
        cos_shift = 1.0;
        sin_shift = 0.0;
      } else {
        cos_shift = -1.0;
        sin_shift = 0.0;
      }

      rbmd::Real k = dihedral_coeffs_k[dihedraltype];
      rbmd::Real x12 = px[dihedralii] - px[dihedraljj];  // i j =vb1
      rbmd::Real y12 = py[dihedralii] - py[dihedraljj];
      rbmd::Real z12 = pz[dihedralii] - pz[dihedraljj];
      MinImageDistance(box, x12, y12, z12);

      rbmd::Real x23 = px[dihedralkk] - px[dihedraljj];  //  k j=vb2
      rbmd::Real y23 = py[dihedralkk] - py[dihedraljj];
      rbmd::Real z23 = pz[dihedralkk] - pz[dihedraljj];
      MinImageDistance(box, x23, y23, z23);

      rbmd::Real x23m = -x23;  // =vb2m
      rbmd::Real y23m = -y23;
      rbmd::Real z23m = -z23;

      rbmd::Real x34 = px[dihedralww] - px[dihedralkk];  // w k   =vb3
      rbmd::Real y34 = py[dihedralww] - py[dihedralkk];
      rbmd::Real z34 = pz[dihedralww] - pz[dihedralkk];
      MinImageDistance(box, x34, y34, z34);
      // c,s calculation

      rbmd::Real ax = y12 * z23m - z12 * y23m;
      rbmd::Real ay = z12 * x23m - x12 * z23m;
      rbmd::Real az = x12 * y23m - y12 * x23m;
      rbmd::Real bx = y34 * z23m - z34 * y23m;
      rbmd::Real by = z34 * x23m - x34 * z23m;
      rbmd::Real bz = x34 * y23m - z34 * x23m;
      rbmd::Real rasq = ax * ax + ay * ay + az * az;
      rbmd::Real rbsq = bx * bx + by * by + bz * bz;
      rbmd::Real rgsq = x23m * x23m + y23m * y23m + z23m * z23m;
      rbmd::Real rg = SQRT(rgsq);

      rbmd::Real rginv, ra2inv, rb2inv;
      rginv = ra2inv = rb2inv = 0.0;
      if (rg > 0) rginv = 1.0 / rg;

      if (rasq > 0) ra2inv = 1.0 / rasq;

      if (rbsq > 0) rb2inv = 1.0 / rbsq;

      rbmd::Real rabinv = SQRT(ra2inv * rb2inv);

      rbmd::Real c = (ax * bx + ay * by + az * bz) * rabinv;
      rbmd::Real s = rg * rabinv * (ax * x34 + ay * y34 + az * z34);

      if (c > 1.0) c = 1.0;
      if (c < -1.0) c = -1.0;

      rbmd::Id m = dihedral_coeffs_multiplicity[dihedraltype];
      rbmd::Real p = 1.0;
      rbmd::Real ddf1, df1;
      ddf1 = df1 = 0.0;
      for (rbmd::Id i = 0; i < m; i++) {
        ddf1 = p * c - df1 * s;
        df1 = p * s + df1 * c;
        p = ddf1;
      }

      p = p * cos_shift + df1 * sin_shift;
      df1 = df1 * cos_shift - ddf1 * sin_shift;
      df1 *= -m;
      p += 1.0;

      if (m == 0) {
        p = 1.0 + cos_shift;
        // p = 1.0 + cos_shift[type];
        df1 = 0.0;
      }

      // energy
      local_energy_dihedral = k * p;

      rbmd::Real fg = x12 * x23m + y12 * y23m + z12 * z23m;
      rbmd::Real hg = x34 * x23m + y34 * y23m + z34 * z23m;

      rbmd::Real fga = fg * ra2inv * rginv;
      rbmd::Real hgb = hg * rb2inv * rginv;
      rbmd::Real gaa = -ra2inv * rg;
      rbmd::Real gbb = rb2inv * rg;

      rbmd::Real dtfx = gaa * ax;
      rbmd::Real dtfy = gaa * ay;
      rbmd::Real dtfz = gaa * az;
      rbmd::Real dtgx = fga * ax - hgb * bx;
      rbmd::Real dtgy = fga * ay - hgb * by;
      rbmd::Real dtgz = fga * az - hgb * bz;
      rbmd::Real dthx = gbb * bx;
      rbmd::Real dthy = gbb * by;
      rbmd::Real dthz = gbb * bz;

      rbmd::Real df = -k * df1;
      // df = -k[type] * df1;
      rbmd::Real sx2 = df * dtgx;
      rbmd::Real sy2 = df * dtgy;
      rbmd::Real sz2 = df * dtgz;

      // force
      rbmd::Real force_dihedrali_x, force_dihedrali_y, force_dihedrali_z;
      rbmd::Real force_dihedralj_x, force_dihedralj_y, force_dihedralj_z;
      rbmd::Real force_dihedralk_x, force_dihedralk_y, force_dihedralk_z;
      rbmd::Real force_dihedralw_x, force_dihedralw_y, force_dihedralw_z;

      force_dihedrali_x = df * dtfx;
      force_dihedrali_y = df * dtfy;
      force_dihedrali_z = df * dtfz;

      force_dihedralj_x = sx2 - force_dihedrali_x;
      force_dihedralj_y = sy2 - force_dihedrali_y;
      force_dihedralj_z = sz2 - force_dihedrali_z;

      force_dihedralw_x = df * dthx;
      force_dihedralw_y = df * dthy;
      force_dihedralw_z = df * dthz;

      force_dihedralk_x = -sx2 - force_dihedralw_x;
      force_dihedralk_y = -sy2 - force_dihedralw_y;
      force_dihedralk_z = -sz2 - force_dihedralw_z;

      atomicAdd(&fx[dihedralii], force_dihedrali_x);
      atomicAdd(&fy[dihedralii], force_dihedrali_y);
      atomicAdd(&fz[dihedralii], force_dihedrali_z);

      atomicAdd(&fx[dihedraljj], force_dihedralj_x);
      atomicAdd(&fy[dihedraljj], force_dihedralj_y);
      atomicAdd(&fz[dihedraljj], force_dihedralj_z);

      atomicAdd(&fx[dihedralww], force_dihedralw_x);
      atomicAdd(&fy[dihedralww], force_dihedralw_y);
      atomicAdd(&fz[dihedralww], force_dihedralw_z);

      atomicAdd(&fx[dihedralkk], force_dihedralk_x);
      atomicAdd(&fy[dihedralkk], force_dihedralk_y);
      atomicAdd(&fz[dihedralkk], force_dihedralk_z);
      // fx[dihedralii] += force_dihedrali_x;  // i atom force
      // fy[dihedralii] += force_dihedrali_y;
      // fz[dihedralii] += force_dihedrali_z;
      //
      // fx[dihedraljj] += force_dihedralj_x;  // j atom force
      // fy[dihedraljj] += force_dihedralj_y;
      // fz[dihedraljj] += force_dihedralj_z;
      //
      // fx[dihedralkk] += force_dihedralk_x;  // k atom force
      // fy[dihedralkk] += force_dihedralk_y;
      // fz[dihedralkk] += force_dihedralk_z;
      //
      // fx[dihedralww] += force_dihedralw_x;  // w atom force
      // fy[dihedralww] += force_dihedralw_y;
      // fz[dihedralww] += force_dihedralw_z;


      rbmd::Real local_virial[6];
      local_virial[0] = 0.25* (x12 * force_dihedrali_x + x23 * force_dihedralk_x +
          (x34 + x23) * force_dihedralw_x);

      local_virial[1] = 0.25* (y12 * force_dihedrali_y + y23 * force_dihedralk_y +
         (y34 + y23) * force_dihedralw_y);

      local_virial[2] = 0.25* (z12 * force_dihedrali_z + z23 * force_dihedralk_z +
         (z34 + z23) * force_dihedralw_z);

      local_virial[3] = 0.25* (x12 * force_dihedrali_y + x23 * force_dihedralk_y +
         (x34 + x23) * force_dihedralw_y);

      local_virial[4] = 0.25* (x12* force_dihedrali_z + x23 * force_dihedralk_z +
         (x34 + x23) * force_dihedralw_z);

      local_virial[5] =  0.25* (y12 * force_dihedrali_z + y23 * force_dihedralk_z +
         (y34 + y23) * force_dihedralw_z);

      // // 将每个 dihedral 的 virial 分量加到相应的原子
      atomicAdd(&flat_virial[0 * num_atoms + dihedralii], local_virial[0]);
      atomicAdd(&flat_virial[1 * num_atoms + dihedralii], local_virial[1]);
      atomicAdd(&flat_virial[2 * num_atoms + dihedralii], local_virial[2]);
      atomicAdd(&flat_virial[3 * num_atoms + dihedralii], local_virial[3]);
      atomicAdd(&flat_virial[4 * num_atoms + dihedralii], local_virial[4]);
      atomicAdd(&flat_virial[5 * num_atoms + dihedralii], local_virial[5]);

      atomicAdd(&flat_virial[0 * num_atoms + dihedraljj], local_virial[0]);
      atomicAdd(&flat_virial[1 * num_atoms + dihedraljj], local_virial[1]);
      atomicAdd(&flat_virial[2 * num_atoms + dihedraljj], local_virial[2]);
      atomicAdd(&flat_virial[3 * num_atoms + dihedraljj], local_virial[3]);
      atomicAdd(&flat_virial[4 * num_atoms + dihedraljj], local_virial[4]);
      atomicAdd(&flat_virial[5 * num_atoms + dihedraljj], local_virial[5]);

      atomicAdd(&flat_virial[0 * num_atoms + dihedralkk], local_virial[0]);
      atomicAdd(&flat_virial[1 * num_atoms + dihedralkk], local_virial[1]);
      atomicAdd(&flat_virial[2 * num_atoms + dihedralkk], local_virial[2]);
      atomicAdd(&flat_virial[3 * num_atoms + dihedralkk], local_virial[3]);
      atomicAdd(&flat_virial[4 * num_atoms + dihedralkk], local_virial[4]);
      atomicAdd(&flat_virial[5 * num_atoms + dihedralkk], local_virial[5]);

      atomicAdd(&flat_virial[0 * num_atoms + dihedralww], local_virial[0]);
      atomicAdd(&flat_virial[1 * num_atoms + dihedralww], local_virial[1]);
      atomicAdd(&flat_virial[2 * num_atoms + dihedralww], local_virial[2]);
      atomicAdd(&flat_virial[3 * num_atoms + dihedralww], local_virial[3]);
      atomicAdd(&flat_virial[4 * num_atoms + dihedralww], local_virial[4]);
      atomicAdd(&flat_virial[5 * num_atoms + dihedralww], local_virial[5]);
      // //
      // for(int i =0;i<6;++i) {
      //   flat_virial[ tid1 * 6 + i ] = local_virial[i];
      // }

    }
    rbmd::Real block_sum =
        BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage)
            .Sum(local_energy_dihedral);

    if (threadIdx.x == 0) {
      atomicAdd(energy_dihedral, block_sum);
    }
  }

  //Imprope
  __global__ void ComputeImproperForce(
    Box box,const rbmd::Id num_impropers,const rbmd::Id* atom_id_to_idx,
    const rbmd::Real* improper_coeffs_k,const rbmd::Real* improper_coeffs_chi,
    const rbmd::Id* improper_type,const rbmd::Id* improperlisti,
    const rbmd::Id* improperlistj,const rbmd::Id* improperlistk,
    const rbmd::Id* improperlistw,const rbmd::Real* px,const rbmd::Real* py,
    const rbmd::Real* pz,rbmd::Real* fx,rbmd::Real* fy,rbmd::Real* fz,
    rbmd::Real* flat_virial,rbmd::Real* energy_improper) {
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage;
    rbmd::Real local_energy_improper = 0;

    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_impropers) {
      rbmd::Id improperi = improperlisti[tid1];
      rbmd::Id improperj = improperlistj[tid1];
      rbmd::Id improperk = improperlistk[tid1];
      rbmd::Id improperw = improperlistw[tid1];

      rbmd::Id improperii = atom_id_to_idx[improperi];
      rbmd::Id improperjj = atom_id_to_idx[improperj];
      rbmd::Id improperkk = atom_id_to_idx[improperk];
      rbmd::Id improperww = atom_id_to_idx[improperw];
      rbmd::Id impropertype = improper_type[tid1];

      rbmd::Real k = improper_coeffs_k[impropertype];
      rbmd::Real chi = improper_coeffs_chi[impropertype];
      rbmd::Real x12 = px[improperii] - px[improperjj];  // i j =vb1
      rbmd::Real y12 = py[improperii] - py[improperjj];
      rbmd::Real z12 = pz[improperii] - pz[improperjj];
      MinImageDistance(box, x12, y12, z12);

      rbmd::Real x23 = px[improperkk] - px[improperjj];  //  k j=vb2
      rbmd::Real y23 = py[improperkk] - py[improperjj];
      rbmd::Real z23 = pz[improperkk] - pz[improperjj];
      MinImageDistance(box, x23, y23, z23);

      rbmd::Real x34 = px[improperww] - px[improperkk];  // w k   =vb3
      rbmd::Real y34 = py[improperww] - py[improperkk];
      rbmd::Real z34 = pz[improperww] - pz[improperkk];
      MinImageDistance(box, x34, y34, z34);

      rbmd::Real ss1 = 1.0 / (x12 * x12 + y12 * y12 + z12 * z12);
      rbmd::Real ss2 = 1.0 / (x23 * x23 + y23 * y23 + z23 * z23);
      rbmd::Real ss3 = 1.0 / (x34 * x34 + y34 * y34 + z34 * z34);

      rbmd::Real r1 = SQRT(ss1);
      rbmd::Real r2 = SQRT(ss2);
      rbmd::Real r3 = SQRT(ss3);

      // sin and cos of angle

      rbmd::Real c0 = (x12 * x34 + y12 * y34 + z12 * z34) * r1 * r3;
      rbmd::Real c1 = (x12 * x23 + y12 * y23 + z12 * z23) * r1 * r2;
      rbmd::Real c2 = -(x34 * x23 + y34 * y23 + z34 * z23) * r3 * r2;

      rbmd::Real  s1 = 1.0 - c1 * c1;
      if (s1 < SMALL) s1 = SMALL;
      s1 = 1.0 / s1;

       rbmd::Real s2 = 1.0 - c2 * c2;
      if (s2 < SMALL) s2 = SMALL;
      s2 = 1.0 / s2;

      rbmd::Real s12 = SQRT(s1 * s2);
      rbmd::Real c = (c1 * c2 + c0) * s12;

      // error check

      if (c > 1.0) c = 1.0;
      if (c < -1.0) c = -1.0;

      rbmd::Real s = SQRT(1.0 - c * c);
      if (s < SMALL) s = SMALL;


      // force & energy
      rbmd::Real domega = ACOS(c) - chi;
      rbmd::Real a = k * domega;
      local_energy_improper = a * domega;

      a = -a * 2.0 / s;
      c = c * a;
      s12 = s12 * a;
      rbmd::Real a11 = c * ss1 * s1;
      rbmd::Real a22 = -ss2 * (2.0 * c0 * s12 - c * (s1 + s2));
      rbmd::Real a33 = c * ss3 * s2;
      rbmd::Real a12 = -r1 * r2 * (c1 * c * s1 + c2 * s12);
      rbmd::Real a13 = -r1 * r3 * s12;
      rbmd::Real a23 = r2 * r3 * (c2 * c * s2 + c1 * s12);

      rbmd::Real sx2 = a22 * x23 + a23 * x34 + a12 * x12;
      rbmd::Real sy2 = a22 * y23 + a23 * y34 + a12 * y12;
      rbmd::Real sz2 = a22 * z23 + a23 * z34 + a12 * z12;

      // force
      rbmd::Real force_improperi_x, force_improperi_y, force_improperi_z;
      rbmd::Real force_improperj_x, force_improperj_y, force_improperj_z;
      rbmd::Real force_improperk_x, force_improperk_y, force_improperk_z;
      rbmd::Real force_improperw_x, force_improperw_y, force_improperw_z;

      force_improperi_x = a12 * x23 + a13 * x34 + a11 * x12;
      force_improperi_y = a12 * y23 + a13 * y34 + a11 * y12;
      force_improperi_z = a12 * z23 + a13 * z34 + a11 * z12;

      force_improperj_x = -sx2 - force_improperi_x;
      force_improperj_y = -sy2 - force_improperi_y;
      force_improperj_z = -sz2 - force_improperi_z;

      force_improperw_x = a23 * x23 + a33 * x34 + a13 * x12;
      force_improperw_y = a23 * y23 + a33 * y34 + a13 * y12;
      force_improperw_z = a23 * z23 + a33 * z34 + a13 * z12;

      force_improperk_x = sx2 - force_improperw_x;
      force_improperk_y = sy2 - force_improperw_y;
      force_improperk_z = sz2 - force_improperw_z;


      atomicAdd(&fx[improperii], force_improperi_x);
      atomicAdd(&fy[improperii], force_improperi_y);
      atomicAdd(&fz[improperii], force_improperi_z);

      atomicAdd(&fx[improperjj], force_improperj_x);
      atomicAdd(&fy[improperjj], force_improperj_y);
      atomicAdd(&fz[improperjj], force_improperj_z);

      atomicAdd(&fx[improperww], force_improperw_x);
      atomicAdd(&fy[improperww], force_improperw_y);
      atomicAdd(&fz[improperww], force_improperw_z);

      atomicAdd(&fx[improperkk], force_improperk_x);
      atomicAdd(&fy[improperkk], force_improperk_y);
      atomicAdd(&fz[improperkk], force_improperk_z);

      rbmd::Real local_virial[6];
      local_virial[0] = (x12 * force_improperi_x + x23 * force_improperk_x +
          (x34 + x23) * force_improperw_x);

      local_virial[1] = (y12 * force_improperi_y+ y23 * force_improperk_y +
         (y34 + y23) * force_improperw_y);

      local_virial[2] = (z12 * force_improperi_z + z23 * force_improperk_z +
         (z34 + z23) * force_improperw_z);

      local_virial[3] = (x12 * force_improperi_y + x23 * force_improperk_y +
         (x34 + x23) * force_improperw_y);

      local_virial[4] = (x12* force_improperi_z + x23 * force_improperk_z +
         (x34 + x23) * force_improperw_z);

      local_virial[5] =  (y12 * force_improperi_z + y23 * force_improperk_z +
         (y34 + y23) * force_improperw_z);

      //
      for(int i =0;i<6;++i) {
        flat_virial[ tid1 * 6 + i ] = local_virial[i];
      }

    }
    rbmd::Real block_sum =
        BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage)
            .Sum(local_energy_improper);

    if (threadIdx.x == 0) {
      atomicAdd(energy_improper, block_sum);
    }
  }

////////////////////////////////////////
//verlet-list:  force  of special lJ_cut and coul_cut
  void SpecialLJCutCoulForceOp<device::DEVICE_GPU>::operator()(
       Box box, ERFTable* erf_table, const rbmd::Real cut_off,
      const rbmd::Id num_atoms, const rbmd::Real alpha, const rbmd::Real qqr2e,
      const rbmd::Id* atoms_type, const rbmd::Id* atoms_id,
      const rbmd::Real* sigma, const rbmd::Real* eps, const rbmd::Id* start_id,
      const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
      const rbmd::Id* special_ids, const rbmd::Real* special_weights,
      const rbmd::Id* special_offset, const rbmd::Id* special_count,
      const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py,
      const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
      rbmd::Real* flat_virial,rbmd::Real* total_evdwl, rbmd::Real* total_ecoul) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(
        ComputeSpecialLJCutCoulForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
            box, erf_table, cut_off, num_atoms, alpha, qqr2e, atoms_type,
            atoms_id, sigma, eps, start_id, end_id, id_verletlist, special_ids,
            special_weights, special_offset, special_count, charge, px, py, pz,
            fx, fy, fz, flat_virial,total_evdwl, total_ecoul));
  }

//RBL:  force  of special lJ_cut and coul_cut
  void SpecialLJCutCoulRBLForceOp<device::DEVICE_GPU>::operator()(
       Box box, ERFTable* erf_table, const rbmd::Real rs, const rbmd::Real rc,
      const rbmd::Id num_atoms, const rbmd::Id neighbor_sample_num,
      const rbmd::Id pice_num, const rbmd::Real alpha, const rbmd::Real qqr2e,
      const rbmd::Id* atoms_type, const rbmd::Id* atoms_id,
      const rbmd::Real* sigma, const rbmd::Real* eps, const rbmd::Id* start_id,
      const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
      const rbmd::Id* id_random_neighbor, const rbmd::Id* random_neighbor_num,
      const rbmd::Id* special_ids, const rbmd::Real* special_weights,
      const rbmd::Id* special_offset, const rbmd::Id* special_count,
      const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py,
      const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(
        ComputeSpecialLJCutCoulRBLForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
            box, erf_table, rs, rc, num_atoms, neighbor_sample_num, pice_num,
            alpha, qqr2e, atoms_type, atoms_id, sigma, eps, start_id, end_id,
            id_verletlist, id_random_neighbor, random_neighbor_num, special_ids,
            special_weights, special_offset, special_count, charge, px, py, pz,
            fx, fy, fz));
  }

  //verlet-list:  Energy of special lJ_cut and coul_cut
  void SpeciaLJCutCoulEnergyOp<device::DEVICE_GPU>::operator()(
     Box box, ERFTable* erf_table, const rbmd::Real cut_off,
    const rbmd::Id num_atoms, const rbmd::Real alpha, const rbmd::Real qqr2e,
    const rbmd::Id* atoms_type, const rbmd::Id* atoms_id,
    const rbmd::Real* sigma, const rbmd::Real* eps, const rbmd::Id* start_id,
    const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
    const rbmd::Id* special_ids, const rbmd::Real* special_weights,
    const rbmd::Id* special_offset, const rbmd::Id* special_count,
    const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz, rbmd::Real* flat_virial,
    rbmd::Real* total_evdwl, rbmd::Real* total_ecoul) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(
        ComputeSpecialLJCutCoulEnergy<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
            box, erf_table, cut_off, num_atoms, alpha, qqr2e, atoms_type,
            atoms_id, sigma, eps, start_id, end_id, id_verletlist, special_ids,
            special_weights, special_offset, special_count, charge, px, py, pz,
            flat_virial,total_evdwl, total_ecoul));
  }

  // force of bond
  void ComputeBondForceOp<device::DEVICE_GPU>::operator()(
       Box box,const rbmd::Id num_atoms,const rbmd::Id num_bonds,
       const rbmd::Id* atom_id_to_idx,const rbmd::Real* bond_coeffs_k,
       const rbmd::Real* bond_coeffs_equilibrium,const rbmd::Id* bond_type,
       const rbmd::Id* bondlisti,const rbmd::Id* bondlistj,
       const rbmd::Real* px, const rbmd::Real* py,const rbmd::Real* pz,
       rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
       rbmd::Real* flat_virial, rbmd::Real* energy_bond) {
    unsigned int blocks_per_grid = (num_bonds + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeBondForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, num_atoms,num_bonds, atom_id_to_idx, bond_coeffs_k, bond_coeffs_equilibrium,
        bond_type, bondlisti, bondlistj, px, py, pz, fx, fy, fz,flat_virial,
        energy_bond));
  }

// force of angle
  void ComputeAngleForceOp<device::DEVICE_GPU>::operator()(
       Box box, const rbmd::Id num_atoms,const rbmd::Id num_anglels,
       const rbmd::Id* _atom_id_to_idx,
      const rbmd::Real* anglel_coeffs_k,
      const rbmd::Real* anglel_coeffs_equilibrium, const rbmd::Id* anglel_type,
      const rbmd::Id* anglelisti, const rbmd::Id* anglelistj,
      const rbmd::Id* anglelistk, const rbmd::Real* px, const rbmd::Real* py,
      const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
      rbmd::Real* flat_virial,rbmd::Real* energy_angle) {
    unsigned int blocks_per_grid = (num_anglels + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeAngleForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, num_atoms,num_anglels, _atom_id_to_idx, anglel_coeffs_k,
        anglel_coeffs_equilibrium, anglel_type, anglelisti, anglelistj,
        anglelistk, px, py, pz, fx, fy, fz, flat_virial,energy_angle));
  }


// force of dihedral
  void ComputeDihedralForceOp<device::DEVICE_GPU>::operator()(
       Box box,const rbmd::Id num_atoms,const rbmd::Id num_dihedrals,
       const rbmd::Id* atom_id_to_idx,const rbmd::Real* dihedral_coeffs_k,
       const rbmd::Id* dihedral_coeffs_sign,
       const rbmd::Id* dihedral_coeffs_multiplicity, const rbmd::Id* dihedral_type,
       const rbmd::Id* dihedrallisti, const rbmd::Id* dihedrallistj,
      const rbmd::Id* dihedrallistk, const rbmd::Id* dihedrallistw,
      const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
      rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,rbmd::Real* flat_virial,
      rbmd::Real* energy_dihedral) {
    unsigned int blocks_per_grid = (num_dihedrals + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeDihedralForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, num_atoms,num_dihedrals, atom_id_to_idx, dihedral_coeffs_k,
        dihedral_coeffs_sign, dihedral_coeffs_multiplicity, dihedral_type,
        dihedrallisti, dihedrallistj, dihedrallistk, dihedrallistw, px, py, pz,
        fx, fy, fz, flat_virial,energy_dihedral));
  }

  // force of improper
  void ComputeImproperForceOp<device::DEVICE_GPU>::operator()(
    Box box,const rbmd::Id num_impropers,const rbmd::Id* atom_id_to_idx,
    const rbmd::Real* improper_coeffs_k,const rbmd::Real* improper_coeffs_chi,
    const rbmd::Id* improper_type,const rbmd::Id* improperlisti,
    const rbmd::Id* improperlistj,const rbmd::Id* improperlistk,
    const rbmd::Id* improperlistw,const rbmd::Real* px,const rbmd::Real* py,
    const rbmd::Real* pz,rbmd::Real* fx,rbmd::Real* fy,rbmd::Real* fz,
    rbmd::Real* flat_virial,rbmd::Real* energy_improper) {
      unsigned int blocks_per_grid = (num_impropers + BLOCK_SIZE - 1) / BLOCK_SIZE;

      CHECK_KERNEL(ComputeImproperForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
          box, num_impropers, atom_id_to_idx, improper_coeffs_k,
          improper_coeffs_chi, improper_type,improperlisti, improperlistj,
          improperlistk, improperlistw, px, py, pz,fx, fy, fz,
          flat_virial,energy_improper));
  }

}

