#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"
#include "../common/erf_table.h"

namespace op {

    template <typename DEVICE>
    struct SpecialLJCutCoulForceOp
    {
      void operator()( Box box,ERFTable* erf_table,
            const rbmd::Real cut_off,
            const rbmd::Id num_atoms,
            const rbmd::Real alpha,
            const rbmd::Real qqr2e,
            const rbmd::Id* atoms_type,
            const rbmd::Id* atoms_id,
            const rbmd::Real* sigma,
            const rbmd::Real* eps,
            const rbmd::Id* start_id,
            const rbmd::Id* end_id,
            const rbmd::Id* id_verletlist,
            const rbmd::Id*  special_ids,
            const rbmd::Real*  special_weights,
            const rbmd::Id*  special_offset,
            const rbmd::Id*  special_count,
            const rbmd::Real* charge,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* fx,
            rbmd::Real* fy,
            rbmd::Real* fz,
            rbmd::Real* flat_virial,
            rbmd::Real* total_evdwl,
            rbmd::Real* total_ecoul);
    };

    template <typename DEVICE>
    struct SpecialLJCutCoulRBLForceOp
    {
      void operator()( Box box,ERFTable* erf_table,
        const rbmd::Real rs,
        const rbmd::Real rc,
        const rbmd::Id num_atoms,
        const rbmd::Id neighbor_sample_num,
        const rbmd::Id pice_num,
        const rbmd::Real alpha,
        const rbmd::Real qqr2e,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Real* sigma,
        const rbmd::Real* eps,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const rbmd::Id* id_random_neighbor,
        const rbmd::Id* random_neighbor_num,
        const rbmd::Id*  special_ids,
        const rbmd::Real*  special_weights,
        const rbmd::Id*  special_offset,
        const rbmd::Id*  special_count,
        const rbmd::Real* charge,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* fx,
        rbmd::Real* fy,
        rbmd::Real* fz);
    };

    template <typename DEVICE>
    struct SpeciaLJCutCoulEnergyOp
    {
      void operator()( Box box,ERFTable* erf_table,
           const rbmd::Real cut_off,
           const rbmd::Id num_atoms,
           const rbmd::Real alpha,
           const rbmd::Real qqr2e,
           const rbmd::Id* atoms_type,
           const rbmd::Id* atoms_id,
           const rbmd::Real* sigma,
           const rbmd::Real* eps,
           const rbmd::Id* start_id,
           const rbmd::Id* end_id,
           const rbmd::Id* id_verletlist,
           const rbmd::Id*  special_ids,
           const rbmd::Real*  special_weights,
           const rbmd::Id*  special_offset,
           const rbmd::Id*  special_count,
           const rbmd::Real* charge,
           const rbmd::Real* px,
           const rbmd::Real* py,
           const rbmd::Real* pz,
           rbmd::Real* flat_virial,
           rbmd::Real* total_evdwl,
           rbmd::Real* total_ecoul);
    };

    template <typename DEVICE>
    struct  ComputeSpecialCoulForceOp
    {
      void  operator()(
         Box box,
        const rbmd::Id num_atoms,
        const rbmd::Real qqr2e,
        const rbmd::Id* atoms_id,
        const rbmd::Id*  atoms_vec,
        const rbmd::Id*  atoms_offset,
        const rbmd::Id*  atom_count,
        const rbmd::Id*  special_ids,
        const rbmd::Real*  special_weights,
        const rbmd::Id*  special_offset,
        const rbmd::Id*  special_count,
        const rbmd::Real* charge,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* fx,
        rbmd::Real* fy,
        rbmd::Real* fz,
        rbmd::Real* flat_virial,
        rbmd::Real* total_especial_coul);
    };

     template <typename DEVICE>
     struct  ComputeBondForceOp
     {
       void  operator()(
          Box box,
          const rbmd::Id num_atoms,
         const rbmd::Id num_bonds,
         const rbmd::Id* atom_id_to_idx,
         const rbmd::Real* bond_coeffs_k,
         const rbmd::Real* bond_coeffs_equilibrium,
         const rbmd::Id* bond_type,
         const rbmd::Id* bondlisti,
         const rbmd::Id* bondlistj,
         const rbmd::Real* px,
         const rbmd::Real* py,
         const rbmd::Real* pz,
         rbmd::Real* fx,
         rbmd::Real* fy,
         rbmd::Real* fz,
         rbmd::Real* flat_virial,
         rbmd::Real* global_virial,
         rbmd::Real* energy_bond);
     };

     template <typename DEVICE>
     struct ComputeAngleForceOp
     {
       void operator()(
        Box box,
        const rbmd::Id num_atoms,
       const rbmd::Id num_anglels,
       const rbmd::Id* atom_id_to_idx,
       const rbmd::Real* anglel_coeffs_k,
       const rbmd::Real* anglel_coeffs_equilibrium,
       const rbmd::Id* anglel_type,
       const rbmd::Id* anglelisti,
       const rbmd::Id* anglelistj,
       const rbmd::Id* anglelistk,
       const rbmd::Real* px,
       const rbmd::Real* py,
       const rbmd::Real* pz,
       rbmd::Real* fx,
       rbmd::Real* fy,
       rbmd::Real* fz,
       rbmd::Real* flat_virial,
       rbmd::Real* global_virial,
       rbmd::Real* energy_angle);
     };

     template <typename DEVICE>
     struct ComputeDihedralForceOp
     {
       void operator()(
        Box box,
        const rbmd::Id num_atoms,
       const rbmd::Id num_dihedrals,
       const rbmd::Id* atom_id_to_idx,
       const rbmd::Real* dihedral_coeffs_k,
       const rbmd::Id* dihedral_coeffs_sign ,
       const rbmd::Id* dihedral_coeffs_multiplicity ,
       const rbmd::Id* dihedral_type,
       const rbmd::Id* dihedrallisti,
       const rbmd::Id* dihedrallistj,
       const rbmd::Id* dihedrallistk,
       const rbmd::Id* dihedrallistw,
       const rbmd::Real* px,
       const rbmd::Real* py,
       const rbmd::Real* pz,
       rbmd::Real* fx,
       rbmd::Real* fy,
       rbmd::Real* fz,
       rbmd::Real* flat_virial,
       rbmd::Real* global_virial,
       rbmd::Real* energy_dihedral);
     };

  template <typename DEVICE>
  struct ComputeDihedralOPLSForceOp
  {
    void operator()(
     Box box,
     const rbmd::Id num_atoms,
    const rbmd::Id num_dihedrals,
    const rbmd::Id* atom_id_to_idx,
    const rbmd::Real* dihedral_coeffs_k1,
    const rbmd::Real* dihedral_coeffs_k2,
    const rbmd::Real* dihedral_coeffs_k3,
    const rbmd::Real* dihedral_coeffs_k4,
    const rbmd::Id* dihedral_type,
    const rbmd::Id* dihedrallisti,
    const rbmd::Id* dihedrallistj,
    const rbmd::Id* dihedrallistk,
    const rbmd::Id* dihedrallistw,
    const rbmd::Real* px,
    const rbmd::Real* py,
    const rbmd::Real* pz,
    rbmd::Real* fx,
    rbmd::Real* fy,
    rbmd::Real* fz,
    rbmd::Real* flat_virial,
    rbmd::Real* global_virial,
    rbmd::Real* energy_dihedral);
  };

    template <typename DEVICE>
    struct ComputeImproperHarmonicForceOp
    {
      void operator()(
         Box box,
         const rbmd::Id num_atoms,
        const rbmd::Id num_impropers,
        const rbmd::Id* atom_id_to_idx,
        const rbmd::Real* improper_coeffs_k,
        const rbmd::Real* improper_coeffs_chi,
        const rbmd::Id* improper_type,
        const rbmd::Id* improperlisti,
        const rbmd::Id* improperlistj,
        const rbmd::Id* improperlistk,
        const rbmd::Id* improperlistw,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* fx,
        rbmd::Real* fy,
        rbmd::Real* fz,
        rbmd::Real* flat_virial,
        rbmd::Real* energy_improper);
  };

 ///////////////////////


  template <>
  struct SpecialLJCutCoulForceOp<device::DEVICE_GPU>
  {
    void operator()( Box box,ERFTable* erf_table,
          const rbmd::Real cut_off,
          const rbmd::Id num_atoms,
          const rbmd::Real alpha,
          const rbmd::Real qqr2e,
          const rbmd::Id* atoms_type,
          const rbmd::Id* atoms_id,
          const rbmd::Real* sigma,
          const rbmd::Real* eps,
          const rbmd::Id* start_id,
          const rbmd::Id* end_id,
          const rbmd::Id* id_verletlist,
          const rbmd::Id*  special_ids,
          const rbmd::Real*  special_weights,
          const rbmd::Id*  special_offset,
          const rbmd::Id*  special_count,
          const rbmd::Real* charge,
          const rbmd::Real* px,
          const rbmd::Real* py,
          const rbmd::Real* pz,
          rbmd::Real* fx,
          rbmd::Real* fy,
          rbmd::Real* fz,
          rbmd::Real* flat_virial,
          rbmd::Real* total_evdwl,
          rbmd::Real* total_ecoul);
  };


  template <>
   struct SpecialLJCutCoulRBLForceOp<device::DEVICE_GPU>
  {
    void operator()( Box box,ERFTable* erf_table,
     const rbmd::Real rs,
     const rbmd::Real rc,
     const rbmd::Id num_atoms,
     const rbmd::Id neighbor_sample_num,
     const rbmd::Id pice_num,
     const rbmd::Real alpha,
     const rbmd::Real qqr2e,
     const rbmd::Id* atoms_type,
     const rbmd::Id* atoms_id,
     const rbmd::Real* sigma,
     const rbmd::Real* eps,
     const rbmd::Id* start_id,
     const rbmd::Id* end_id,
     const rbmd::Id* id_verletlist,
     const rbmd::Id* id_random_neighbor,
     const rbmd::Id* random_neighbor_num,
     const rbmd::Id*  special_ids,
     const rbmd::Real*  special_weights,
     const rbmd::Id*  special_offset,
     const rbmd::Id*  special_count,
     const rbmd::Real* charge,
     const rbmd::Real* px,
     const rbmd::Real* py,
     const rbmd::Real* pz,
     rbmd::Real* fx,
     rbmd::Real* fy,
     rbmd::Real* fz);
  };



  template <>
  struct SpeciaLJCutCoulEnergyOp<device::DEVICE_GPU>
  {
    void operator()( Box box,ERFTable* erf_table,
         const rbmd::Real cut_off,
         const rbmd::Id num_atoms,
         const rbmd::Real alpha,
         const rbmd::Real qqr2e,
         const rbmd::Id* atoms_type,
         const rbmd::Id* atoms_id,
         const rbmd::Real* sigma,
         const rbmd::Real* eps,
         const rbmd::Id* start_id,
         const rbmd::Id* end_id,
         const rbmd::Id* id_verletlist,
         const rbmd::Id*  special_ids,
         const rbmd::Real*  special_weights,
         const rbmd::Id*  special_offset,
         const rbmd::Id*  special_count,
         const rbmd::Real* charge,
         const rbmd::Real* px,
         const rbmd::Real* py,
         const rbmd::Real* pz,
         rbmd::Real* flat_virial,
         rbmd::Real* total_evdwl,
         rbmd::Real* total_ecoul);
  };


  template <>
  struct  ComputeSpecialCoulForceOp<device::DEVICE_GPU>
  {
    void  operator()(
       Box box,
      const rbmd::Id num_atoms,
      const rbmd::Real qqr2e,
      const rbmd::Id* atoms_id,
      const rbmd::Id*  atoms_vec,
      const rbmd::Id*  atoms_offset,
      const rbmd::Id*  atom_count,
      const rbmd::Id*  special_ids,
      const rbmd::Real*  special_weights,
      const rbmd::Id*  special_offset,
      const rbmd::Id*  special_count,
      const rbmd::Real* charge,
      const rbmd::Real* px,
      const rbmd::Real* py,
      const rbmd::Real* pz,
      rbmd::Real* fx,
      rbmd::Real* fy,
      rbmd::Real* fz,
      rbmd::Real* flat_virial,
      rbmd::Real* total_especial_coul);
  };

  template <>
  struct  ComputeBondForceOp<device::DEVICE_GPU>
  {
    void  operator()(
       Box box,
       const rbmd::Id num_atoms,
      const rbmd::Id num_bonds,
      const rbmd::Id* atom_id_to_idx,
      const rbmd::Real* bond_coeffs_k,
      const rbmd::Real* bond_coeffs_equilibrium,
      const rbmd::Id* bond_type,
      const rbmd::Id* bondlisti,
      const rbmd::Id* bondlistj,
      const rbmd::Real* px,
      const rbmd::Real* py,
      const rbmd::Real* pz,
      rbmd::Real* fx,
      rbmd::Real* fy,
      rbmd::Real* fz,
      rbmd::Real* flat_virial,
      rbmd::Real* global_virial,
      rbmd::Real* energy_bond);
  };


  template <>
  struct ComputeAngleForceOp<device::DEVICE_GPU>
  {
    void operator()(
     Box box,
     const rbmd::Id num_atoms,
    const rbmd::Id num_anglels,
    const rbmd::Id* atom_id_to_idx,
    const rbmd::Real* anglel_coeffs_k,
    const rbmd::Real* anglel_coeffs_equilibrium,
    const rbmd::Id* anglel_type,
    const rbmd::Id* anglelisti,
    const rbmd::Id* anglelistj,
    const rbmd::Id* anglelistk,
    const rbmd::Real* px,
    const rbmd::Real* py,
    const rbmd::Real* pz,
    rbmd::Real* fx,
    rbmd::Real* fy,
    rbmd::Real* fz,
    rbmd::Real* flat_virial,
    rbmd::Real* global_virial,
    rbmd::Real* energy_angle);
  };

  template <>
  struct ComputeDihedralForceOp<device::DEVICE_GPU>
  {
    void operator()(
       Box box,
       const rbmd::Id num_atoms,
      const rbmd::Id num_dihedrals,
      const rbmd::Id* atom_id_to_idx,
      const rbmd::Real* dihedral_coeffs_k,
      const rbmd::Id* dihedral_coeffs_sign ,
      const rbmd::Id* dihedral_coeffs_multiplicity ,
      const rbmd::Id* dihedral_type,
      const rbmd::Id* dihedrallisti,
      const rbmd::Id* dihedrallistj,
      const rbmd::Id* dihedrallistk,
      const rbmd::Id* dihedrallistw,
      const rbmd::Real* px,
      const rbmd::Real* py,
      const rbmd::Real* pz,
      rbmd::Real* fx,
      rbmd::Real* fy,
      rbmd::Real* fz,
      rbmd::Real* flat_virial,
      rbmd::Real* global_virial,
      rbmd::Real* energy_dihedral);
  };

  template <>
  struct ComputeDihedralOPLSForceOp<device::DEVICE_GPU>
  {
    void operator()(
       Box box,
       const rbmd::Id num_atoms,
      const rbmd::Id num_dihedrals,
      const rbmd::Id* atom_id_to_idx,
      const rbmd::Real* dihedral_coeffs_k1,
      const rbmd::Real* dihedral_coeffs_k2,
      const rbmd::Real* dihedral_coeffs_k3,
      const rbmd::Real* dihedral_coeffs_k4,
      const rbmd::Id* dihedral_type,
      const rbmd::Id* dihedrallisti,
      const rbmd::Id* dihedrallistj,
      const rbmd::Id* dihedrallistk,
      const rbmd::Id* dihedrallistw,
      const rbmd::Real* px,
      const rbmd::Real* py,
      const rbmd::Real* pz,
      rbmd::Real* fx,
      rbmd::Real* fy,
      rbmd::Real* fz,
      rbmd::Real* flat_virial,
      rbmd::Real* global_virial,
      rbmd::Real* energy_dihedral);
  };

  template <>
  struct ComputeImproperHarmonicForceOp<device::DEVICE_GPU>
  {
    void operator()(
       Box box,
       const rbmd::Id num_atoms,
      const rbmd::Id num_impropers,
      const rbmd::Id* atom_id_to_idx,
      const rbmd::Real* improper_coeffs_k,
      const rbmd::Real* improper_coeffs_chi,
      const rbmd::Id* improper_type,
      const rbmd::Id* improperlisti,
      const rbmd::Id* improperlistj,
      const rbmd::Id* improperlistk,
      const rbmd::Id* improperlistw,
      const rbmd::Real* px,
      const rbmd::Real* py,
      const rbmd::Real* pz,
      rbmd::Real* fx,
      rbmd::Real* fy,
      rbmd::Real* fz,
      rbmd::Real* flat_virial,
      rbmd::Real* energy_improper);
  };





}