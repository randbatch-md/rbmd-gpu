#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"
#include "../common/erf_table.h"

namespace op
{

        template <typename DEVICE>
        struct LJForceOp
        {
          void operator()(Box* box,
                          const rbmd::Real cut_off,
                          const rbmd::Id num_atoms,
                          const rbmd::Id* atoms_type,
                          const rbmd::Id* molecular_type,
                          const rbmd::Real* sigma,
                          const rbmd::Real* eps,
                          const rbmd::Id* start_id,
                          const rbmd::Id* end_id,
                          const rbmd::Id* id_verletlist,
                          const rbmd::Real* px,
                          const rbmd::Real* py,
                          const rbmd::Real* pz,
                          rbmd::Real* fx,
                          rbmd::Real* fy,
                          rbmd::Real* fz,
                          rbmd::Real* total_evdwl);
        };

	template <typename DEVICE>
	struct LJForceVirialOp
	{
		void operator()(Box* box,
			const rbmd::Real cut_off,
			const rbmd::Id num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Id* molecular_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz,
			rbmd::Real* virial_xx,
			rbmd::Real* virial_yy,
			rbmd::Real* virial_zz,
			rbmd::Real* virial_xy,
			rbmd::Real* virial_xz,
			rbmd::Real* virial_yz,
			rbmd::Real* total_evdwl);
	};

	template <typename DEVICE>
	struct LJRBLForceOp
	{
		void operator()(Box* box,
			const rbmd::Real rs,
			const rbmd::Real rc,
			const rbmd::Id num_atoms,
			const rbmd::Id neighbor_sample_num,
			const rbmd::Id pice_num,
			const rbmd::Id* atoms_type,
			const rbmd::Id* molecular_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Id* id_random_neighbor,
			const rbmd::Id* random_neighbor_num,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

	template <typename DEVICE>
	struct FixRBLForceOp
	{
		void operator()(
			const rbmd::Id num_atoms,
			const rbmd::Real corr_value_x,
			const rbmd::Real corr_value_y,
			const rbmd::Real corr_value_z,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

	template <typename DEVICE>
	struct LJEnergyOp
	{
		void operator()(Box* box,
			const rbmd::Real cut_off,
			const rbmd::Id num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Id* molecular_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* total_evdwl);
	};

        template <typename DEVICE>
        struct LJCutCoulForceOp
        {
          void operator()(Box* box,ERFTable* erf_table,
            const rbmd::Real cut_off,
            const rbmd::Id num_atoms,
            const rbmd::Real alpha,
            const rbmd::Real qqr2e,
            const rbmd::Id* atoms_type,
            const rbmd::Id* molecular_type,
            const rbmd::Real* sigma,
            const rbmd::Real* eps,
            const rbmd::Id* start_id,
            const rbmd::Id* end_id,
            const rbmd::Id* id_verletlist,
            const rbmd::Real* charge,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* fx,
            rbmd::Real* fy,
            rbmd::Real* fz,
            rbmd::Real* total_evdwl,
            rbmd::Real* total_ecoul);
        };

       template <typename DEVICE>
       struct LJCutCoulRBLForceOp
       {
         void operator()(Box* box,
                  const rbmd::Real rs,
                  const rbmd::Real rc,
                  const rbmd::Id num_atoms,
                  const rbmd::Id neighbor_sample_num,
                  const rbmd::Id pice_num,
                  const rbmd::Real alpha,
                  const rbmd::Real qqr2e,
                  const rbmd::Id* atoms_type,
                  const rbmd::Id* molecular_type,
                  const rbmd::Real* sigma,
                  const rbmd::Real* eps,
                  const rbmd::Id* start_id,
                  const rbmd::Id* end_id,
                  const rbmd::Id* id_verletlist,
                  const rbmd::Id* id_random_neighbor,
                  const rbmd::Id* random_neighbor_num,
                  const rbmd::Real* charge,
                  const rbmd::Real* px,
                  const rbmd::Real* py,
                  const rbmd::Real* pz,
                  rbmd::Real* fx,
                  rbmd::Real* fy,
                  rbmd::Real* fz);
       };

        template <typename DEVICE>
        struct LJCutCoulEnergyOp
        {
          void operator()(Box* box,ERFTable* erf_table,
               const rbmd::Real cut_off,
               const rbmd::Id num_atoms,
               const rbmd::Real alpha,
               const rbmd::Real qqr2e,
               const rbmd::Id* atoms_type,
               const rbmd::Id* molecular_type,
               const rbmd::Real* sigma,
               const rbmd::Real* eps,
               const rbmd::Id* start_id,
               const rbmd::Id* end_id,
               const rbmd::Id* id_verletlist,
               const rbmd::Real* charge,
               const rbmd::Real* px,
               const rbmd::Real* py,
               const rbmd::Real* pz,
               rbmd::Real* total_evdwl,
               rbmd::Real* total_ecoul);
        };

	template <typename DEVICE>
	struct ComputeChargeStructureFactorOp
	{
		void operator()(
			const rbmd::Id num_atoms,
			const Real3 K,
			const rbmd::Real* charge,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* density_real,
			rbmd::Real* density_imag);
	};

	template <typename DEVICE>
	struct ComputeEwaldForceOp
	{
		void operator()(
			Box* box,
			const rbmd::Id num_atoms,
			const rbmd::Id  Kmax,
			const rbmd::Real alpha,
			const rbmd::Real qqr2e,
			const rbmd::Real* real_array,
			const rbmd::Real* imag_array,
			const rbmd::Real* charge,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

       template<typename DEVICE>
       struct SqchargeOp
        {
         void operator()(
           const rbmd::Id num_atoms,
           const rbmd::Real* charge,
           rbmd::Real* sq_charge);
       };


        template<typename DEVICE>
        struct GenerateIndexArrayOp
        {
          void operator()(
          const rbmd::Id  num_atoms,
          const rbmd::Id  RBE_P,
          rbmd::Id* psample_key);
        };

       //RBE
       template <typename DEVICE>
       struct ComputePnumberChargeStructureFactorOp
       {
         void operator()(
            Box* box,
            const rbmd::Id num_atoms,
            const rbmd::Id p_number,
            const rbmd::Real* charge,
            const rbmd::Real* p_sample_x,
            const rbmd::Real* p_sample_y,
            const rbmd::Real* p_sample_z,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* density_real,
            rbmd::Real* density_imag);
       };

      template <typename DEVICE>
      struct ComputeRBEForceOp
      {
        void operator()(
           Box* box,
           const rbmd::Id num_atoms,
           const rbmd::Id  p_number,
           const rbmd::Real alpha,
           const rbmd::Real qqr2e,
           const rbmd::Real* real_array,
           const rbmd::Real* imag_array,
           const rbmd::Real* charge,
           const rbmd::Real* p_sample_x,
           const rbmd::Real* p_sample_y,
           const rbmd::Real* p_sample_z,
           const rbmd::Real* px,
           const rbmd::Real* py,
           const rbmd::Real* pz,
           rbmd::Real* fx,
           rbmd::Real* fy,
           rbmd::Real* fz);
      };

     template <typename DEVICE>
     struct  ComputeBondForceOp
     {
       void  operator()(
         Box* box,
         const rbmd::Id num_bonds,
         const rbmd::Id num_atoms,
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
         rbmd::Real* energy_bond);
     };

     template <typename DEVICE>
     struct ComputeAngleForceOp
     {
       void operator()(
       Box* box,
       const rbmd::Id num_anglels,
       const rbmd::Id num_atoms,
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
       rbmd::Real* energy_angle);
     };

     template <typename DEVICE>
     struct ComputeDihedralForceOp
     {
       void operator()(
         Box* box,
         const rbmd::Id num_dihedrals,
         const rbmd::Id num_atoms,
         const rbmd::Real* dihedral_coeffs_k,
         const rbmd::Id* dihedral_coeffs_sign ,
         const rbmd::Id* dihedral_coeffs_multiplicity ,
         const rbmd::Id* dihedral_type,
         const int4* dihedrallist,
         const rbmd::Real* px,
         const rbmd::Real* py,
         const rbmd::Real* pz,
         rbmd::Real* fx,
         rbmd::Real* fy,
         rbmd::Real* fz,
         rbmd::Real* energy_dihedral);
     };

     template <typename DEVICE>
     struct  ComputeSpecialCoulForceOp
      {
       void  operator()(Box* box,
         const rbmd::Id num_atoms,
         const rbmd::Id*  group_vec,
         const rbmd::Id*  special_ids,
         const rbmd::Id*  special_weights,
         const rbmd::Real* charge,
         const rbmd::Real* px,
         const rbmd::Real* py,
         const rbmd::Real* pz,
         rbmd::Real* fx,
         rbmd::Real* fy,
         rbmd::Real* fz);
     };

        // // // // // // // // // // // // // // // // // // //
        template <>
        struct LJForceOp<device::DEVICE_GPU>
        {
          void operator()(Box* box,
                          const rbmd::Real cut_off,
                          const rbmd::Id num_atoms,
                          const rbmd::Id* atoms_type,
                          const rbmd::Id* molecular_type,
                          const rbmd::Real* sigma,
                          const rbmd::Real* eps,
                          const rbmd::Id* start_id,
                          const rbmd::Id* end_id,
                          const rbmd::Id* id_verletlist,
                          const rbmd::Real* px,
                          const rbmd::Real* py,
                          const rbmd::Real* pz,
                          rbmd::Real* fx,
                          rbmd::Real* fy,
                          rbmd::Real* fz,
                          rbmd::Real* total_evdwl);
        };

	template <>
	struct LJForceVirialOp<device::DEVICE_GPU>
	{
		void operator()(Box* box,
			const rbmd::Real cut_off,
			const rbmd::Id num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Id* molecular_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz,
			rbmd::Real* virial_xx,
			rbmd::Real* virial_yy,
			rbmd::Real* virial_zz,
			rbmd::Real* virial_xy,
			rbmd::Real* virial_xz,
			rbmd::Real* virial_yz,
			rbmd::Real* total_evdwl);
	};

	template <>
	struct LJRBLForceOp<device::DEVICE_GPU>
	{
		void operator()(Box* box,
			const rbmd::Real rs,
			const rbmd::Real rc,
			const rbmd::Id num_atoms,
			const rbmd::Id neighbor_sample_num,
			const rbmd::Id pice_num,
			const rbmd::Id* atoms_type,
			const rbmd::Id* molecular_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Id* id_random_neighbor,
			const rbmd::Id* random_neighbor_num,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

	template <>
	struct FixRBLForceOp<device::DEVICE_GPU>
	{
		void operator()(
			const rbmd::Id num_atoms,
			const rbmd::Real corr_value_x,
			const rbmd::Real corr_value_y,
			const rbmd::Real corr_value_z,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

	template <>
	struct LJEnergyOp<device::DEVICE_GPU>
	{
		void operator()(Box* box,
			const rbmd::Real cut_off,
			const rbmd::Id num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Id* molecular_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* total_evdwl);
	};

        template <>
        struct LJCutCoulForceOp<device::DEVICE_GPU>
        {
          void operator()(Box* box,ERFTable* erf_table,
            const rbmd::Real cut_off,
            const rbmd::Id num_atoms,
            const rbmd::Real alpha,
            const rbmd::Real qqr2e,
            const rbmd::Id* atoms_type,
            const rbmd::Id* molecular_type,
            const rbmd::Real* sigma,
            const rbmd::Real* eps,
            const rbmd::Id* start_id,
            const rbmd::Id* end_id,
            const rbmd::Id* id_verletlist,
            const rbmd::Real* charge,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* fx,
            rbmd::Real* fy,
            rbmd::Real* fz,
            rbmd::Real* total_evdwl,
            rbmd::Real* total_ecoul);
        };

     template <>
     struct LJCutCoulRBLForceOp<device::DEVICE_GPU>
     {
       void operator()(Box* box,ERFTable* erf_table,
                const rbmd::Real rs,
                const rbmd::Real rc,
                const rbmd::Id num_atoms,
                const rbmd::Id neighbor_sample_num,
                const rbmd::Id pice_num,
                const rbmd::Real alpha,
                const rbmd::Real qqr2e,
                const rbmd::Id* atoms_type,
                const rbmd::Id* molecular_type,
                const rbmd::Real* sigma,
                const rbmd::Real* eps,
                const rbmd::Id* start_id,
                const rbmd::Id* end_id,
                const rbmd::Id* id_verletlist,
                const rbmd::Id* id_random_neighbor,
                const rbmd::Id* random_neighbor_num,
                const rbmd::Real* charge,
                const rbmd::Real* px,
                const rbmd::Real* py,
                const rbmd::Real* pz,
                rbmd::Real* fx,
                rbmd::Real* fy,
                rbmd::Real* fz);
     };

       template <>
       struct LJCutCoulEnergyOp<device::DEVICE_GPU>
       {
         void operator()(Box* box,ERFTable* erf_table,
              const rbmd::Real cut_off,
              const rbmd::Id num_atoms,
              const rbmd::Real alpha,
              const rbmd::Real qqr2e,
              const rbmd::Id* atoms_type,
              const rbmd::Id* molecular_type,
              const rbmd::Real* sigma,
              const rbmd::Real* eps,
              const rbmd::Id* start_id,
              const rbmd::Id* end_id,
              const rbmd::Id* id_verletlist,
              const rbmd::Real* charge,
              const rbmd::Real* px,
              const rbmd::Real* py,
              const rbmd::Real* pz,
              rbmd::Real* total_evdwl,
              rbmd::Real* total_ecoul);
       };

	template <>
	struct ComputeChargeStructureFactorOp<device::DEVICE_GPU>
	{
		void operator()(
			const rbmd::Id num_atoms,
			const Real3 K,
			const rbmd::Real* charge,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* density_real,
			rbmd::Real* density_imag);
	};

	template <>
	struct ComputeEwaldForceOp<device::DEVICE_GPU>
	{
		void operator()(
			Box* box,
			const rbmd::Id num_atoms,
			const rbmd::Id  Kmax,
			const rbmd::Real alpha,
			const rbmd::Real qqr2e,
			const rbmd::Real* real_array,
			const rbmd::Real* imag_array,
			const rbmd::Real* charge,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

        template<>
        struct SqchargeOp<device::DEVICE_GPU>
        {
          void operator()(
            const rbmd::Id num_atoms,
            const rbmd::Real* charge,
            rbmd::Real* sq_charge);
        };


        //
       template<>
       struct GenerateIndexArrayOp<device::DEVICE_GPU>
       {
         void operator()(
         const rbmd::Id  num_atoms,
         const rbmd::Id  RBE_P,
         rbmd::Id* psample_key);
       };

        //RBE
       template <>
       struct ComputePnumberChargeStructureFactorOp<device::DEVICE_GPU>
       {
         void operator()(
            Box* box,
            const rbmd::Id num_atoms,
            const rbmd::Id p_number,
            const rbmd::Real* charge,
            const rbmd::Real* p_sample_x,
            const rbmd::Real* p_sample_y,
            const rbmd::Real* p_sample_z,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* density_real,
            rbmd::Real* density_imag);
       };

      template <>
      struct ComputeRBEForceOp<device::DEVICE_GPU>
      {
        void operator()(
           Box* box,
           const rbmd::Id num_atoms,
           const rbmd::Id  p_number,
           const rbmd::Real alpha,
           const rbmd::Real qqr2e,
           const rbmd::Real* real_array,
           const rbmd::Real* imag_array,
           const rbmd::Real* charge,
           const rbmd::Real* p_sample_x,
           const rbmd::Real* p_sample_y,
           const rbmd::Real* p_sample_z,
           const rbmd::Real* px,
           const rbmd::Real* py,
           const rbmd::Real* pz,
           rbmd::Real* fx,
           rbmd::Real* fy,
           rbmd::Real* fz);
      };

    template <>
    struct  ComputeSpecialCoulForceOp<device::DEVICE_GPU>
   {
        void  operator()(Box* box,
          const rbmd::Id num_atoms,
          const rbmd::Id*  group_vec,
          const rbmd::Id*  special_ids,
          const rbmd::Id*  special_weights,
          const rbmd::Real* charge,
          const rbmd::Real* px,
          const rbmd::Real* py,
          const rbmd::Real* pz,
          rbmd::Real* fx,
          rbmd::Real* fy,
          rbmd::Real* fz);
    };



     template <>
     struct  ComputeBondForceOp<device::DEVICE_GPU>
     {
       void  operator()(
         Box* box,
         const rbmd::Id num_bonds,
         const rbmd::Id num_atoms,
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
         rbmd::Real* energy_bond);
     };

     template <>
     struct ComputeAngleForceOp<device::DEVICE_GPU>
    {
       void operator()(
       Box* box,
       const rbmd::Id num_anglels,
       const rbmd::Id num_atoms,
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
       rbmd::Real* energy_angle);
     };

    template <>
    struct ComputeDihedralForceOp<device::DEVICE_GPU>
    {
      void operator()(
        Box* box,
        const rbmd::Id num_dihedrals,
        const rbmd::Id num_atoms,
        const rbmd::Real* dihedral_coeffs_k,
        const rbmd::Id* dihedral_coeffs_sign ,
        const rbmd::Id* dihedral_coeffs_multiplicity ,
        const rbmd::Id* dihedral_type,
        const int4* dihedrallist,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* fx,
        rbmd::Real* fy,
        rbmd::Real* fz,
        rbmd::Real* energy_dihedral);
    };



}// namespace op