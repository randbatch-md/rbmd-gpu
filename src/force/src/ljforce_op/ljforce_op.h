#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"

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
	struct FixLJRBLForceOp
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
          void operator()(Box* box,
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
	struct FixLJRBLForceOp<device::DEVICE_GPU>
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
          void operator()(Box* box,
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

}// namespace op