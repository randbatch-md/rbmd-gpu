#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"
#include "../common/erf_table.h"

namespace op
{
        template <typename DEVICE>
        struct LJCutCoulForceOp
        {
          void operator()( Box box,ERFTable* erf_table,
                const rbmd::Real cut_off,
                const rbmd::Id num_atoms,
                const rbmd::Real alpha,
                const rbmd::Real qqr2e,
                const rbmd::Id* atoms_type,
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
                rbmd::Real* flat_virial,
                rbmd::Real* total_evdwl,
                rbmd::Real* total_ecoul);
        };


       template <typename DEVICE>
       struct LJCutCoulRBLForceOp
       {
         void operator()( Box box,
                  const rbmd::Real rs,
                  const rbmd::Real rc,
                  const rbmd::Id num_atoms,
                  const rbmd::Id neighbor_sample_num,
                  const rbmd::Id pice_num,
                  const rbmd::Real alpha,
                  const rbmd::Real qqr2e,
                  const rbmd::Id* atoms_type,
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
          void operator()( Box box,ERFTable* erf_table,
               const rbmd::Real cut_off,
               const rbmd::Id num_atoms,
               const rbmd::Real alpha,
               const rbmd::Real qqr2e,
               const rbmd::Id* atoms_type,
               const rbmd::Real* sigma,
               const rbmd::Real* eps,
               const rbmd::Id* start_id,
               const rbmd::Id* end_id,
               const rbmd::Id* id_verletlist,
               const rbmd::Real* charge,
               const rbmd::Real* px,
               const rbmd::Real* py,
               const rbmd::Real* pz,
               rbmd::Real*  flat_virial,
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
			 Box box,
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
                        rbmd::Real* fz,
                        rbmd::Real* flat_virial);
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
             Box box,
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
            Box box,
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
           rbmd::Real* fz,
           rbmd::Real* flat_virial);
      };

      template <typename DEVICE>
      struct AddForceOp
      {
        void operator()(
                const rbmd::Id num_atoms,
                const rbmd::Real* input_fx,
                const rbmd::Real* input_fy,
                const rbmd::Real* input_fz,
                rbmd::Real* fx,
                rbmd::Real* fy,
                rbmd::Real* fz);
      };

    template <typename DEVICE>
    struct EikOp
    {
      void operator()(
      const rbmd::Id num_atoms,const rbmd::Real gsqmx,Real3 unitk,
      const rbmd::Id kmax,Int3 kmax_array,const rbmd::Real* px,
      const rbmd::Real* py,const rbmd::Real* pz,const rbmd::Real* charge,
      rbmd::Real* cs, rbmd::Real* sn,rbmd::Real* sfacrl, rbmd::Real* sfacim);
    };

    template <typename DEVICE>
    struct EwaldForceFixOp
    {
      void operator()(
      const rbmd::Id num_atoms,const rbmd::Id kcount,const rbmd::Id k_index,
      const rbmd::Real qqr2e,Int3 kmax_vec3D,const rbmd::Real* eg,
      const rbmd::Real* cs,const rbmd::Real* sn,const rbmd::Real* charge,
      const rbmd::Real* qfactor_real,const rbmd::Real* qfactor_image,
      rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz);

    };
        // // // // // // // // // // // // // // // // // // //
        template <>
        struct LJCutCoulForceOp<device::DEVICE_GPU>
        {
          void operator()( Box box,ERFTable* erf_table,
                const rbmd::Real cut_off,
                const rbmd::Id num_atoms,
                const rbmd::Real alpha,
                const rbmd::Real qqr2e,
                const rbmd::Id* atoms_type,
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
                rbmd::Real* flat_virial,
                rbmd::Real* total_evdwl,
                rbmd::Real* total_ecoul);
        };

     template <>
     struct LJCutCoulRBLForceOp<device::DEVICE_GPU>
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
         void operator()( Box box,ERFTable* erf_table,
              const rbmd::Real cut_off,
              const rbmd::Id num_atoms,
              const rbmd::Real alpha,
              const rbmd::Real qqr2e,
              const rbmd::Id* atoms_type,
              const rbmd::Real* sigma,
              const rbmd::Real* eps,
              const rbmd::Id* start_id,
              const rbmd::Id* end_id,
              const rbmd::Id* id_verletlist,
              const rbmd::Real* charge,
              const rbmd::Real* px,
              const rbmd::Real* py,
              const rbmd::Real* pz,
              rbmd::Real*  flat_virial,
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
			 Box box,
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
                        rbmd::Real* fz,
                        rbmd::Real* flat_virial);
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
             Box box,
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
            Box box,
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
           rbmd::Real* fz,
           rbmd::Real* flat_virial);
      };

    template <>
    struct AddForceOp<device::DEVICE_GPU>
    {
      void operator()(
              const rbmd::Id num_atoms,
              const rbmd::Real* input_fx,
              const rbmd::Real* input_fy,
              const rbmd::Real* input_fz,
              rbmd::Real* fx,
              rbmd::Real* fy,
              rbmd::Real* fz);
    };

  template <>
  struct EikOp<device::DEVICE_GPU>
  {
    void operator()(
    const rbmd::Id num_atoms,const rbmd::Real gsqmx,Real3 unitk,
    const rbmd::Id kmax,Int3 kmax_array,const rbmd::Real* px,
    const rbmd::Real* py,const rbmd::Real* pz,const rbmd::Real* charge,
    rbmd::Real* cs, rbmd::Real* sn,rbmd::Real* sfacrl, rbmd::Real* sfacim);
  };

  template <>
  struct EwaldForceFixOp<device::DEVICE_GPU>
  {
    void operator()(
    const rbmd::Id num_atoms,const rbmd::Id kcount, const rbmd::Id k_index,
    const rbmd::Real qqr2e, Int3 kmax_vec3D,const rbmd::Real* eg,
    const rbmd::Real* cs,const rbmd::Real* sn,const rbmd::Real* charge,
    const rbmd::Real* qfactor_real,const rbmd::Real* qfactor_image,
    rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz);

  };


}// namespace op