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
        struct LJCutCoulForceUserOp
        {
          void operator()( Box box, const rbmd::Real cut_off, const rbmd::Id num_atoms,const  rbmd::Real qqr2e,
         const rbmd::Real rbsog_sigma,const rbmd::Real rbsog_b, const rbmd::Id rbsog_mmax,
         const rbmd::Real rbsog_w0,const rbmd::Real* taylor_coeff,
         const rbmd::Id* atoms_type,
         const rbmd::Real* sigma, const rbmd::Real* eps,
         const rbmd::Id* start_id, const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
         const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
         rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
         rbmd::Real* flat_virial, rbmd::Real* total_evdwl, rbmd::Real* total_ecoul);
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
      struct LJCutCoulForceUserOp<device::DEVICE_GPU>
      {
        void operator()( Box box, const rbmd::Real cut_off, const rbmd::Id num_atoms,const  rbmd::Real qqr2e,
       const rbmd::Real rbsog_sigma,const rbmd::Real rbsog_b, const rbmd::Id rbsog_mmax,
       const rbmd::Real rbsog_w0,const rbmd::Real* taylor_coeff,
       const rbmd::Id* atoms_type,
       const rbmd::Real* sigma, const rbmd::Real* eps,
       const rbmd::Id* start_id, const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
       const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
       rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
       rbmd::Real* flat_virial, rbmd::Real* total_evdwl, rbmd::Real* total_ecoul);
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

}// namespace op