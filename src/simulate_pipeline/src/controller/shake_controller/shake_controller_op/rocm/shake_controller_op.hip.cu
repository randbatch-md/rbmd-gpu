#include <hip/hip_runtime.h>

#include "rbmd_define.h"
#include "shake_controller_op.h"
#include <math.h>
namespace op {
#define THREADS_PER_BLOCK 256

__device__ rbmd::Real Dot(const Real3& p_1,
                          const Real3& p_2)
    {
        return p_1.x * p_2.x + p_1.y * p_2.y + p_1.z * p_2.z;
    }

    __device__ bool IsNan(const rbmd::Real& value) {
        //return value != value;  // TODO: isnan
        return (isnan(value) != 0);
}

__global__ void ShakeA(const rbmd::Id num_angle,
                       const rbmd::Real dt,
                       const rbmd::Real fmt2v,
                       Box* box,
                       const rbmd::Id* atom_id_to_idx,
                       const rbmd::Real* mass,
                       const rbmd::Id* atoms_type,
                       const Id3* angle_id_vec,
                       rbmd::Real* shake_px,
                       rbmd::Real* shake_py,
                       rbmd::Real* shake_pz,
                       rbmd::Real* shake_vx,
                       rbmd::Real* shake_vy,
                       rbmd::Real* shake_vz,
                       const rbmd::Real* fx,
                       const rbmd::Real* fy,
                       const rbmd::Real* fz,
                       rbmd::Id* flag_px,
                       rbmd::Id* flag_py,
                       rbmd::Id* flag_pz) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid < num_angle) {

        rbmd::Id id_0x = angle_id_vec[tid].x;
        rbmd::Id id_1x = angle_id_vec[tid].y;
        rbmd::Id id_2x = angle_id_vec[tid].z;

    rbmd::Id id_0 = atom_id_to_idx[id_0x];
    rbmd::Id id_1 = atom_id_to_idx[id_1x];
    rbmd::Id id_2 = atom_id_to_idx[id_2x];

    Real3 shake_position_0, shake_position_1, shake_position_2;
    shake_position_0.x = shake_px[id_0] + shake_vx[id_0] * dt + 0.5 * dt * dt * fx[id_0] / mass[atoms_type[id_0]] * fmt2v;
    shake_position_0.y = shake_py[id_0] + shake_vy[id_0] * dt + 0.5 * dt * dt * fy[id_0] / mass[atoms_type[id_0]] * fmt2v;
    shake_position_0.z = shake_pz[id_0] + shake_vz[id_0] * dt + 0.5 * dt * dt * fz[id_0] / mass[atoms_type[id_0]] * fmt2v;

    shake_position_1.x = shake_px[id_1] + shake_vx[id_1] * dt + 0.5 * dt * dt * fx[id_1] / mass[atoms_type[id_1]] * fmt2v;
    shake_position_1.y = shake_py[id_1] + shake_vy[id_1] * dt + 0.5 * dt * dt * fy[id_1] / mass[atoms_type[id_1]] * fmt2v;
    shake_position_1.z = shake_pz[id_1] + shake_vz[id_1] * dt + 0.5 * dt * dt * fz[id_1] / mass[atoms_type[id_1]] * fmt2v;

    shake_position_2.x = shake_px[id_2] + shake_vx[id_2] * dt + 0.5 * dt * dt * fx[id_2] / mass[atoms_type[id_2]] * fmt2v;
    shake_position_2.y = shake_py[id_2] + shake_vy[id_2] * dt + 0.5 * dt * dt * fy[id_2] / mass[atoms_type[id_2]] * fmt2v;
    shake_position_2.z = shake_pz[id_2] + shake_vz[id_2] * dt + 0.5 * dt * dt * fz[id_2] / mass[atoms_type[id_2]] * fmt2v;

        //printf("tid为：%d----输出结果为：%f, %f, %f \n",tid,shake_position_0.x,shake_position_1.x,shake_position_2.x);


    rbmd::Real bond1 = 1.0;
    rbmd::Real bond2 = 1.0;
    rbmd::Real bond12 = SQRT(bond1 * bond1 + bond2 * bond2 - 2.0 *
      bond1 * bond2 * COS((109.4700 / 180.0) * M_PI));
    
    // minimum image
    Real3 r01 = MinImageDistanceVec(shake_px[id_1], shake_py[id_1],
      shake_pz[id_1], shake_px[id_0], shake_py[id_0], shake_pz[id_0], box);
    Real3 r12 = MinImageDistanceVec(shake_px[id_2], shake_py[id_2],
      shake_pz[id_2], shake_px[id_1], shake_py[id_1], shake_pz[id_1], box);
    Real3 r20 = MinImageDistanceVec(shake_px[id_0], shake_py[id_0],
      shake_pz[id_0], shake_px[id_2], shake_py[id_2], shake_pz[id_2], box);

    // s01,s02,s12 = distance vec after unconstrained update, with PBC
    Real3 s10 = MinImageDistanceVec(shake_position_0.x, shake_position_0.y,
      shake_position_0.z, shake_position_1.x, shake_position_1.y, shake_position_1.z, box);
    Real3 s21 = MinImageDistanceVec(shake_position_1.x, shake_position_1.y,
      shake_position_1.z, shake_position_2.x, shake_position_2.y, shake_position_2.z, box);
    Real3 s02 = MinImageDistanceVec(shake_position_2.x, shake_position_2.y,
      shake_position_2.z, shake_position_0.x, shake_position_0.y, shake_position_0.z, box);

    // scalar distances between atoms
    rbmd::Real r01sq = Dot(r01, r01);
    rbmd::Real r02sq = Dot(r20, r20);
    rbmd::Real r12sq = Dot(r12, r12);
    rbmd::Real s01sq = Dot(s10, s10);
    rbmd::Real s02sq = Dot(s02, s02);
    rbmd::Real s12sq = Dot(s21, s21);

    // matrix coeffs and rhs for lamda equations
    rbmd::Real invmass0 = 1 / mass[atoms_type[id_0]];
    rbmd::Real invmass1 = 1 / mass[atoms_type[id_1]];
    rbmd::Real invmass2 = 1 / mass[atoms_type[id_2]];
    rbmd::Real a11 = 2 * (invmass0 + invmass1) * Dot(s10, r01);
    rbmd::Real a12 = -2 * invmass1 * Dot(s10, r12);
    rbmd::Real a13 = -2 * invmass0 * Dot(s10, r20);
    rbmd::Real a21 = -2 * invmass1 * Dot(s21, r01);
    rbmd::Real a22 = 2 * (invmass1 + invmass2) * Dot(s21, r12);
    rbmd::Real a23 = -2 * invmass2 * Dot(s21, r20);
    rbmd::Real a31 = -2 * invmass0 * Dot(s02, r01);
    rbmd::Real a32 = -2 * invmass2 * Dot(s02, r12);
    rbmd::Real a33 = 2 * (invmass0 + invmass2) * (Dot(s02, r20));

      // inverse of matrix

      rbmd::Real determ = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 -
        a11 * a23 * a32 - a12 * a21 * a33 - a13 * a22 * a31;
      if (ABS(determ) < 0.0001)
      {
          printf("Shake determinant = 0.0");
      }


      rbmd::Real determinv = 1 / determ;

      rbmd::Real a11inv = determinv * (a22 * a33 - a23 * a32);
      rbmd::Real a12inv = -determinv * (a12 * a33 - a13 * a32);
      rbmd::Real a13inv = determinv * (a12 * a23 - a13 * a22);
      rbmd::Real a21inv = -determinv * (a21 * a33 - a23 * a31);
      rbmd::Real a22inv = determinv * (a11 * a33 - a13 * a31);
      rbmd::Real a23inv = -determinv * (a11 * a23 - a13 * a21);
      rbmd::Real a31inv = determinv * (a21 * a32 - a22 * a31);
      rbmd::Real a32inv = -determinv * (a11 * a32 - a12 * a31);
      rbmd::Real a33inv = determinv * (a11 * a22 - a12 * a21);


      rbmd::Real r0120 = Dot(r01, r20);
      rbmd::Real r0112 = Dot(r01, r12);
      rbmd::Real r2012 = Dot(r20, r12);

      rbmd::Real quad1_0101 = (invmass0 + invmass1) * (invmass0 + invmass1) * r01sq;
      rbmd::Real quad1_1212 = invmass1 * invmass1 * r12sq;
      rbmd::Real quad1_2020 = invmass0 * invmass0 * r02sq;
      rbmd::Real quad1_0120 = -2 * (invmass0 + invmass1) * invmass0 * r0120;
      rbmd::Real quad1_0112 = -2 * (invmass0 + invmass1) * invmass1 * r0112;
      rbmd::Real quad1_2012 = 2 * invmass0 * invmass1 * r2012;

      rbmd::Real quad2_0101 = invmass1 * invmass1 * r01sq;
      rbmd::Real quad2_1212 = (invmass1 + invmass2) * (invmass1 + invmass2) * r12sq;
      rbmd::Real quad2_2020 = invmass2 * invmass2 * r02sq;
      rbmd::Real quad2_0120 = 2 * invmass1 * invmass2 * r0120;
      rbmd::Real quad2_0112 = -2 * (invmass1 + invmass2) * invmass1 * r0112;
      rbmd::Real quad2_2012 = -2 * (invmass1 + invmass2) * invmass2 * r2012;

      rbmd::Real quad3_0101 = invmass0 * invmass0 * r01sq;
      rbmd::Real quad3_1212 = invmass2 * invmass2 * r12sq;
      rbmd::Real quad3_2020 = (invmass0 + invmass2) * (invmass0 + invmass2) * r02sq;
      rbmd::Real quad3_0120 = -2 * (invmass0 + invmass2) * invmass0 * r0120;
      rbmd::Real quad3_0112 = 2 * invmass0 * invmass2 * r0112;
      rbmd::Real quad3_2012 = -2 * (invmass0 + invmass2) * invmass2 * r2012;

        // iterate until converged
        rbmd::Real tolerance = 0.00001;     // original 0.001
        rbmd::Id max_iter = 5000; // original: 100

        rbmd::Real lamda01 = 0.0;
        rbmd::Real lamda20 = 0.0;
        rbmd::Real lamda12 = 0.0;
        rbmd::Id niter = 0;
        rbmd::Id done = 0;
        rbmd::Id flag_overflow = 0;
        rbmd::Real quad1, quad2, quad3, b1, b2, b3, lamda01_new, lamda20_new, lamda12_new;

        while (!done && niter < max_iter)
        {

            quad1 = quad1_0101 * lamda01 * lamda01 + quad1_2020 * lamda20 * lamda20 +
                    quad1_1212 * lamda12 * lamda12 + quad1_0120 * lamda01 * lamda20 +
                    quad1_0112 * lamda01 * lamda12 + quad1_2012 * lamda20 * lamda12;

            quad2 = quad2_0101 * lamda01 * lamda01 + quad2_2020 * lamda20 * lamda20 +
                    quad2_1212 * lamda12 * lamda12 + quad2_0120 * lamda01 * lamda20 +
                    quad2_0112 * lamda01 * lamda12 + quad2_2012 * lamda20 * lamda12;

            quad3 = quad3_0101 * lamda01 * lamda01 + quad3_2020 * lamda20 * lamda20 +
                    quad3_1212 * lamda12 * lamda12 + quad3_0120 * lamda01 * lamda20 +
                    quad3_0112 * lamda01 * lamda12 + quad3_2012 * lamda20 * lamda12;

            b1 = bond1 * bond1 - s01sq - quad1;
            b2 = bond2 * bond2 - s12sq - quad2;
            b3 = bond12 * bond12 - s02sq - quad3;


            lamda01_new = a11inv * b1 + a12inv * b2 + a13inv * b3;
            lamda12_new = a21inv * b1 + a22inv * b2 + a23inv * b3;
            lamda20_new = a31inv * b1 + a32inv * b2 + a33inv * b3;


            done = 1;
            if (ABS(lamda01_new - lamda01) > tolerance)
                done = 0;
            if (ABS(lamda20_new - lamda20) > tolerance)
                done = 0;
            if (ABS(lamda12_new - lamda12) > tolerance)
                done = 0;

            lamda01 = lamda01_new;
            lamda20 = lamda20_new;
            lamda12 = lamda12_new;

            if (IsNan(lamda01) || IsNan(lamda20) || IsNan(lamda12) ||
                ABS(lamda01) > 1e20 || ABS(lamda20) > 1e20 || ABS(lamda12) > 1e20)
            {
                done = 1;
                flag_overflow = 1;
            }
            niter++;
        }

        Real3 position_constraint_i0;
        Real3 position_constraint_i1;
        Real3 position_constraint_i2;

        position_constraint_i0.x = lamda01 * r01.x * invmass0 - lamda20 * r20.x * invmass0;
        position_constraint_i0.y = lamda01 * r01.y * invmass0 - lamda20 * r20.y * invmass0;
        position_constraint_i0.z = lamda01 * r01.z * invmass0 - lamda20 * r20.z * invmass0;

        position_constraint_i1.x = lamda12 * r12.x * invmass1 - lamda01 * r01.x * invmass1;
        position_constraint_i1.y = lamda12 * r12.y * invmass1 - lamda01 * r01.y * invmass1;
        position_constraint_i1.z = lamda12 * r12.z * invmass1 - lamda01 * r01.z * invmass1;

        position_constraint_i2.x = lamda20 * r20.x * invmass2 - lamda12 * r12.x * invmass2;
        position_constraint_i2.y = lamda20 * r20.y * invmass2 - lamda12 * r12.y * invmass2;
        position_constraint_i2.z = lamda20 * r20.z * invmass2 - lamda12 * r12.z * invmass2;


        Real3 velocity_constraint_i0,velocity_constraint_i1, velocity_constraint_i2;

        velocity_constraint_i0.x = position_constraint_i0.x / dt;
        velocity_constraint_i0.y = position_constraint_i0.y / dt;
        velocity_constraint_i0.z = position_constraint_i0.z / dt;

        velocity_constraint_i1.x = position_constraint_i1.x / dt;
        velocity_constraint_i1.y = position_constraint_i1.y / dt;
        velocity_constraint_i1.z = position_constraint_i1.z / dt;

        velocity_constraint_i2.x = position_constraint_i2.x / dt;
        velocity_constraint_i2.y = position_constraint_i2.y / dt;
        velocity_constraint_i2.z = position_constraint_i2.z / dt;

        Real3 shake_velocity_0,shake_velocity_1,shake_velocity_2;
        shake_velocity_0.x = shake_vx[id_0] + 0.5 * dt * fx[id_0]/mass[atoms_type[id_0]] * fmt2v;
        shake_velocity_0.y = shake_vy[id_0] + 0.5 * dt * fy[id_0]/mass[atoms_type[id_0]] * fmt2v;
        shake_velocity_0.z = shake_vz[id_0] + 0.5 * dt * fz[id_0]/mass[atoms_type[id_0]] * fmt2v;

        shake_velocity_1.x = shake_vx[id_1] + 0.5 * dt * fx[id_1]/mass[atoms_type[id_1]] * fmt2v;
        shake_velocity_1.y = shake_vy[id_1] + 0.5 * dt * fy[id_1]/mass[atoms_type[id_1]] * fmt2v;
        shake_velocity_1.z = shake_vz[id_1] + 0.5 * dt * fz[id_1]/mass[atoms_type[id_1]] * fmt2v;

        shake_velocity_2.x = shake_vx[id_2] + 0.5 * dt * fx[id_2]/mass[atoms_type[id_2]] * fmt2v;
        shake_velocity_2.y = shake_vy[id_2] + 0.5 * dt * fy[id_2]/mass[atoms_type[id_2]] * fmt2v;
        shake_velocity_2.z = shake_vz[id_2] + 0.5 * dt * fz[id_2]/mass[atoms_type[id_2]] * fmt2v;

        // velocity
        shake_vx[id_0] = shake_velocity_0.x + velocity_constraint_i0.x;
        shake_vy[id_0] = shake_velocity_0.y + velocity_constraint_i0.y;
        shake_vz[id_0] = shake_velocity_0.z + velocity_constraint_i0.z;

        shake_vx[id_1] = shake_velocity_1.x + velocity_constraint_i1.x;
        shake_vy[id_1] = shake_velocity_1.y + velocity_constraint_i1.y;
        shake_vz[id_1] = shake_velocity_1.z + velocity_constraint_i1.z;

        shake_vx[id_2] = shake_velocity_2.x + velocity_constraint_i2.x;
        shake_vy[id_2] = shake_velocity_2.y + velocity_constraint_i2.y;
        shake_vz[id_2] = shake_velocity_2.z + velocity_constraint_i2.z;

        // position
        rbmd::Real shake_px_0 = shake_position_0.x + position_constraint_i0.x;
        rbmd::Real shake_py_0 = shake_position_0.y + position_constraint_i0.y;
        rbmd::Real shake_pz_0 = shake_position_0.z + position_constraint_i0.z;

        rbmd::Real shake_px_1= shake_position_1.x + position_constraint_i1.x;
        rbmd::Real shake_py_1= shake_position_1.y + position_constraint_i1.y;
        rbmd::Real shake_pz_1 = shake_position_1.z + position_constraint_i1.z;

        rbmd::Real shake_px_2= shake_position_2.x + position_constraint_i2.x;
        rbmd::Real shake_py_2 = shake_position_2.y + position_constraint_i2.y;
        rbmd::Real shake_pz_2 = shake_position_2.z + position_constraint_i2.z;

        // pbc
        Real3 whole_position_pbc_0, whole_position_pbc_1, whole_position_pbc_2;

        // position
        // whole_position_pbc_0.x = shake_px[id_0];
        // whole_position_pbc_0.y = shake_py[id_0];
        // whole_position_pbc_0.z = shake_pz[id_0];
        //
        // whole_position_pbc_1.x = shake_px[id_1];
        // whole_position_pbc_1.y = shake_py[id_1];
        // whole_position_pbc_1.z = shake_pz[id_1];
        //
        // whole_position_pbc_2.x = shake_px[id_2];
        // whole_position_pbc_2.y = shake_py[id_2];
        // whole_position_pbc_2.z = shake_pz[id_2];

      whole_position_pbc_0.x = shake_px_0;
      whole_position_pbc_0.y = shake_py_0;
      whole_position_pbc_0.z = shake_pz_0;

      whole_position_pbc_1.x = shake_px_1;
      whole_position_pbc_1.y = shake_py_1;
      whole_position_pbc_1.z = shake_pz_1;

      whole_position_pbc_2.x = shake_px_2;
      whole_position_pbc_2.y = shake_py_2;
      whole_position_pbc_2.z = shake_pz_2;
        //position_flag
        Id3 whole_pts_flag_pbc_0, whole_pts_flag_pbc_1, whole_pts_flag_pbc_2;

        whole_pts_flag_pbc_0.x = flag_px[id_0];
        whole_pts_flag_pbc_0.y = flag_py[id_0];
        whole_pts_flag_pbc_0.z = flag_pz[id_0];

        whole_pts_flag_pbc_1.x = flag_px[id_1];
        whole_pts_flag_pbc_1.y = flag_py[id_1];
        whole_pts_flag_pbc_1.z = flag_pz[id_1];

        whole_pts_flag_pbc_2.x = flag_px[id_2];
        whole_pts_flag_pbc_2.y = flag_py[id_2];
        whole_pts_flag_pbc_2.z = flag_pz[id_2];

        // anglelist[0]
      ApplyPBC(box,whole_position_pbc_0.x,whole_position_pbc_0.y,
        whole_position_pbc_0.z,whole_pts_flag_pbc_0.x,
        whole_pts_flag_pbc_0.y,whole_pts_flag_pbc_0.z);
      // anglelist1]
      ApplyPBC(box,whole_position_pbc_1.x,whole_position_pbc_1.y,
  whole_position_pbc_1.z,whole_pts_flag_pbc_1.x,
  whole_pts_flag_pbc_1.y,whole_pts_flag_pbc_1.z);
      // anglelist[2]
      ApplyPBC(box,whole_position_pbc_2.x,whole_position_pbc_2.y,
  whole_position_pbc_2.z,whole_pts_flag_pbc_2.x,
  whole_pts_flag_pbc_2.y,whole_pts_flag_pbc_2.z);

      shake_px[id_0] = whole_position_pbc_0.x;
      shake_py[id_0] = whole_position_pbc_0.y;
      shake_pz[id_0] = whole_position_pbc_0.z;

      shake_px[id_1] = whole_position_pbc_1.x;
      shake_py[id_1] = whole_position_pbc_1.y;
      shake_pz[id_1] = whole_position_pbc_1.z;

      shake_px[id_2] = whole_position_pbc_2.x;
      shake_py[id_2] = whole_position_pbc_2.y;
      shake_pz[id_2] = whole_position_pbc_2.z;

      //flag
      flag_px[id_0] = whole_pts_flag_pbc_0.x;
      flag_py[id_0] = whole_pts_flag_pbc_0.y;
      flag_pz[id_0] = whole_pts_flag_pbc_0.z;

      flag_px[id_1] = whole_pts_flag_pbc_1.x;
      flag_py[id_1] = whole_pts_flag_pbc_1.y;
      flag_pz[id_1] = whole_pts_flag_pbc_1.z;

      flag_px[id_2] = whole_pts_flag_pbc_2.x;
      flag_py[id_2] = whole_pts_flag_pbc_2.y;
      flag_pz[id_2] = whole_pts_flag_pbc_2.z;

      // atomicAdd(&shake_px[id_0], whole_position_pbc_0.x);
      // atomicAdd(&shake_py[id_0], whole_position_pbc_0.y);
      // atomicAdd(&shake_pz[id_0], whole_position_pbc_0.z);
      //
      // atomicAdd(&shake_px[id_1], whole_position_pbc_1.x);
      // atomicAdd(&shake_py[id_1], whole_position_pbc_1.y);
      // atomicAdd(&shake_pz[id_1], whole_position_pbc_1.z);
      //
      // atomicAdd(&shake_px[id_2], whole_position_pbc_2.x);
      // atomicAdd(&shake_py[id_2], whole_position_pbc_2.y);
      // atomicAdd(&shake_pz[id_2], whole_position_pbc_2.z);

  }
}

__global__ void ShakeB(const rbmd::Id num_angle,
                       const rbmd::Real dt,
                       const rbmd::Real fmt2v,
                       Box* box,
                       const rbmd::Id* atom_id_to_idx,
                       const rbmd::Real* mass,
                       const rbmd::Id* atoms_type,
                       const Id3* angle_id_vec,
                       const rbmd::Real* px,
                       const rbmd::Real* py,
                       const rbmd::Real* pz,
                       rbmd::Real* shake_vx,
                       rbmd::Real* shake_vy,
                       rbmd::Real* shake_vz,
                       const rbmd::Real* fx,
                       const rbmd::Real* fy,
                       const rbmd::Real* fz) {

        int tid = threadIdx.x + blockIdx.x * blockDim.x;

        if (tid < num_angle) {
            rbmd::Id id_0x = angle_id_vec[tid].x;
            rbmd::Id id_1x = angle_id_vec[tid].y;
            rbmd::Id id_2x = angle_id_vec[tid].z;

            rbmd::Id id_0 = atom_id_to_idx[id_0x];
            rbmd::Id id_1 = atom_id_to_idx[id_1x];
            rbmd::Id id_2 = atom_id_to_idx[id_2x];

            Real3 shake_velocity_0,shake_velocity_1,shake_velocity_2;
            shake_velocity_0.x = shake_vx[id_0] + 0.5 * dt * fx[id_0]/mass[atoms_type[id_0]] * fmt2v;
            shake_velocity_0.y = shake_vy[id_0] + 0.5 * dt * fy[id_0]/mass[atoms_type[id_0]] * fmt2v;
            shake_velocity_0.z = shake_vz[id_0] + 0.5 * dt * fz[id_0]/mass[atoms_type[id_0]] * fmt2v;

            shake_velocity_1.x = shake_vx[id_1] + 0.5 * dt * fx[id_1]/mass[atoms_type[id_1]] * fmt2v;
            shake_velocity_1.y = shake_vy[id_1] + 0.5 * dt * fy[id_1]/mass[atoms_type[id_1]] * fmt2v;
            shake_velocity_1.z = shake_vz[id_1] + 0.5 * dt * fz[id_1]/mass[atoms_type[id_1]] * fmt2v;

            shake_velocity_2.x = shake_vx[id_2] + 0.5 * dt * fx[id_2]/mass[atoms_type[id_2]] * fmt2v;
            shake_velocity_2.y = shake_vy[id_2] + 0.5 * dt * fy[id_2]/mass[atoms_type[id_2]] * fmt2v;
            shake_velocity_2.z = shake_vz[id_2] + 0.5 * dt * fz[id_2]/mass[atoms_type[id_2]] * fmt2v;

            // minimum image
            Real3 r01 = MinImageDistanceVec(px[id_1], py[id_1],
              pz[id_1], px[id_0], py[id_0], pz[id_0], box);
            Real3 r12 = MinImageDistanceVec(px[id_2], py[id_2],
              pz[id_2], px[id_1], py[id_1], pz[id_1], box);
            Real3 r20 = MinImageDistanceVec(px[id_0], py[id_0],
              pz[id_0], px[id_2], py[id_2], pz[id_2], box);

            //rbmd::Real sv10[3],sv21[3],sv02[3];
            Real3 sv10, sv21,sv02;
            sv10.x = shake_velocity_0.x - shake_velocity_1.x;
            sv10.y = shake_velocity_0.y - shake_velocity_1.y;
            sv10.z = shake_velocity_0.z - shake_velocity_1.z;

            sv21.x = shake_velocity_1.x - shake_velocity_2.x;
            sv21.y = shake_velocity_1.y - shake_velocity_2.y;
            sv21.z = shake_velocity_1.z - shake_velocity_2.z;

            sv02.x = shake_velocity_2.x - shake_velocity_0.x;
            sv02.y = shake_velocity_2.y - shake_velocity_0.y;
            sv02.z = shake_velocity_2.z - shake_velocity_0.z;

            rbmd::Real invmass0 = 1 / mass[atoms_type[id_0]];
            rbmd::Real invmass1 = 1 / mass[atoms_type[id_1]];
            rbmd::Real invmass2 = 1 / mass[atoms_type[id_2]];

            Real3 c, l;
            Real3 a_0, a_1, a_2;

            // setup matrix
            a_0.x = (invmass1 + invmass0) * Dot(r01, r01);
            a_0.y = -invmass1 * Dot(r01, r12);
            a_0.z = (-invmass0) * Dot(r01, r20);
            a_1.x = a_0.y;
            a_1.y = (invmass1 + invmass2) * Dot(r12, r12);
            a_1.z = -(invmass2) * Dot(r20, r12);
            a_2.x = a_0.z;
            a_2.y = a_1.z;
            a_2.z = (invmass0 + invmass2) * Dot(r20, r20);

            // sestup RHS
            c.x = -Dot(sv10, r01);
            c.y = -Dot(sv21, r12);
            c.z = -Dot(sv02, r20);

            Real3 ai_0, ai_1, ai_2;
            rbmd::Real determ, determinv = 0.0;

            // calculate the determinant of the matrix
            determ = a_0.x * a_1.y * a_2.z + a_0.y * a_1.z * a_2.x + a_0.z *
              a_1.x * a_2.y - a_0.x * a_1.z * a_2.y - a_0.y * a_1.x * a_2.z -
                a_0.z * a_1.y * a_2.x;

            // check if matrix is actually invertible
            if (ABS(determ) < 0.0001)
                printf(" Error: Rattle determinant = 0.0 ");

            // calculate the inverse 3x3 matrix: A^(-1) = (ai_jk)
            determinv = 1 / determ;
            ai_0.x =  determinv * (a_1.y * a_2.z - a_1.z * a_2.y);
            ai_0.y = -determinv * (a_0.y * a_2.z - a_0.z * a_2.y);
            ai_0.z =  determinv * (a_0.y * a_1.z - a_0.z * a_1.y);
            ai_1.x = -determinv * (a_1.x * a_2.z - a_1.z * a_2.x);
            ai_1.y =  determinv * (a_0.x * a_2.z - a_0.z * a_2.x);
            ai_1.z = -determinv * (a_0.x * a_1.z - a_0.z * a_1.x);
            ai_2.x =  determinv * (a_1.x * a_2.y - a_1.y * a_2.x);
            ai_2.y = -determinv * (a_0.x * a_2.y - a_0.y * a_2.x);
            ai_2.z =  determinv * (a_0.x * a_1.y - a_0.y * a_1.x);

            // calculate the solution:  (l01, l02, l12)^T = A^(-1) * c
            l.x = 0;
            l.y = 0;
            l.z = 0;

            l.x += ai_0.x * c.x;
            l.x += ai_0.y * c.y;
            l.x += ai_0.z * c.z;

            l.y += ai_1.x * c.x;
            l.y += ai_1.y * c.y;
            l.y += ai_1.z * c.z;

            l.z += ai_2.x * c.x;
            l.z += ai_2.y * c.y;
            l.z += ai_2.z * c.z;

            // [l01,l02,l12]^T = [lamda12,lamda23,lamda31]^T
            Real3 velocity_constraint_i0, velocity_constraint_i1, velocity_constraint_i2;
            velocity_constraint_i0.x = l.x * r01.x * invmass0 - l.z * r20.x * invmass0;
            velocity_constraint_i0.y = l.x * r01.y * invmass0 - l.z * r20.y * invmass0;
            velocity_constraint_i0.z = l.x * r01.z * invmass0 - l.z * r20.z * invmass0;

            velocity_constraint_i1.x = l.y * r12.x * invmass0 - l.x * r01.x * invmass0;
            velocity_constraint_i1.y = l.y * r12.y * invmass0 - l.x * r01.y * invmass0;
            velocity_constraint_i1.z = l.y * r12.z * invmass0 - l.x * r01.z * invmass0;

            velocity_constraint_i2.x = l.z * r20.x * invmass0 - l.y * r12.x * invmass0;
            velocity_constraint_i2.y = l.z * r20.y * invmass0 - l.y * r12.y * invmass0;
            velocity_constraint_i2.z = l.z * r20.z * invmass0 - l.y * r12.z * invmass0;

            shake_vx[id_0] = shake_velocity_0.x + velocity_constraint_i0.x;
            shake_vy[id_0] = shake_velocity_0.y + velocity_constraint_i0.y;
            shake_vz[id_0] = shake_velocity_0.z + velocity_constraint_i0.z;

            shake_vx[id_1] = shake_velocity_1.x + velocity_constraint_i1.x;
            shake_vy[id_1] = shake_velocity_1.y + velocity_constraint_i1.y;
            shake_vz[id_1] = shake_velocity_1.z + velocity_constraint_i1.z;

            shake_vx[id_2] = shake_velocity_2.x + velocity_constraint_i2.x;
            shake_vy[id_2] = shake_velocity_2.y + velocity_constraint_i2.y;
            shake_vz[id_2] = shake_velocity_2.z + velocity_constraint_i2.z;

          // rbmd::Real shake_vx_0 = shake_velocity_0.x + velocity_constraint_i0.x;
          // rbmd::Real shake_vy_0 = shake_velocity_0.y + velocity_constraint_i0.y;
          // rbmd::Real shake_vz_0 = shake_velocity_0.z + velocity_constraint_i0.z;
          //
          // rbmd::Real shake_vx_1 = shake_velocity_1.x + velocity_constraint_i1.x;
          // rbmd::Real shake_vy_1 = shake_velocity_1.y + velocity_constraint_i1.y;
          // rbmd::Real shake_vz_1 = shake_velocity_1.z + velocity_constraint_i1.z;
          //
          // rbmd::Real shake_vx_2 = shake_velocity_2.x + velocity_constraint_i2.x;
          // rbmd::Real shake_vy_2 = shake_velocity_2.y + velocity_constraint_i2.y;
          // rbmd::Real shake_vz_2 = shake_velocity_2.z + velocity_constraint_i2.z;
          // atomicAdd(&shake_vx[id_0], shake_vx_0);
          // atomicAdd(&shake_vy[id_0], shake_vy_0);
          // atomicAdd(&shake_vz[id_0], shake_vz_0);
          //
          // atomicAdd(&shake_vx[id_1], shake_vx_1);
          // atomicAdd(&shake_vy[id_1], shake_vy_1);
          // atomicAdd(&shake_vz[id_1], shake_vz_1);
          //
          // atomicAdd(&shake_vx[id_2], shake_vx_2);
          // atomicAdd(&shake_vy[id_2], shake_vy_2);
          // atomicAdd(&shake_vz[id_2], shake_vz_2);

            Real3 a00, a01, a02, a10, a11, a12, a20, a21, a22, cal_v12t, cal_v23t, cal_v31t;
            a00.x=(invmass1 + invmass0) * r01.x;
            a00.y=(invmass1 + invmass0) * r01.y;
            a00.z=(invmass1 + invmass0) * r01.z;

            a01.x = -invmass1 * r12.x;
            a01.y = -invmass1 * r12.y;
            a01.z = -invmass1 * r12.z;

            a02.x = (-invmass0) * r20.x;
            a02.y = (-invmass0) * r20.y;
            a02.z = (-invmass0) * r20.z;

            a10.x = -invmass1 * r01.x;
            a10.y = -invmass1 * r01.y;
            a10.z = -invmass1 * r01.z;

            a11.x = (invmass1 + invmass2) * r12.x;
            a11.y = (invmass1 + invmass2) * r12.y;
            a11.z = (invmass1 + invmass2) * r12.z;

            a12.x = -(invmass2) * r20.x;
            a12.y = -(invmass2) * r20.y;
            a12.z = -(invmass2) * r20.z;

            a20.x = (-invmass0) * r01.x;
            a20.y = (-invmass0) * r01.y;
            a20.z = (-invmass0) * r01.z;

            a21.x = -(invmass2) * r12.x;
            a21.y = -(invmass2) * r12.y;
            a21.z = -(invmass2) * r12.z;

            a22.x = (invmass0 + invmass2) * r20.x;
            a22.y = (invmass0 + invmass2) * r20.y;
            a22.z = (invmass0 + invmass2) * r20.z;

            cal_v12t.x = sv10.x + (a00.x * l.x + a01.x * l.y + a02.x * l.z);
            cal_v12t.y = sv10.y + (a00.y * l.x + a01.y * l.y + a02.y * l.z);
            cal_v12t.z = sv10.z + (a00.z * l.x + a01.z * l.y + a02.z * l.z);

            cal_v23t.x = sv21.x + (a10.x * l.x + a11.x * l.y + a12.x * l.z);
            cal_v23t.y = sv21.y + (a10.y * l.x + a11.y * l.y + a12.y * l.z);
            cal_v23t.z = sv21.z + (a10.z * l.x + a11.z * l.y + a12.z * l.z);

            cal_v31t.x = sv02.x + (a20.x * l.x + a21.x * l.y + a22.x * l.z);
            cal_v31t.y = sv02.y + (a20.y * l.x + a21.y * l.y + a22.y * l.z);
            cal_v31t.z = sv02.z + (a20.z * l.x + a21.z * l.y + a22.z * l.z);


            rbmd::Real cal_dv1 = Dot(r01, cal_v12t);
            rbmd::Real cal_dv2 = Dot(r20, cal_v31t);
            rbmd::Real cal_dv12 = Dot(r12, cal_v23t);

            Real3 velocity01, velocity20, velocity12;
            velocity01.x = shake_vx[id_1] - shake_vx[id_0];
            velocity01.y = shake_vy[id_1] - shake_vy[id_0];
            velocity01.z = shake_vz[id_1] - shake_vz[id_0];

            velocity20.x = shake_vx[id_0] - shake_vx[id_2];
            velocity20.y = shake_vy[id_0] - shake_vy[id_2];
            velocity20.z = shake_vz[id_0] - shake_vz[id_2];

            velocity12.x = shake_vx[id_2] - shake_vx[id_1];
            velocity12.y = shake_vy[id_2] - shake_vy[id_1];
            velocity12.z = shake_vz[id_2] - shake_vz[id_1];

            rbmd::Real dv1 = Dot(r01, velocity01);
            rbmd::Real dv2 = Dot(r20, velocity20);
            rbmd::Real dv12 = Dot(r12, velocity12);

            if (ABS(dv1) > 0.1 || ABS(dv2) > 0.1 || ABS(dv12) > 0.1)
            {
                printf("i0 = %d, i1 = %d, i2 = %d\n", id_0, id_1, id_2);
                printf("dv1 = %f, dv2 = %f, dv12 = %f\n", dv1, dv2, dv12);
                printf("cal_dv1 = %f, cal_dv2 = %f, cal_dv12 = %f\n", cal_dv1, cal_dv2, cal_dv12);
                printf("velocity_i0 = [%f,%f,%f], velocity_i1 = [%f,%f,%f], velocity_i2 = [%f,%f,%f]\n",
                       shake_vx[id_0], shake_vy[id_0], shake_vz[id_0],
                       shake_vx[id_1], shake_vy[id_1], shake_vz[id_1],
                       shake_vx[id_2], shake_vy[id_2], shake_vz[id_2]);
                printf("\n");
            }
        }
    }

void ShakeAOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_angle,
                                              const rbmd::Real dt,
                                              const rbmd::Real fmt2v,
                                              Box* box,
                                              const rbmd::Id* atom_id_to_idx,
                                              const rbmd::Real* mass,
                                              const rbmd::Id* atoms_type,
                                              const Id3* angle_id_vec,
                                              rbmd::Real* shake_px,
                                              rbmd::Real* shake_py,
                                              rbmd::Real* shake_pz,
                                              rbmd::Real* shake_vx,
                                              rbmd::Real* shake_vy,
                                              rbmd::Real* shake_vz,
                                              const rbmd::Real* fx,
                                              const rbmd::Real* fy,
                                              const rbmd::Real* fz,
                                              rbmd::Id* flag_px,
                                              rbmd::Id* flag_py,
                                              rbmd::Id* flag_pz) 
{
  unsigned int blocks_per_grid = (num_angle + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(ShakeA<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_angle, dt, fmt2v, box, atom_id_to_idx,mass, atoms_type,
      angle_id_vec,shake_px, shake_py, shake_pz, shake_vx, shake_vy,
      shake_vz, fx, fy, fz, flag_px, flag_py, flag_pz));
}

    void ShakeBOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_angle,
                                                  const rbmd::Real dt,
                                                  const rbmd::Real fmt2v,
                                                  Box* box,
                                                  const rbmd::Id* atom_id_to_idx,
                                                  const rbmd::Real* mass,
                                                  const rbmd::Id* atoms_type,
                                                  const Id3* angle_id_vec,
                                                  const rbmd::Real* px,
                                                  const rbmd::Real* py,
                                                  const rbmd::Real* pz,
                                                  rbmd::Real* shake_vx,
                                                  rbmd::Real* shake_vy,
                                                  rbmd::Real* shake_vz,
                                                  const rbmd::Real* fx,
                                                  const rbmd::Real* fy,
                                                  const rbmd::Real* fz)
    {
        unsigned int blocks_per_grid = (num_angle + BLOCK_SIZE - 1) / BLOCK_SIZE;
        CHECK_KERNEL(ShakeB<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>
          (num_angle, dt, fmt2v, box, atom_id_to_idx, mass, atoms_type,
           angle_id_vec,px, py, pz, shake_vx, shake_vy, shake_vz, fx, fy, fz));
    }
}  // namespace op
