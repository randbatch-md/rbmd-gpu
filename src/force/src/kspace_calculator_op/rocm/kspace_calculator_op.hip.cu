#include "../common/rbmd_define.h"
#include "force/src/kspace_calculator_op/kspace_calculator_op.h"

#include "model/box.h"

static const double MY_PIS=1.77245385090551602729;
static const double MY_PI2=1.57079632679489661923;
static const int  OffSet_warpSize = 32;
namespace op {

  //  Warp-level
  __device__ __forceinline__ rbmd::Real warpReduceSum(rbmd::Real val) {
    // #pragma unroll
    for (int offset = OffSet_warpSize / 2; offset > 0; offset /= 2) {
      val += SHFL_DOWN_SYNC(0xFFFFFFFF, val, offset,OffSet_warpSize);
    }
    return val;
  }

  //  Block-level
  __device__ __forceinline__ rbmd::Real blockReduceSum(rbmd::Real val) {

    __shared__ rbmd::Real shared[OffSet_warpSize];

    int lane = threadIdx.x % OffSet_warpSize;
    int wid = threadIdx.x / OffSet_warpSize;

    // 1.
    val = warpReduceSum(val);

    // 2.
    if (lane == 0) {
      shared[wid] = val;
    }
    __syncthreads(); //

    // 3.
    val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0.0;
    if (wid == 0) {
      val = warpReduceSum(val);
    }

    return val;
  }

  __device__ inline void block_reduce_sum(rbmd::Real& sum) {
    __shared__ rbmd::Real s_data[BLOCK_SIZE];
    s_data[threadIdx.x] = sum;
    __syncthreads();

    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
      if (threadIdx.x < offset) {
        s_data[threadIdx.x] += s_data[threadIdx.x + offset];
      }
      __syncthreads();
    }
    sum = s_data[0]; // Result is in thread 0
  }

  // EwaldForce
  __device__ void EwaldForce( Box box, const rbmd::Real alpha, const Int3 M,
                             const rbmd::Real qqr2e, const rbmd::Real rhok_real_i,
                             const rbmd::Real rhok_imag_i,
                             const rbmd::Real charge, const rbmd::Real px,
                             const rbmd::Real py, const rbmd::Real pz,
                             rbmd::Real& force_ewald_single,
                             rbmd::Real& force_ewald_x, rbmd::Real& force_ewald_y,
                             rbmd::Real& force_ewald_z) {
    rbmd::Real force_ewald;
    rbmd::Real volume = box._length[0] * box._length[1] * box._length[2];
    Real3 K = make_Real3(2 * M_PI * M.x / box._length[0],
                         2 * M_PI * M.y / box._length[1],
                         2 * M_PI * M.z / box._length[2]);

    rbmd::Real range_K_2 = K.x * K.x + K.y * K.y + K.z * K.z;
    rbmd::Real dot_product = (-K.x)* px + (-K.y) * py + (-K.z) * pz;
    rbmd::Real alpha_inv = 1 / alpha;

    rbmd::Real factor_a = -4 * M_PI * charge; //+
    rbmd::Real factor_b = EXP(-0.25 * range_K_2 * alpha_inv);
    rbmd::Real factor_c = COS(dot_product) * rhok_imag_i;
    rbmd::Real factor_d = SIN(dot_product) * rhok_real_i;

     force_ewald =
         factor_a / (volume * range_K_2) * factor_b * (factor_c + factor_d);
    force_ewald *= qqr2e;
    force_ewald_x = force_ewald * K.x;
    force_ewald_y = force_ewald * K.y;
    force_ewald_z = force_ewald * K.z;
    //
    force_ewald_single = force_ewald;
  }

  // RBEForce
  __device__ void RBEForce( Box box, const Real3 M, const rbmd::Real qqr2e,
                           const rbmd::Real rhok_real_i,
                           const rbmd::Real rhok_imag_i, const rbmd::Real charge,
                           const rbmd::Real px, const rbmd::Real py,
                           const rbmd::Real pz, rbmd::Real& force_rbe_single,
                           rbmd::Real& force_rbe_x,rbmd::Real& force_rbe_y,
                           rbmd::Real& force_rbe_z) {
    rbmd::Real force_rbe;
    rbmd::Real volume = box._length[0] * box._length[1] * box._length[2];
    Real3 K = make_Real3(2 * M_PI * M.x / box._length[0],
                         2 * M_PI * M.y / box._length[1],
                         2 * M_PI * M.z / box._length[2]);

    rbmd::Real range_K_2 = K.x * K.x + K.y * K.y + K.z * K.z;
    rbmd::Real dot_product = K.x * px + K.y * py + K.z * pz;

    rbmd::Real factor_a = -4 * M_PI * charge;
    rbmd::Real factor_b = COS(dot_product) * rhok_imag_i;
    rbmd::Real factor_c = SIN(dot_product) * rhok_real_i;

    force_rbe = (factor_a / (volume * range_K_2)) * (factor_b - factor_c);
    force_rbe *= qqr2e;
    force_rbe_x = force_rbe * K.x;
    force_rbe_y = force_rbe * K.y;
    force_rbe_z = force_rbe * K.z;
    //
    force_rbe_single = force_rbe;
  }

  __device__ void ComputeS( Box box, const rbmd::Real alpha, rbmd::Real& S) {
    Real3 H{0.0, 0.0, 0.0};
    for (rbmd::Id i = 0; i < 3; ++i) {
      const rbmd::Real factor = -(alpha * box._length[i] * box._length[i]);

      for (rbmd::Id m = -10; m <= 10; m++) {
        rbmd::Real expx = m * m * factor;
         REAL_DATA(H)[i] += EXP(expx);
      }
      REAL_DATA(H)[i] *= SQRT(-(factor) / M_PI);
    }

    rbmd::Real factor_3 = REAL_DATA(H)[0] * REAL_DATA(H)[1] * REAL_DATA(H)[2];
    S = factor_3 - 1;
  }

  template <typename Func>
  __device__ void ExecuteOnKmax(const Int3& k_maxconst, Func& function) {
    rbmd::Id indexEwald = 0;
    for (rbmd::Id i = -INT_DATA(k_maxconst)[0]; i <= INT_DATA(k_maxconst)[0]; i++) {
      for (rbmd::Id j = -INT_DATA(k_maxconst)[1]; j <= INT_DATA(k_maxconst)[1]; j++) {
        for (rbmd::Id k = -INT_DATA(k_maxconst)[2]; k <= INT_DATA(k_maxconst)[2]; k++) {
          if(!(i == 0 && j == 0 && k == 0)){
            Int3 M = make_Int3(i, j, k); //fix
            function(M, indexEwald);
            indexEwald++;
          }
        }
      }
    }
  }



// StructureFactor
__global__ void ComputeChargeStructureFactor(
    const rbmd::Id num_atoms, const Real3 K, const rbmd::Real* charge,
    const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* density_real, rbmd::Real* density_imag) {
    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_atoms) {
      rbmd::Real local_charge = charge[tid1];
      rbmd::Real dot_product = K.x * px[tid1] + K.y * py[tid1] + K.z * pz[tid1];

      density_real[tid1] = local_charge * COS(dot_product);
      density_imag[tid1] = local_charge * SIN(dot_product);
    }
  }

// Charge Structure  Factor on Pnumber
__global__ void ComputePnumberChargeStructureFactor(
    Box   box, const rbmd::Id num_atoms, const rbmd::Id p_number,
    const rbmd::Real* __restrict__ charge, const rbmd::Real* __restrict__ p_sample_x,
    const rbmd::Real* __restrict__ p_sample_y, const rbmd::Real* __restrict__ p_sample_z,
    const rbmd::Real* __restrict__ px, const rbmd::Real* __restrict__ py, const rbmd::Real* __restrict__ pz,
    rbmd::Real* __restrict__ density_real, rbmd::Real* __restrict__ density_imag) {
    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    __shared__ rbmd::Real shared_px[BLOCK_SIZE];
    __shared__ rbmd::Real shared_py[BLOCK_SIZE];
    __shared__ rbmd::Real shared_pz[BLOCK_SIZE];
    __shared__ rbmd::Real shared_charge[BLOCK_SIZE];

    if (tid1 < num_atoms) {
      if (threadIdx.x < BLOCK_SIZE) {
        shared_px[threadIdx.x] = __ldg(&px[tid1]);
        shared_py[threadIdx.x] = __ldg(&py[tid1]);
        shared_pz[threadIdx.x] = __ldg(&pz[tid1]);
        shared_charge[threadIdx.x] = __ldg(&charge[tid1]);
      }
      __syncthreads();
      rbmd::Real chargei = shared_charge[threadIdx.x];
      rbmd::Real p_x = shared_px[threadIdx.x];
      rbmd::Real p_y = shared_py[threadIdx.x];
      rbmd::Real p_z = shared_pz[threadIdx.x];

      for (rbmd::Id i = 0; i < p_number; i++) {
        rbmd::Id index = tid1 + i * num_atoms;

        rbmd::Real k_x = __ldg(&p_sample_x[i]);
        rbmd::Real k_y = __ldg(&p_sample_y[i]);
        rbmd::Real k_z = __ldg(&p_sample_z[i]);
        k_x = 2 * M_PI * k_x / box._length[0];
        k_y = 2 * M_PI * k_y / box._length[1];
        k_z = 2 * M_PI * k_z / box._length[2];

        rbmd::Real dot_product = k_x * p_x + k_y * p_y + k_z * p_z;
        density_real[index] = chargei * COS(dot_product);
        density_imag[index] = chargei * SIN(dot_product);
      }
    }
  }

__global__ void ComputePnumberChargeStructureFactorOpt(
    Box box, const rbmd::Id num_atoms, const rbmd::Id p_number,
    const rbmd::Real* __restrict__ charge,
    const rbmd::Real* __restrict__ p_sample_x,
    const rbmd::Real* __restrict__ p_sample_y,
    const rbmd::Real* __restrict__ p_sample_z,
    const rbmd::Real* __restrict__ px,
    const rbmd::Real* __restrict__ py,
    const rbmd::Real* __restrict__ pz,
    rbmd::Real* __restrict__ rhok_real,
    rbmd::Real* __restrict__ rhok_image)
{
  const rbmd::Id P_index = blockIdx.x;
  if (P_index >= p_number) return;

  //
  rbmd::Real m_x = __ldg(&p_sample_x[P_index]);
  rbmd::Real m_y = __ldg(&p_sample_y[P_index]);
  rbmd::Real m_z = __ldg(&p_sample_z[P_index]);

  rbmd::Real k_x = 2 * M_PI * m_x / box._length[0];
  rbmd::Real k_y = 2 * M_PI * m_y / box._length[1];
  rbmd::Real k_z = 2 * M_PI * m_z / box._length[2];

  rbmd::Real local_real_sum = 0.0;
  rbmd::Real local_image_sum = 0.0;

  // Grid-stride
  for (rbmd::Id N_index = threadIdx.x; N_index < num_atoms; N_index += blockDim.x) {
    if (N_index >= num_atoms) continue;

    rbmd::Real chargei = __ldg(&charge[N_index]);
    rbmd::Real p_x = __ldg(&px[N_index]);
    rbmd::Real p_y = __ldg(&py[N_index]);
    rbmd::Real p_z = __ldg(&pz[N_index]);

    rbmd::Real dot_product = k_x * p_x + k_y * p_y + k_z * p_z;

    local_real_sum += chargei * COS(dot_product);
    local_image_sum += chargei * SIN(dot_product);
  }

  //
  rbmd::Real total_real = blockReduceSum(local_real_sum);
  rbmd::Real total_image = blockReduceSum(local_image_sum);

  //
  if (threadIdx.x == 0) {
    rhok_real[P_index] = total_real;
    rhok_image[P_index] = total_image;
  }
}


__global__ void EwaldForceFix(const rbmd::Id num_atoms,const rbmd::Id kcount,
  const rbmd::Id k_index,const rbmd::Real qqr2e,Int3 kmax_vec3D,
  const rbmd::Real* eg,const rbmd::Real* cs,const rbmd::Real* sn,
  const rbmd::Real* charge,const rbmd::Real* qfactor_real,
  const rbmd::Real* qfactor_image,rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz)
{
    rbmd::Real sum_fx = 0;
    rbmd::Real sum_fy = 0;
    rbmd::Real sum_fz = 0;

    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_atoms)
    {
      rbmd::Real charge_i = charge[tid1];
      rbmd::Real cos_couple_yz,sin_couple_yz;
      rbmd::Real position_phase_real,position_phase_image;
      rbmd::Real partial;

      rbmd::Id kx_index_0 = kmax_vec3D.x  * (3*num_atoms) + 0 * num_atoms + tid1;
      rbmd::Id ky_index_1 = kmax_vec3D.y  * (3*num_atoms) + 1 * num_atoms + tid1;
      rbmd::Id kz_index_2 = kmax_vec3D.z  * (3*num_atoms) + 2 * num_atoms + tid1;

      cos_couple_yz = cs[ky_index_1] *cs[kz_index_2] -sn[ky_index_1] *sn[kz_index_2];
      sin_couple_yz = sn[ky_index_1] *cs[kz_index_2] +cs[ky_index_1] *sn[kz_index_2];

      position_phase_real  = cs[kx_index_0] *cos_couple_yz - sn[kx_index_0] *sin_couple_yz;
      position_phase_image = sn[kx_index_0] *cos_couple_yz + cs[kx_index_0] * sin_couple_yz;

       partial = position_phase_real * qfactor_real[k_index] -
                  position_phase_image * qfactor_image[k_index];
       sum_fx += partial * eg[k_index * kcount * 3 + 0];
       sum_fy += partial * eg[k_index * kcount * 3 + 1];
       sum_fz += partial * eg[k_index * kcount * 3 + 2];
      //printf("force: %f %f  %f\n", sum_fx,sum_fy,sum_fz);
      fx[tid1] = qqr2e * charge_i *sum_fx;
      fy[tid1] = qqr2e * charge_i *sum_fy;
      fz[tid1] = qqr2e * charge_i *sum_fz;
    }
}


  // EwaldForce
  __global__ void ComputeEwaldForce(
       Box box, const rbmd::Id num_atoms, const Int3 Kmax,
      const rbmd::Real alpha, const rbmd::Real qqr2e,
      const rbmd::Real* real_array, const rbmd::Real* imag_array,
      const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py,
      const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
      rbmd::Real* flat_virial) {
    rbmd::Real sum_fx = 0;
    rbmd::Real sum_fy = 0;
    rbmd::Real sum_fz = 0;
    //virial init
    rbmd::Real sum_virial[6];
    for (int i = 0; i < 6; ++i)
    {
      sum_virial[i] = 0.0;
    }

    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_atoms) {
      rbmd::Real p_x = px[tid1];
      rbmd::Real p_y = py[tid1];
      rbmd::Real p_z = pz[tid1];
      rbmd::Real charge_i = charge[tid1];

      auto function = [&](const Int3& M, const rbmd::Id& indexEwald) {
        const rbmd::Real rhok_real_i = real_array[indexEwald];
        const rbmd::Real rhok_imag_i = imag_array[indexEwald];

        rbmd::Real force_Ewald_x, force_Ewald_y, force_Ewald_z;
        rbmd::Real force_Ewald_single;
        EwaldForce(box, alpha, M, qqr2e, rhok_real_i, rhok_imag_i, charge_i, p_x,
                   p_y, p_z, force_Ewald_single,force_Ewald_x, force_Ewald_y, force_Ewald_z);

        sum_fx += force_Ewald_x;
        sum_fy += force_Ewald_y;
        sum_fz += force_Ewald_z;

        //compute virial_ewald
        //force_Ewald_single = -0.5*force_Ewald_single;
        Real3 K;
        K.x = 2.0 * M_PI * M.x / box._length[0];
        K.y = 2.0 * M_PI * M.y / box._length[1];
        K.z = 2.0 * M_PI * M.z / box._length[2];
        sum_virial[0] += (p_x * K.x + p_x * K.x) * force_Ewald_single; //_xx
        sum_virial[1] += (p_y * K.y + p_y * K.y) * force_Ewald_single; // yy
        sum_virial[2] += (p_z * K.z + p_z * K.z) * force_Ewald_single; // zz
        sum_virial[3] += (p_x * K.y + p_y * K.x) * force_Ewald_single; // xy
        sum_virial[4] += (p_x * K.z + p_z * K.x) * force_Ewald_single; // xz
        sum_virial[5] += (p_y * K.z + p_z * K.y) * force_Ewald_single; // yz
      };
      ExecuteOnKmax(Kmax, function);

      fx[tid1] = sum_fx;
      fy[tid1] = sum_fy;
      fz[tid1] = sum_fz;
      //
      for(int i =0;i<6;++i) {
        flat_virial[  i * num_atoms + tid1] = sum_virial[i];
      }
    }
  }



__global__ void ComputeSqCharge(const rbmd::Id num_atoms,
                                const rbmd::Real* charge,
                                rbmd::Real* sq_charge) {
    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_atoms) {
      rbmd::Real chargei = charge[tid1];
      sq_charge[tid1] = chargei * chargei;
    }
  }

// index
__global__ void GenerateIndexArray(const rbmd::Id num_atoms,
                                   const rbmd::Id RBE_P,
                                   rbmd::Id* psample_key) {
    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_atoms * RBE_P) {
      psample_key[tid1] = tid1 / num_atoms;
    }
  }

__global__ void ComputeRBEVirial(
    const Box box,
    const rbmd::Id P,
    const rbmd::Real alpha,
    const rbmd::Real qqrd2e,
    const rbmd::Real qsqsum,
    const rbmd::Real qsum,
    const rbmd::Real* real_array,
    const rbmd::Real* imag_array,
    const rbmd::Real* p_sample_x,
    const rbmd::Real* p_sample_y,
    const rbmd::Real* p_sample_z,
    rbmd::Real* virial_tensor,
    rbmd::Real* energy
) {

    //
    //__shared__ rbmd::Real shared_virial[6];
    __shared__ rbmd::Real shared_energy;
    __shared__ rbmd::Real shared_virial[6][BLOCK_SIZE];
    //
  // 初始化共享内存
  if (threadIdx.x < 6) {
    for (int i = 0; i < BLOCK_SIZE; i++) {
      shared_virial[threadIdx.x][i] = 0.0;
    }
  }
  __syncthreads();

    // if (threadIdx.x < 6) {
    //     shared_virial[threadIdx.x] = 0.0;
    // }
    if (threadIdx.x == 0) {
      shared_energy = 0.0;
    }
    __syncthreads();

    //
    rbmd::Real S, V;
    ComputeS(box, alpha, S);
    V = box._length[0] * box._length[1] * box._length[2];

    //
    const rbmd::Real energy_coef = (2 * M_PI * S) / (P * V);
    const rbmd::Real virial_coef = - (S / P) * M_PI / V;

    for (rbmd::Id i = threadIdx.x; i < P; i += blockDim.x){
        //
        const rbmd::Real m_x = __ldg(&p_sample_x[i]);
        const rbmd::Real m_y = __ldg(&p_sample_y[i]);
        const rbmd::Real m_z = __ldg(&p_sample_z[i]);
        //printf("test---  %f %f  %f ", m_x ,m_y,m_z);

        // (2πm/L)
        const rbmd::Real k_x = 2 * M_PI * m_x / box._length[0];
        const rbmd::Real k_y = 2 * M_PI * m_y / box._length[1];
        const rbmd::Real k_z = 2 * M_PI * m_z / box._length[2];

        // |k|²
        const rbmd::Real k2 = k_x * k_x + k_y * k_y + k_z * k_z;

        //
        const rbmd::Real rho_real = __ldg(&real_array[i]);
        const rbmd::Real rho_imag = __ldg(&imag_array[i]);

        // |ρ(k)|²
        const rbmd::Real rho2 = rho_real * rho_real + rho_imag * rho_imag;

      //
        const rbmd::Real energy_contribution = energy_coef * rho2 / k2;
        atomicAdd(&shared_energy, energy_contribution);

        //printf("test---  %f %f  %f %f", k2 ,rho2, S,V);
        const rbmd::Real base = virial_coef * rho2 / k2;


        // (1/4α² + 1/|k|²)
        const rbmd::Real inner_coef = (1.0 / (4 * alpha)) + (1.0 / k2);

        //  (δ_βγ - 2k_βk_γ inner_coef)
        const rbmd::Real delta_term = 1.0; // δ_βγ
        const rbmd::Real kterm_xx = 2 * k_x * k_x * inner_coef;
        const rbmd::Real kterm_yy = 2 * k_y * k_y * inner_coef;
        const rbmd::Real kterm_zz = 2 * k_z * k_z * inner_coef;
        const rbmd::Real kterm_xy = 2 * k_x * k_y * inner_coef;
        const rbmd::Real kterm_xz = 2 * k_x * k_z * inner_coef;
        const rbmd::Real kterm_yz = 2 * k_y * k_z * inner_coef;

        //
        // atomicAdd(&shared_virial[0], base * (delta_term - kterm_xx)); // xx
        // atomicAdd(&shared_virial[1], base * (delta_term - kterm_yy)); // yy
        // atomicAdd(&shared_virial[2], base * (delta_term - kterm_zz)); // zz
        // atomicAdd(&shared_virial[3], base * (-kterm_xy));              // xy
        // atomicAdd(&shared_virial[4], base * (-kterm_xz));              // xz
        // atomicAdd(&shared_virial[5], base * (-kterm_yz));              // yz

      shared_virial[0][threadIdx.x] += base * (delta_term - kterm_xx);         // xx
      shared_virial[1][threadIdx.x] += base * (delta_term - kterm_yy);       //yy
      shared_virial[2][threadIdx.x] += base * (delta_term - kterm_zz);  // zz
      shared_virial[3][threadIdx.x] += base * (-kterm_xy);        // xy
      shared_virial[4][threadIdx.x] += base * (-kterm_xz);        // xz
      shared_virial[5][threadIdx.x] += base * (-kterm_yz);        // yz
    }
    __syncthreads();


    // 树状归约
    for (int s = BLOCK_SIZE / 2; s > 0; s >>= 1) {
      if (threadIdx.x < s) {
        for (int j = 0; j < 6; j++) {
          shared_virial[j][threadIdx.x] += shared_virial[j][threadIdx.x + s];
        }
      }
      __syncthreads();
    }

    //
    if (threadIdx.x == 0) {
      // 1.
      shared_energy -= SQRT(alpha) * qsqsum / MY_PIS + MY_PI2 * qsum * qsum / (alpha * V);
      *energy = shared_energy * qqrd2e;

      // 2.
      for (int j = 0; j < 6; j++) {
       // virial_tensor[j] = shared_virial[j] * qqrd2e;
        virial_tensor[j] = shared_virial[j][0] * qqrd2e;
      }
    }
}

__global__ void ComputeRBEVirial0(
    const Box box,
    const rbmd::Id P,
    const rbmd::Real alpha,
    const rbmd::Real qqrd2e,
    const rbmd::Real qsqsum,
    const rbmd::Real qsum,
    const rbmd::Real* real_array,
    const rbmd::Real* imag_array,
    const rbmd::Real* p_sample_x,
    const rbmd::Real* p_sample_y,
    const rbmd::Real* p_sample_z,
    rbmd::Real* virial_tensor,
    rbmd::Real* energy
) {

    //
    __shared__ rbmd::Real shared_virial[6];
    __shared__ rbmd::Real shared_energy;

    if (threadIdx.x < 6) {
        shared_virial[threadIdx.x] = 0.0;
    }
    if (threadIdx.x == 0) {
      shared_energy = 0.0;
    }
    __syncthreads();

    //
    rbmd::Real S, V;
    ComputeS(box, alpha, S);
    V = box._length[0] * box._length[1] * box._length[2];

    //
    const rbmd::Real energy_coef = (2 * M_PI * S) / (P * V);
    const rbmd::Real virial_coef = - (S / P) * M_PI / V;

    for (rbmd::Id i = threadIdx.x; i < P; i += blockDim.x){
        //
        const rbmd::Real m_x = __ldg(&p_sample_x[i]);
        const rbmd::Real m_y = __ldg(&p_sample_y[i]);
        const rbmd::Real m_z = __ldg(&p_sample_z[i]);
        //printf("test---  %f %f  %f ", m_x ,m_y,m_z);

        // (2πm/L)
        const rbmd::Real k_x = 2 * M_PI * m_x / box._length[0];
        const rbmd::Real k_y = 2 * M_PI * m_y / box._length[1];
        const rbmd::Real k_z = 2 * M_PI * m_z / box._length[2];

        // |k|²
        const rbmd::Real k2 = k_x * k_x + k_y * k_y + k_z * k_z;

        //
        const rbmd::Real rho_real = __ldg(&real_array[i]);
        const rbmd::Real rho_imag = __ldg(&imag_array[i]);

        // |ρ(k)|²
        const rbmd::Real rho2 = rho_real * rho_real + rho_imag * rho_imag;

      //
        const rbmd::Real energy_contribution = energy_coef * rho2 / k2;
        atomicAdd(&shared_energy, energy_contribution);

        //printf("test---  %f %f  %f %f", k2 ,rho2, S,V);
        const rbmd::Real base = virial_coef * rho2 / k2;


        // (1/4α² + 1/|k|²)
        const rbmd::Real inner_coef = (1.0 / (4 * alpha)) + (1.0 / k2);

        //  (δ_βγ - 2k_βk_γ inner_coef)
        const rbmd::Real delta_term = 1.0; // δ_βγ
        const rbmd::Real kterm_xx = 2 * k_x * k_x * inner_coef;
        const rbmd::Real kterm_yy = 2 * k_y * k_y * inner_coef;
        const rbmd::Real kterm_zz = 2 * k_z * k_z * inner_coef;
        const rbmd::Real kterm_xy = 2 * k_x * k_y * inner_coef;
        const rbmd::Real kterm_xz = 2 * k_x * k_z * inner_coef;
        const rbmd::Real kterm_yz = 2 * k_y * k_z * inner_coef;

        //
        atomicAdd(&shared_virial[0], base * (delta_term - kterm_xx)); // xx
        atomicAdd(&shared_virial[1], base * (delta_term - kterm_yy)); // yy
        atomicAdd(&shared_virial[2], base * (delta_term - kterm_zz)); // zz
        atomicAdd(&shared_virial[3], base * (-kterm_xy));              // xy
        atomicAdd(&shared_virial[4], base * (-kterm_xz));              // xz
        atomicAdd(&shared_virial[5], base * (-kterm_yz));              // yz
    }
    __syncthreads();

    //
    if (threadIdx.x == 0) {
      // 1.
      shared_energy -= SQRT(alpha) * qsqsum / MY_PIS + MY_PI2 * qsum * qsum / (alpha * V);
      *energy = shared_energy * qqrd2e;

      // 2.
      for (int j = 0; j < 6; j++) {
       virial_tensor[j] = shared_virial[j] * qqrd2e;
      }
    }
}


  // RBEForce
  __global__ void ComputeRBEForce(
      Box   box, const rbmd::Id num_atoms, const rbmd::Id p_number,
      const rbmd::Real alpha, const rbmd::Real qqr2e,
      const rbmd::Real* __restrict__ real_array, const rbmd::Real* __restrict__ imag_array,
      const rbmd::Real* __restrict__ charge, const rbmd::Real* __restrict__ p_sample_x,
      const rbmd::Real* __restrict__ p_sample_y, const rbmd::Real* __restrict__ p_sample_z,
      const rbmd::Real* __restrict__ px, const rbmd::Real* __restrict__ py, const rbmd::Real* __restrict__ pz,
      rbmd::Real* __restrict__ fx, rbmd::Real* __restrict__ fy, rbmd::Real* __restrict__ fz,
      rbmd::Real* flat_virial) {
    rbmd::Real sum_fx = 0;
    rbmd::Real sum_fy = 0;
    rbmd::Real sum_fz = 0;

    __shared__ rbmd::Real shared_px[BLOCK_SIZE];
    __shared__ rbmd::Real shared_py[BLOCK_SIZE];
    __shared__ rbmd::Real shared_pz[BLOCK_SIZE];
    __shared__ rbmd::Real shared_charge[BLOCK_SIZE];

    //virial init
    rbmd::Real sum_virial[6];
    for (int i = 0; i < 6; ++i)
    {
      sum_virial[i] = 0.0;
    }

    unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid1 < num_atoms) {
      if (threadIdx.x < BLOCK_SIZE) {
        shared_px[threadIdx.x] = __ldg(&px[tid1]);
        shared_py[threadIdx.x] = __ldg(&py[tid1]);
        shared_pz[threadIdx.x] = __ldg(&pz[tid1]);
        shared_charge[threadIdx.x] = __ldg(&charge[tid1]);
      }
      __syncthreads();
      rbmd::Real p_x =  shared_px[threadIdx.x];
      rbmd::Real p_y = shared_py[threadIdx.x];
      rbmd::Real p_z = shared_pz[threadIdx.x];
      rbmd::Real charge_i = shared_charge[threadIdx.x];
      //
      for (rbmd::Id i = 0; i < p_number; i++) {
        const Real3 M = make_Real3(__ldg(&p_sample_x[i]), __ldg(&p_sample_y[i]), __ldg(&p_sample_z[i]));

        const rbmd::Real rhok_real_i = __ldg(&real_array[i]);
        const rbmd::Real rhok_imag_i = __ldg(&imag_array[i]);
        //printf(" test--  %f %f",rhok_real_i,rhok_imag_i);

        rbmd::Real force_rbe_x, force_rbe_y, force_rbe_z;
        rbmd::Real force_rbe_single;
        RBEForce(box, M, qqr2e, rhok_real_i, rhok_imag_i, charge_i, p_x, p_y, p_z,
                 force_rbe_single,force_rbe_x, force_rbe_y, force_rbe_z);

        sum_fx += force_rbe_x;
        sum_fy += force_rbe_y;
        sum_fz += force_rbe_z;

        //compute virial_rbe
        Real3 K;
        K.x = 2.0 * M_PI * M.x / box._length[0];
        K.y = 2.0 * M_PI * M.y / box._length[1];
        K.z = 2.0 * M_PI * M.z / box._length[2];
        sum_virial[0] += (p_x * K.x + p_x * K.x) * force_rbe_single; //_xx
        sum_virial[1] += (p_y * K.y + p_y * K.y) * force_rbe_single; // yy
        sum_virial[2] += (p_z * K.z + p_z * K.z) * force_rbe_single; // zz
        sum_virial[3] += (p_x * K.y + p_y * K.x) * force_rbe_single; // xy
        sum_virial[4] += (p_x * K.z + p_z * K.x) * force_rbe_single; // xz
        sum_virial[5] += (p_y * K.z + p_z * K.y) * force_rbe_single; // yz
      }
      //
      rbmd::Real sum_gauss;
      ComputeS(box, alpha, sum_gauss);
      // printf("--------test---sum_gauss:%f\n",sum_gauss);

      sum_fx = sum_fx * sum_gauss / p_number;
      sum_fy = sum_fy * sum_gauss / p_number;
      sum_fz = sum_fz * sum_gauss / p_number;

      fx[tid1] = sum_fx;
      fy[tid1] = sum_fy;
      fz[tid1] = sum_fz;

      //virial
      for(int i =0;i<6;++i) {
        sum_virial[i] = sum_virial[i] * sum_gauss / p_number;
      }
      for(int i =0;i<6;++i) {
        flat_virial[  i* num_atoms + tid1] = sum_virial[i];
      }
    }
  }




__device__ inline rbmd::Real Gaussian_Fourier_Plus(
    rbmd::Real Kx, rbmd::Real Ky, rbmd::Real Kz,
    rbmd::Real sigma, rbmd::Real b, int Mmax,
    const rbmd::Real* __restrict__ d_coef)
{
    rbmd::Real k2 = Kx * Kx + Ky * Ky + Kz * Kz;
    rbmd::Real sum = 0.00;
    rbmd::Real sigma2 = (-sigma * sigma*0.5);
    rbmd::Real b2 = b * b;

    for (int i = 0; i < Mmax; i++)
    {
       rbmd::Real b_2i = POW(b2,i);
       rbmd::Real  mid = b_2i* sigma2;
        mid = mid * k2;
       rbmd::Real exp = EXP(mid);
        sum = sum + d_coef[i] * exp;
    }
    return sum;
}

__device__ inline rbmd::Real Gaussian_Fourier_Plus_modify(
    rbmd::Real Kx, rbmd::Real Ky, rbmd::Real Kz,
    rbmd::Real sigma, rbmd::Real b, int Mmax,
    const rbmd::Real* __restrict__ d_sl,
    const rbmd::Real* __restrict__ d_coef_npt)
{
    rbmd::Real k2 = Kx * Kx + Ky * Ky + Kz * Kz;
    rbmd::Real sum = 0.00;
    rbmd::Real sigma2 = sigma * sigma;
    rbmd::Real b2 = b * b;

    for (int i = 0; i < Mmax; i++)
    {
        sum += d_coef_npt[i] * EXP(-(d_sl[i] * d_sl[i] * 0.5) * k2);
    }
    return sum;
}


//----------------------------------------------------------------------
// 2. ComputeRBSOGFactorsOp
//----------------------------------------------------------------------

__global__ void ComputeRBSOGFactorsKernel(
  const rbmd::Id P,const rbmd::Real* d_K_x, const rbmd::Real* d_K_y, const rbmd::Real* d_K_z,
  const rbmd::Real* d_K_npt_x, const rbmd::Real* d_K_npt_y, const rbmd::Real* d_K_npt_z,
  Box box,const rbmd::Real sigma, const rbmd::Real b, const rbmd::Id Mmax,
  const rbmd::Real* d_coef, const rbmd::Real* d_coef_npt,
  const rbmd::Real L_ratio, const rbmd::Real S_ratio,
  const rbmd::Real S_npt_ratio, rbmd::Real* d_fac, rbmd::Real* d_fac_npt)
{
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= P) return;

    // --- Force Factor (fac) ---
    rbmd::Real Kx = d_K_x[tid];
    rbmd::Real Ky = d_K_y[tid];
    rbmd::Real Kz = d_K_z[tid];

    rbmd::Real Kx0 = Kx * L_ratio;
    rbmd::Real Ky0 = Ky * L_ratio;
    rbmd::Real Kz0 = Kz * L_ratio;

    rbmd::Real K2 = Kx * Kx + Ky * Ky + Kz * Kz;
    rbmd::Real K2_0 = Kx0 * Kx0 + Ky0 * Ky0 + Kz0 * Kz0;

   //Gaussian_Fourier_Plus
    rbmd::Real sum = 0.00;
    rbmd::Real sum_mid0 = 0.00;

    rbmd::Real sigma2 = (-sigma * sigma*0.5);
    rbmd::Real b2 = b * b;

    for (int i = 0; i < Mmax; i++)
    {
      rbmd::Real b_2i = POW(b2,i);
      rbmd::Real  mid = b_2i* sigma2;
      rbmd::Real  mid0 = b_2i* sigma2;

      mid = mid * K2;
      mid0 = mid0 * K2_0;

      rbmd::Real exp = EXP(mid);
      rbmd::Real exp0 = EXP(mid0);

      sum = sum + d_coef[i] * exp;
      sum_mid0 = sum_mid0 + d_coef[i] * exp0;
    }
    rbmd::Real mid = sum* K2;
    rbmd::Real mid0 = sum_mid0* K2_0;
    // rbmd::Real mid = Gaussian_Fourier_Plus(Kx, Ky, Kz, sigma, b, Mmax, d_coef) * K2;
    // rbmd::Real mid0 = Gaussian_Fourier_Plus(Kx0, Ky0, Kz0, sigma, b, Mmax, d_coef) * K2_0;

    d_fac[tid] = S_ratio * mid / mid0;

    // --- NPT  (fac_npt) ---
    rbmd::Real Kx_npt = d_K_npt_x[tid];
    rbmd::Real Ky_npt = d_K_npt_y[tid];
    rbmd::Real Kz_npt = d_K_npt_z[tid];

    rbmd::Real Kx_npt0 = Kx_npt * L_ratio;
    rbmd::Real Ky_npt0 = Ky_npt * L_ratio;
    rbmd::Real Kz_npt0 = Kz_npt * L_ratio;

    rbmd::Real K2_npt = Kx_npt * Kx_npt + Ky_npt * Ky_npt + Kz_npt * Kz_npt;
    rbmd::Real K2_npt_0 = Kx_npt0 * Kx_npt0 + Ky_npt0 * Ky_npt0 + Kz_npt0 * Kz_npt0;
    rbmd::Real K4_npt = K2_npt * K2_npt;
    rbmd::Real K4_npt_0 = K2_npt_0 * K2_npt_0;

    rbmd::Real sum_npt = 0.00;
    rbmd::Real sum_npt0 = 0.00;

  for (int i = 0; i < Mmax; i++)
  {
    rbmd::Real b_2i = POW(b2,i);
    rbmd::Real  mid_npt = b_2i* sigma2;
    rbmd::Real  mid_npt0 = b_2i* sigma2;

    mid_npt = mid_npt * K2_npt;
    mid_npt0 = mid_npt0 * K2_npt_0;

    rbmd::Real exp_npt = EXP(mid_npt);
    rbmd::Real exp_npt0 = EXP(mid_npt0);

    sum_npt = sum_npt + d_coef[i] * exp_npt;
    sum_npt0 = sum_mid0 + d_coef_npt[i] * exp_npt0;
  }
    rbmd::Real mid_npt = sum_npt* K4_npt;
    rbmd::Real mid0_npt0 = sum_npt0* K4_npt_0;
    // rbmd::Real mid_npt = Gaussian_Fourier_Plus(Kx_npt, Ky_npt, Kz_npt, sigma, b, Mmax, d_coef_npt) * K4_npt;
    // rbmd::Real mid0_npt = Gaussian_Fourier_Plus(Kx_npt0, Ky_npt0, Kz_npt0, sigma, b, Mmax, d_coef_npt) * K4_npt_0;

    d_fac_npt[tid] = S_npt_ratio * mid_npt / mid0_npt0;
}

void ComputeRBSOGFactorsOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id P,const rbmd::Real* d_K_x, const rbmd::Real* d_K_y, const rbmd::Real* d_K_z,
    const rbmd::Real* d_K_npt_x, const rbmd::Real* d_K_npt_y, const rbmd::Real* d_K_npt_z,
    Box box,const rbmd::Real sigma, const rbmd::Real b, const rbmd::Id Mmax,
    const rbmd::Real* d_coef, const rbmd::Real* d_coef_npt,
    const rbmd::Real L_ratio, const rbmd::Real S_ratio,
    const rbmd::Real S_npt_ratio,rbmd::Real* d_fac, rbmd::Real* d_fac_npt)
{
    unsigned int blocks_per_grid = (P + BLOCK_SIZE - 1) / BLOCK_SIZE;
    CHECK_KERNEL(ComputeRBSOGFactorsKernel<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        P, d_K_x, d_K_y, d_K_z, d_K_npt_x, d_K_npt_y, d_K_npt_z,
        box, sigma, b, Mmax, d_coef, d_coef_npt,
        L_ratio, S_ratio, S_npt_ratio,d_fac, d_fac_npt
    ));
}

__global__ void ComputePnumberChargeStructureFactorSOG(
    Box box, const rbmd::Id num_atoms, const rbmd::Id p_number,
    const rbmd::Real* __restrict__ charge,
    const rbmd::Real* __restrict__ p_sample_x,
    const rbmd::Real* __restrict__ p_sample_y,
    const rbmd::Real* __restrict__ p_sample_z,
    const rbmd::Real* __restrict__ px,
    const rbmd::Real* __restrict__ py,
    const rbmd::Real* __restrict__ pz,
    rbmd::Real* __restrict__ rhok_real,
    rbmd::Real* __restrict__ rhok_image)
{
  const rbmd::Id P_index = blockIdx.x;
  if (P_index >= p_number) return;

  //
  rbmd::Real k_x = __ldg(&p_sample_x[P_index]);
  rbmd::Real k_y = __ldg(&p_sample_y[P_index]);
  rbmd::Real k_z = __ldg(&p_sample_z[P_index]);

  rbmd::Real local_real_sum = 0.0;
  rbmd::Real local_image_sum = 0.0;

  // Grid-stride
  for (rbmd::Id N_index = threadIdx.x; N_index < num_atoms; N_index += blockDim.x) {
    if (N_index >= num_atoms) continue;

    rbmd::Real chargei = __ldg(&charge[N_index]);
    rbmd::Real p_x = __ldg(&px[N_index]);
    rbmd::Real p_y = __ldg(&py[N_index]);
    rbmd::Real p_z = __ldg(&pz[N_index]);

    rbmd::Real dot_product = k_x * p_x + k_y * p_y + k_z * p_z;

    local_real_sum += chargei * COS(dot_product);
    local_image_sum += chargei * SIN(dot_product);
  }

  //
  rbmd::Real total_real = blockReduceSum(local_real_sum);
  rbmd::Real total_image = blockReduceSum(local_image_sum);

  //
  if (threadIdx.x == 0) {
    rhok_real[P_index] = total_real;
    rhok_image[P_index] = total_image;
  }
}

//----------------------------------------------------------------------
// 3. ComputeDirectChargeStructureFactorOp
//----------------------------------------------------------------------

__global__ void ComputeDirectChargeStructureFactorKernel(
    const rbmd::Id num_atoms, const rbmd::Id num_k_direct,
    const Real3* __restrict__ d_k_direct,const rbmd::Real* __restrict__ charge,
    const rbmd::Real* __restrict__ px,const rbmd::Real* __restrict__ py,
    const rbmd::Real* __restrict__ pz,rbmd::Real* __restrict__ d_rho_direct_real,
    rbmd::Real* __restrict__ d_rho_direct_imag)
{
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_atoms) return;

    // Load atom data
    const rbmd::Real chargei = __ldg(&charge[tid]);
    const rbmd::Real p_x = __ldg(&px[tid]);
    const rbmd::Real p_y = __ldg(&py[tid]);
    const rbmd::Real p_z = __ldg(&pz[tid]);

    // Loop over all direct K-vectors
    for (int k_idx = 0; k_idx < num_k_direct; k_idx++) {
        const Real3 K = d_k_direct[k_idx];

        rbmd::Real dot = K.x * p_x + K.y * p_y + K.z * p_z;
        rbmd::Real real_val = chargei * COS(dot);
        rbmd::Real imag_val = chargei * SIN(dot);

        // Atomically add this atom's contribution to the total for this K-vector
        atomicAdd(&d_rho_direct_real[k_idx], real_val);
        atomicAdd(&d_rho_direct_imag[k_idx], imag_val);
    }
}

__global__ void ComputeDirectChargeStructureFactorKernel_opt(
  const rbmd::Id num_atoms, const rbmd::Id num_k_direct,
  const rbmd::Real*  d_k_direct_x,  const rbmd::Real*  d_k_direct_y,
  const rbmd::Real*  d_k_direct_z,const rbmd::Real* charge,
  const rbmd::Real* px,const rbmd::Real* py,
  const rbmd::Real* pz,rbmd::Real* d_rho_direct_real,
  rbmd::Real* d_rho_direct_imag)
{
  const rbmd::Id k_index = blockIdx.x;
  if (k_index >= num_k_direct) return;

  rbmd::Real k_x =d_k_direct_x[k_index];
  rbmd::Real k_y =d_k_direct_y[k_index];
  rbmd::Real k_z =d_k_direct_z[k_index];

  rbmd::Real local_real_sum = 0.0;
  rbmd::Real local_image_sum = 0.0;

  // Grid-stride
  for (rbmd::Id N_index = threadIdx.x; N_index < num_atoms; N_index += blockDim.x) {
    if (N_index >= num_atoms) continue;

    rbmd::Real chargei = __ldg(&charge[N_index]);
    rbmd::Real p_x = __ldg(&px[N_index]);
    rbmd::Real p_y = __ldg(&py[N_index]);
    rbmd::Real p_z = __ldg(&pz[N_index]);

    rbmd::Real dot_product = k_x * p_x + k_y * p_y + k_z * p_z;

    local_real_sum += chargei * COS(dot_product);
    local_image_sum += chargei * SIN(dot_product);
  }

  //
  rbmd::Real total_real = blockReduceSum(local_real_sum);
  rbmd::Real total_image = blockReduceSum(local_image_sum);

  //
  if (threadIdx.x == 0) {
    d_rho_direct_real[k_index] = total_real;
    d_rho_direct_imag[k_index] = total_image;
  }
}

void ComputeDirectChargeStructureFactorOp<device::DEVICE_GPU>::operator()(
  const rbmd::Id num_atoms, const rbmd::Id num_k_direct,
  const rbmd::Real*  d_k_direct_x,  const rbmd::Real*  d_k_direct_y,
  const rbmd::Real*  d_k_direct_z,const rbmd::Real* charge,
  const rbmd::Real* px,const rbmd::Real* py,
  const rbmd::Real* pz,rbmd::Real* d_rho_direct_real,
  rbmd::Real* d_rho_direct_imag)
{
    // unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    // CHECK_KERNEL(ComputeDirectChargeStructureFactorKernel<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
    //     num_atoms, num_k_direct, d_k_direct, charge, px, py, pz,d_rho_direct_real, d_rho_direct_imag
    // ));

  unsigned int blocks_per_grid = num_k_direct;
  CHECK_KERNEL(ComputeDirectChargeStructureFactorKernel_opt<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, num_k_direct, d_k_direct_x,d_k_direct_y, d_k_direct_z,charge, px, py, pz,d_rho_direct_real, d_rho_direct_imag
  ));
}


//----------------------------------------------------------------------
// 4. ComputeRBSOGSampleForceOp (Force/Virial Kernel + Energy Kernel)
//----------------------------------------------------------------------

__global__ void ComputeRBSOGSampleForceKernel(
    Box box, const rbmd::Id num_atoms, const rbmd::Id P,
    rbmd::Real qqr2e, rbmd::Real S0_sample, rbmd::Real S_npt_sample,
    const rbmd::Real* __restrict__ d_K_x,const rbmd::Real* __restrict__ d_K_y,
    const rbmd::Real* __restrict__ d_K_z,const rbmd::Real* __restrict__ d_K_npt_x,
    const rbmd::Real* __restrict__ d_K_npt_y,const rbmd::Real* __restrict__ d_K_npt_z,
    const rbmd::Id* __restrict__ d_idx_npt_all,const rbmd::Real* __restrict__ d_fac,
    const rbmd::Real* __restrict__ d_fac_npt,const rbmd::Real* __restrict__ d_rho_real,
    const rbmd::Real* __restrict__ d_rho_imag,const rbmd::Real* __restrict__ charge,
    const rbmd::Real* __restrict__ px,const rbmd::Real* __restrict__ py,
    const rbmd::Real* __restrict__ pz,rbmd::Real* __restrict__ fx,
    rbmd::Real* __restrict__ fy,rbmd::Real* __restrict__ fz)
{
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_atoms) return;

    const rbmd::Real V = box._length[0] * box._length[1] * box._length[2];

    // Pre-calculate constants from rbsog_intel.cpp
    const rbmd::Real MIDTERM_Force = - (S0_sample / P) * qqr2e / V;
    const rbmd::Real Moment_Term_Virial = (S0_sample / P) * qqr2e / (2.0 * V);
    const rbmd::Real Moment_Term_NPT_Virial = - (S_npt_sample / P) * qqr2e / (4.0 * V);

    // Load atom data
    const rbmd::Real chargei = __ldg(&charge[tid]);
    const rbmd::Real p_x = __ldg(&px[tid]);
    const rbmd::Real p_y = __ldg(&py[tid]);
    const rbmd::Real p_z = __ldg(&pz[tid]);

    rbmd::Real sum_fx = 0.0;
    rbmd::Real sum_fy = 0.0;
    rbmd::Real sum_fz = 0.0;
    rbmd::Real sum_virial[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Loop over all P sampled K-vectors
    for (int j = 0; j < P; j++) {
        // --- Load Sample Data ---
        const rbmd::Real Kx = __ldg(&d_K_x[j]);
        const rbmd::Real Ky = __ldg(&d_K_y[j]);
        const rbmd::Real Kz = __ldg(&d_K_z[j]);
        const rbmd::Real K2 = Kx * Kx + Ky * Ky + Kz * Kz;

        const rbmd::Real Kx_npt = __ldg(&d_K_npt_x[j]);
        const rbmd::Real Ky_npt = __ldg(&d_K_npt_y[j]);
        const rbmd::Real Kz_npt = __ldg(&d_K_npt_z[j]);
        const rbmd::Real K2_npt = Kx_npt * Kx_npt + Ky_npt * Ky_npt + Kz_npt * Kz_npt;
        const rbmd::Real K4_npt = K2_npt * K2_npt;

        const rbmd::Real fac_j = __ldg(&d_fac[j]);
        const rbmd::Real fac_npt_j = __ldg(&d_fac_npt[j]);

        const rbmd::Real rho_real_j = __ldg(&d_rho_real[j]);
        const rbmd::Real rho_imag_j = __ldg(&d_rho_imag[j]);
        const rbmd::Real rho_sq_j = rho_real_j * rho_real_j + rho_imag_j * rho_imag_j;

        // --- Force Calculation ---
        const rbmd::Real dot = Kx * p_x + Ky * p_y + Kz * p_z;
        rbmd::Real  s_dot = SIN(dot);
        rbmd::Real c_dot = COS(dot);

        const rbmd::Real force_mid_term = (c_dot * rho_imag_j - s_dot * rho_real_j) * (MIDTERM_Force * fac_j) / K2;

        sum_fx += chargei * force_mid_term * Kx;
        sum_fy += chargei * force_mid_term * Ky;
        sum_fz += chargei * force_mid_term * Kz;

        // --- Virial Calculation ---
     // const rbmd::Id indx = __ldg(&d_idx_npt_all[j]); // Get the gathered index
     // const rbmd::Real rho_real_indx = __ldg(&d_rho_real[indx]);
     // const rbmd::Real rho_imag_indx = __ldg(&d_rho_imag[indx]);
     // rbmd::Real coef1 = fac_j * (rho_real_j* rho_real_j + rho_imag_j * rho_imag_j)/ K2;
     // rbmd::Real coef2 = fac_npt_j * (rho_real_indx* rho_real_indx + rho_imag_indx* rho_imag_indx)/ K4_npt;

     //um_virial[0] += coef1 * Moment_Term_Virial + coef2 * Moment_Term_Virial *   Kx_npt * Kx_npt;
     //um_virial[1] += coef1 * Moment_Term_Virial + coef2 * Moment_Term_NPT_Virial * Ky_npt * Ky_npt;
     //um_virial[2] += coef1 * Moment_Term_Virial + coef2 * Moment_Term_NPT_Virial * Kz_npt * Kz_npt;
     //um_virial[3] += coef2 * Moment_Term_NPT_Virial * Kx_npt * Ky_npt;
     //um_virial[4] += coef2 * Moment_Term_NPT_Virial * Kx_npt * Kz_npt;
     //um_virial[5] += coef2 * Moment_Term_NPT_Virial * Ky_npt * Kz_npt;
    }

    // Write Force
    fx[tid] = sum_fx;
    fy[tid] = sum_fy;
    fz[tid] = sum_fz;

}

__global__ void ComputeRBSOGSampleEnergyKernel(
    const rbmd::Id P, Box box, rbmd::Real qqr2e, rbmd::Real S0_sample,
    const rbmd::Real* __restrict__ d_K_x,
    const rbmd::Real* __restrict__ d_K_y,
    const rbmd::Real* __restrict__ d_K_z,
    const rbmd::Real* __restrict__ d_fac,
    const rbmd::Real* __restrict__ d_rho_real,
    const rbmd::Real* __restrict__ d_rho_imag,
    rbmd::Real* d_energy_parts) // Output (index 0)
{
    rbmd::Real sum_energy = 0.0;
    rbmd::Real V = box._length[0]*box._length[1]*box._length[2];

    const rbmd::Real P_real = (rbmd::Real)P;
    const rbmd::Real Moment_Term_Energy = (S0_sample / (2.0 * P_real * V)) * qqr2e;

    // Parallel reduction over P
    for (int j = threadIdx.x; j < P; j += blockDim.x) {
        const rbmd::Real Kx = __ldg(&d_K_x[j]);
        const rbmd::Real Ky = __ldg(&d_K_y[j]);
        const rbmd::Real Kz = __ldg(&d_K_z[j]);
        const rbmd::Real K2 = Kx * Kx + Ky * Ky + Kz * Kz;

        const rbmd::Real fac_j = __ldg(&d_fac[j]);
        const rbmd::Real rho_real_j = __ldg(&d_rho_real[j]);
        const rbmd::Real rho_imag_j = __ldg(&d_rho_imag[j]);
        const rbmd::Real rho_sq_j = rho_real_j * rho_real_j + rho_imag_j * rho_imag_j;

        sum_energy += fac_j * Moment_Term_Energy * rho_sq_j / K2;
    }

    block_reduce_sum(sum_energy);

    if (threadIdx.x == 0) {
        d_energy_parts[0] = sum_energy;
    }
}

__global__ void ComputeRBSOGSampleEnergyVirialKernel(
    const rbmd::Id P, Box box, const rbmd::Real qqr2e,
    const rbmd::Real S0_sample, const rbmd::Real S_npt_sample,
    const rbmd::Real* __restrict__ d_K_x, const rbmd::Real* __restrict__ d_K_y, const rbmd::Real* __restrict__ d_K_z,
    const rbmd::Real* __restrict__ d_K_npt_x, const rbmd::Real* __restrict__ d_K_npt_y, const rbmd::Real* __restrict__ d_K_npt_z,
    const rbmd::Id* __restrict__ d_idx_npt_all,
    const rbmd::Real* __restrict__ d_fac, const rbmd::Real* __restrict__ d_fac_npt,
    const rbmd::Real* __restrict__ d_rho_real, const rbmd::Real* __restrict__ d_rho_imag,
    rbmd::Real* __restrict__ global_virial, rbmd::Real* __restrict__ d_energy_parts)
{
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_e;
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_v[6];

    rbmd::Real sum_e = 0.0;
    rbmd::Real sum_v[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    const rbmd::Real V = box._length[0] * box._length[1] * box._length[2];
    const rbmd::Real P_real = (rbmd::Real)P;

    // (rbsog_intel.cpp: compute(), eflag_global)
    const rbmd::Real Moment_Term_Energy = (S0_sample / P_real) * qqr2e / (2.0 * V);
    // (rbsog_intel.cpp: compute(), vflag_global)
    const rbmd::Real Moment_Term_Virial     = (S0_sample / P_real) * qqr2e / (2.0 * V);
    const rbmd::Real Moment_Term_NPT_Virial = - (S_npt_sample / P_real) * qqr2e / (4.0 * V);

    // 使用单块(single-block)网格步长循环 (grid-stride loop)
    for (int j = threadIdx.x; j < P; j += blockDim.x) {
        const rbmd::Real Kx = __ldg(&d_K_x[j]);
        const rbmd::Real Ky = __ldg(&d_K_y[j]);
        const rbmd::Real Kz = __ldg(&d_K_z[j]);
        const rbmd::Real K2 = Kx * Kx + Ky * Ky + Kz * Kz;

        const rbmd::Real Kx_npt = __ldg(&d_K_npt_x[j]);
        const rbmd::Real Ky_npt = __ldg(&d_K_npt_y[j]);
        const rbmd::Real Kz_npt = __ldg(&d_K_npt_z[j]);
        const rbmd::Real K2_npt = Kx_npt * Kx_npt + Ky_npt * Ky_npt + Kz_npt * Kz_npt;
        const rbmd::Real K4_npt = K2_npt * K2_npt;

        const rbmd::Real fac_j = __ldg(&d_fac[j]);
        const rbmd::Real fac_npt_j = __ldg(&d_fac_npt[j]);

        const rbmd::Real rho_real_j = __ldg(&d_rho_real[j]);
        const rbmd::Real rho_imag_j = __ldg(&d_rho_imag[j]);
        const rbmd::Real rho_sq_j = rho_real_j * rho_real_j + rho_imag_j * rho_imag_j;

        // --- Virial & Energy Gather ---
        const rbmd::Id indx = __ldg(&d_idx_npt_all[j]);
        const rbmd::Real rho_real_indx = __ldg(&d_rho_real[indx]);
        const rbmd::Real rho_imag_indx = __ldg(&d_rho_imag[indx]);
        const rbmd::Real rho_sq_indx = rho_real_indx * rho_real_indx + rho_imag_indx * rho_imag_indx;

        // --- Energy Calculation ---
        const rbmd::Real coef1_e = fac_j * rho_sq_j / K2;
        sum_e += coef1_e * Moment_Term_Energy;

        // --- Virial Calculation ---
        const rbmd::Real coef1_v = fac_j * rho_sq_j / K2;
        const rbmd::Real coef2_v = fac_npt_j * rho_sq_indx / K4_npt;

        // (!! 修复了您代码中的数学错误 !!)
        sum_v[0] += coef1_v * Moment_Term_Virial + coef2_v * Moment_Term_NPT_Virial * Kx_npt * Kx_npt;
        sum_v[1] += coef1_v * Moment_Term_Virial + coef2_v * Moment_Term_NPT_Virial * Ky_npt * Ky_npt;
        sum_v[2] += coef1_v * Moment_Term_Virial + coef2_v * Moment_Term_NPT_Virial * Kz_npt * Kz_npt;
        sum_v[3] += coef2_v * Moment_Term_NPT_Virial * Kx_npt * Ky_npt;
        sum_v[4] += coef2_v * Moment_Term_NPT_Virial * Kx_npt * Kz_npt;
        sum_v[5] += coef2_v * Moment_Term_NPT_Virial * Ky_npt * Kz_npt;
    }

    // --- 最终规约 (Final Reduction) ---
    rbmd::Real block_sum_e = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage_e).Sum(sum_e);
    if (threadIdx.x == 0) {
        atomicAdd(&d_energy_parts[0], block_sum_e);
    }

    for (int i = 0; i < 6; ++i) {
        rbmd::Real block_sum_v = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage_v[i]).Sum(sum_v[i]);
        if (threadIdx.x == 0) {
            atomicAdd(&global_virial[i], block_sum_v);
        }
    }
}

__global__ void ComputeRBSOGSampleVirialKernel(
    const rbmd::Id P, Box box, rbmd::Real qqr2e, rbmd::Real S0_sample,rbmd::Real S_npt_sample,
    const rbmd::Real* __restrict__ d_K_x,const rbmd::Real* __restrict__ d_K_y,
    const rbmd::Real* __restrict__ d_K_z,const rbmd::Real* __restrict__ d_K_npt_x,
    const rbmd::Real* __restrict__ d_K_npt_y,const rbmd::Real* __restrict__ d_K_npt_z,
    const rbmd::Id* __restrict__ d_idx_npt_all,const rbmd::Real* __restrict__ d_fac,
    const rbmd::Real* __restrict__ d_fac_npt,const rbmd::Real* __restrict__ d_rho_real,
    const rbmd::Real* __restrict__ d_rho_imag,rbmd::Real* global_virial)
{
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_v[6];

    rbmd::Real V = box._length[0]*box._length[1]*box._length[2];
    const rbmd::Real Moment_Term_Virial = (S0_sample / P) * qqr2e / (2.0 * V);
    const rbmd::Real Moment_Term_NPT_Virial = - (S_npt_sample / P) * qqr2e / (4.0 * V);

    rbmd::Real sum_virial[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // Parallel reduction over P
    for (int j = threadIdx.x; j < P; j += blockDim.x) {
      const rbmd::Real Kx = __ldg(&d_K_x[j]);
      const rbmd::Real Ky = __ldg(&d_K_y[j]);
      const rbmd::Real Kz = __ldg(&d_K_z[j]);
      const rbmd::Real K2 = Kx * Kx + Ky * Ky + Kz * Kz;

      const rbmd::Real Kx_npt = __ldg(&d_K_npt_x[j]);
      const rbmd::Real Ky_npt = __ldg(&d_K_npt_y[j]);
      const rbmd::Real Kz_npt = __ldg(&d_K_npt_z[j]);
      const rbmd::Real K2_npt = Kx_npt * Kx_npt + Ky_npt * Ky_npt + Kz_npt * Kz_npt;
      const rbmd::Real K4_npt = K2_npt * K2_npt;

      const rbmd::Real fac_j = __ldg(&d_fac[j]);
      const rbmd::Real fac_npt_j = __ldg(&d_fac_npt[j]);

      const rbmd::Real rho_real_j = __ldg(&d_rho_real[j]);
      const rbmd::Real rho_imag_j = __ldg(&d_rho_imag[j]);

      // --- Virial Calculation ---
      const rbmd::Id indx = __ldg(&d_idx_npt_all[j]); // Get the gathered index
      const rbmd::Real rho_real_indx = __ldg(&d_rho_real[indx]);
      const rbmd::Real rho_imag_indx = __ldg(&d_rho_imag[indx]);
      rbmd::Real coef1 = fac_j * (rho_real_j* rho_real_j + rho_imag_j * rho_imag_j)/ K2;
      rbmd::Real coef2 = fac_npt_j * (rho_real_indx* rho_real_indx + rho_imag_indx* rho_imag_indx)/ K4_npt;

      sum_virial[0] += coef1 * Moment_Term_Virial + coef2 * Moment_Term_Virial *   Kx_npt * Kx_npt;
      sum_virial[1] += coef1 * Moment_Term_Virial + coef2 * Moment_Term_NPT_Virial * Ky_npt * Ky_npt;
      sum_virial[2] += coef1 * Moment_Term_Virial + coef2 * Moment_Term_NPT_Virial * Kz_npt * Kz_npt;
      sum_virial[3] += coef2 * Moment_Term_NPT_Virial * Kx_npt * Ky_npt;
      sum_virial[4] += coef2 * Moment_Term_NPT_Virial * Kx_npt * Kz_npt;
      sum_virial[5] += coef2 * Moment_Term_NPT_Virial * Ky_npt * Kz_npt;
    }

    //
    for (int i = 0; i < 6; ++i) {
      rbmd::Real block_sum_v = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage_v[i]).Sum(sum_virial[i]);
      if (threadIdx.x == 0) {
        atomicAdd(&global_virial[i], block_sum_v);
      }
    }
}



void ComputeRBSOGSampleForceOp<device::DEVICE_GPU>::operator()(
    Box box, const rbmd::Id num_atoms, const rbmd::Id P,
    const rbmd::Real qqr2e, const rbmd::Real S0_sample, const rbmd::Real S_npt_sample,
    const rbmd::Real* d_K_x, const rbmd::Real* d_K_y, const rbmd::Real* d_K_z,
    const rbmd::Real* d_K_npt_x, const rbmd::Real* d_K_npt_y, const rbmd::Real* d_K_npt_z,
    const rbmd::Id* d_idx_npt_all,const rbmd::Real* d_fac, const rbmd::Real* d_fac_npt,
    const rbmd::Real* d_rho_real, const rbmd::Real* d_rho_imag,
    const rbmd::Real* charge,const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,rbmd::Real* global_virial,
    rbmd::Real* d_energy_parts)
{
    // --- Kernel 1: Force and Virial ---
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    CHECK_KERNEL(ComputeRBSOGSampleForceKernel<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, num_atoms, P, qqr2e, S0_sample, S_npt_sample,
        d_K_x, d_K_y, d_K_z, d_K_npt_x, d_K_npt_y, d_K_npt_z,
        d_idx_npt_all, d_fac, d_fac_npt,
        d_rho_real, d_rho_imag,
        charge, px, py, pz,
        fx, fy, fz));

    // --- Kernel 2: Energy ---
    CHECK_KERNEL(ComputeRBSOGSampleEnergyKernel<<<1, BLOCK_SIZE, 0, 0>>>(
        P, box, qqr2e, S0_sample,
        d_K_x, d_K_y, d_K_z,
        d_fac, d_rho_real, d_rho_imag,
        d_energy_parts // Output to index 0
    ));

    //--- Kernel 3:
    CHECK_KERNEL(ComputeRBSOGSampleVirialKernel<<<1, BLOCK_SIZE, 0, 0>>>(
    P, box, qqr2e, S0_sample,S_npt_sample,
    d_K_x, d_K_y, d_K_z,d_K_npt_x, d_K_npt_y, d_K_npt_z,
    d_idx_npt_all, d_fac, d_fac_npt,
     d_rho_real, d_rho_imag,
    global_virial ));
}


//----------------------------------------------------------------------
// 5. ComputeRBSOGDirectForceOp (Force/Virial Kernel + Energy Kernel)
//----------------------------------------------------------------------

__global__ void ComputeRBSOGDirectForceKernel(
  Box box, const rbmd::Id num_atoms, const rbmd::Id num_k_direct,
  const rbmd::Real qqr2e,    const rbmd::Real*  d_k_direct_x,
  const rbmd::Real*  d_k_direct_y,const rbmd::Real*  d_k_direct_z,
  const rbmd::Real* d_f_b_sigma,const rbmd::Real* d_f_b_sigma_npt,
  const rbmd::Real* d_rho_direct_real,const rbmd::Real* d_rho_direct_imag,
  const rbmd::Real* charge,const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
  rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz)
{
    rbmd::Real sum_fx = 0.0;
    rbmd::Real sum_fy = 0.0;
    rbmd::Real sum_fz = 0.0;
    rbmd::Real sum_virial[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_atoms) return;

    const rbmd::Real V = box._length[0] * box._length[1] * box._length[2];

    // Pre-calculate constants
    const rbmd::Real MID_Force = -qqr2e / V;
    // Again, virial is global in rbsog_intel. Skipping per-atom virial.
    const rbmd::Real c1_Virial = qqr2e / (2.0 * V);
    const rbmd::Real c2_Virial = -qqr2e / (4.0 * V);

    // Load atom data
    const rbmd::Real chargei = __ldg(&charge[tid]);
    const rbmd::Real p_x = __ldg(&px[tid]);
    const rbmd::Real p_y = __ldg(&py[tid]);
    const rbmd::Real p_z = __ldg(&pz[tid]);

    // Loop over all direct K-vectors
    for (int k_idx = 0; k_idx < num_k_direct; k_idx++) {
        const rbmd::Real  Kx = d_k_direct_x[k_idx];
        const rbmd::Real  Ky = d_k_direct_y[k_idx];
        const rbmd::Real  Kz = d_k_direct_z[k_idx];
        const Real3 K = {Kx , Ky, Kz};

        const rbmd::Real f_b_sigma = __ldg(&d_f_b_sigma[k_idx]);
        const rbmd::Real f_b_sigma_npt = __ldg(&d_f_b_sigma_npt[k_idx]); // For virial

        const rbmd::Real rho_real = __ldg(&d_rho_direct_real[k_idx]);
        const rbmd::Real rho_imag = __ldg(&d_rho_direct_imag[k_idx]);
        const rbmd::Real rho_sq = rho_real * rho_real + rho_imag * rho_imag; // For virial/energy

        // --- Force ---
        const rbmd::Real dot = K.x * p_x + K.y * p_y + K.z * p_z;

        rbmd::Real  s_dot = SIN(dot);
        rbmd::Real c_dot = COS(dot);
        const rbmd::Real force_mid_term = (c_dot * rho_imag - s_dot * rho_real) * (f_b_sigma * MID_Force);

        sum_fx += chargei * force_mid_term * K.x;
        sum_fy += chargei * force_mid_term * K.y;
        sum_fz += chargei * force_mid_term * K.z;

       //// --- Virial---
       //sum_virial[0] += c1_Virial * f_b_sigma * rho_sq + c2_Virial * f_b_sigma_npt * rho_sq * Kx * Kx;
       //sum_virial[1] += c1_Virial * f_b_sigma * rho_sq + c2_Virial * f_b_sigma_npt * rho_sq * Ky * Ky;
       //sum_virial[2] += c1_Virial * f_b_sigma * rho_sq + c2_Virial * f_b_sigma_npt * rho_sq * Kz * Kz;

       //sum_virial[3] += c2_Virial * f_b_sigma_npt *   rho_sq  * Kx * Ky;
       //sum_virial[4] += c1_Virial * f_b_sigma_npt *   rho_sq  * Kx * Kz;
       //sum_virial[5] += c1_Virial * f_b_sigma_npt *   rho_sq  * Ky * Kz;

    }

    // Atomically ADD to existing forces
    fx[tid] = sum_fx;
    fy[tid] = sum_fy;
    fz[tid] = sum_fz;

}

__global__ void ComputeRBSOGDirectEnergyKernel(
    const rbmd::Id num_k_direct,   Box box, rbmd::Real qqr2e,
    const rbmd::Real* __restrict__ d_f_b_sigma,
    const rbmd::Real* __restrict__ d_rho_direct_real,
    const rbmd::Real* __restrict__ d_rho_direct_imag,
    rbmd::Real* d_energy_parts) // Output (index 1)
{
    rbmd::Real sum_energy = 0.0;
    rbmd::Real V = box._length[0] * box._length[1] * box._length[2];
    const rbmd::Real c1_Energy = qqr2e / (2.0 * V);

    // Parallel reduction over num_k_direct
    for (int k_idx = threadIdx.x; k_idx < num_k_direct; k_idx += blockDim.x) {
        const rbmd::Real f_b = __ldg(&d_f_b_sigma[k_idx]);
        const rbmd::Real rho_real = __ldg(&d_rho_direct_real[k_idx]);
        const rbmd::Real rho_imag = __ldg(&d_rho_direct_imag[k_idx]);
        const rbmd::Real rho_sq = rho_real * rho_real + rho_imag * rho_imag;

        sum_energy += c1_Energy * f_b * rho_sq;
    }

    block_reduce_sum(sum_energy);

    if (threadIdx.x == 0) {
        // Atomically ADD to the energy buffer (in case sample part is also writing)
        atomicAdd(&d_energy_parts[1], sum_energy);
    }
}

__global__ void ComputeRBSOGDirectVirialKernel(
    const rbmd::Id num_k_direct,   Box box, rbmd::Real qqr2e,
    const rbmd::Real*  d_k_direct_x,const rbmd::Real*  d_k_direct_y,
    const rbmd::Real*  d_k_direct_z,
    const rbmd::Real* __restrict__ d_f_b_sigma,const rbmd::Real* d_f_b_sigma_npt,
    const rbmd::Real* __restrict__ d_rho_direct_real,
    const rbmd::Real* __restrict__ d_rho_direct_imag,
    rbmd::Real* global_virial)
{
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_v[6];

    rbmd::Real sum_virial[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    rbmd::Real V = box._length[0] * box._length[1] * box._length[2];
    const rbmd::Real c1_Virial = qqr2e / (2.0 * V);
    const rbmd::Real c2_Virial = -qqr2e / (4.0 * V);


    // Parallel reduction over num_k_direct
    for (int k_idx = threadIdx.x; k_idx < num_k_direct; k_idx += blockDim.x) {
      const rbmd::Real  Kx = d_k_direct_x[k_idx];
      const rbmd::Real  Ky = d_k_direct_y[k_idx];
      const rbmd::Real  Kz = d_k_direct_z[k_idx];

      const rbmd::Real rho_real = __ldg(&d_rho_direct_real[k_idx]);
      const rbmd::Real rho_imag = __ldg(&d_rho_direct_imag[k_idx]);
      const rbmd::Real rho_sq = rho_real * rho_real + rho_imag * rho_imag;

      const rbmd::Real f_b_sigma = __ldg(&d_f_b_sigma[k_idx]);
      const rbmd::Real f_b_sigma_npt = __ldg(&d_f_b_sigma_npt[k_idx]); // For virial


      // --- Virial---
      sum_virial[0] += c1_Virial * f_b_sigma * rho_sq + c2_Virial * f_b_sigma_npt * rho_sq * Kx * Kx;
      sum_virial[1] += c1_Virial * f_b_sigma * rho_sq + c2_Virial * f_b_sigma_npt * rho_sq * Ky * Ky;
      sum_virial[2] += c1_Virial * f_b_sigma * rho_sq + c2_Virial * f_b_sigma_npt * rho_sq * Kz * Kz;

      sum_virial[3] += c2_Virial * f_b_sigma_npt *   rho_sq  * Kx * Ky;
      sum_virial[4] += c1_Virial * f_b_sigma_npt *   rho_sq  * Kx * Kz;
      sum_virial[5] += c1_Virial * f_b_sigma_npt *   rho_sq  * Ky * Kz;

    }

    for (int i = 0; i < 6; ++i) {
      rbmd::Real block_sum_v = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage_v[i]).Sum(sum_virial[i]);
      if (threadIdx.x == 0) {
        atomicAdd(&global_virial[i], block_sum_v);
      }
    }

}

void ComputeRBSOGDirectForceOp<device::DEVICE_GPU>::operator()(
  Box box, const rbmd::Id num_atoms, const rbmd::Id num_k_direct,
  const rbmd::Real qqr2e,    const rbmd::Real*  d_k_direct_x,
  const rbmd::Real*  d_k_direct_y,const rbmd::Real*  d_k_direct_z,
  const rbmd::Real* d_f_b_sigma,const rbmd::Real* d_f_b_sigma_npt,
  const rbmd::Real* d_rho_direct_real,const rbmd::Real* d_rho_direct_imag,
  const rbmd::Real* charge,const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
  rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,rbmd::Real* global_virial,
  rbmd::Real* d_energy_parts)
{
    // --- Kernel 1: Force and Virial ---
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    CHECK_KERNEL(ComputeRBSOGDirectForceKernel<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, num_atoms, num_k_direct, qqr2e,
        d_k_direct_x,d_k_direct_y, d_k_direct_z,d_f_b_sigma, d_f_b_sigma_npt,
        d_rho_direct_real, d_rho_direct_imag,
        charge, px, py, pz,fx, fy, fz));

    // // --- Kernel 2: Energy ---
    CHECK_KERNEL(ComputeRBSOGDirectEnergyKernel<<<1, BLOCK_SIZE, 0, 0>>>(
        num_k_direct, box, qqr2e,
        d_f_b_sigma, d_rho_direct_real, d_rho_direct_imag,
        d_energy_parts // Output to index 1
    ));

    // // --- Kernel 3   ---
    CHECK_KERNEL(ComputeRBSOGDirectVirialKernel<<<1, BLOCK_SIZE, 0, 0>>>(
    num_k_direct, box, qqr2e,d_k_direct_x,d_k_direct_y, d_k_direct_z,
    d_f_b_sigma,d_f_b_sigma_npt, d_rho_direct_real, d_rho_direct_imag,
    global_virial ));
}


//////////////////////////////


// Charge Structure Factor
void ComputeChargeStructureFactorOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const Real3 K, const rbmd::Real* charge,
    const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* density_real, rbmd::Real* density_imag) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(
        ComputeChargeStructureFactor<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
            num_atoms, K, charge, px, py, pz, density_real, density_imag));
  }

// EwaldForce
void ComputeEwaldForceOp<device::DEVICE_GPU>::operator()(
     Box box, const rbmd::Id num_atoms, const Int3 Kmax,
    const rbmd::Real alpha, const rbmd::Real qqr2e,
    const rbmd::Real* real_array, const rbmd::Real* imag_array,
    const rbmd::Real* charge, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeEwaldForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, num_atoms, Kmax, alpha, qqr2e, real_array, imag_array, charge, px,
        py, pz, fx, fy, fz,flat_virial));
  }


// sq_charge
void SqchargeOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_atoms,
                                                const rbmd::Real* charge,
                                                rbmd::Real* sq_charge) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeSqCharge<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        num_atoms, charge, sq_charge));
  }

//Generate Index Array
void GenerateIndexArrayOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Id RBE_P, rbmd::Id* psample_key) {
    unsigned int blocks_per_grid =
        ((num_atoms * RBE_P) + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(GenerateIndexArray<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        num_atoms, RBE_P, psample_key));
  }

// RBE: Charge Structure Factor
void ComputePnumberChargeStructureFactorOp<device::DEVICE_GPU>::operator()(
     Box box, const rbmd::Id num_atoms, const rbmd::Id p_number,
    const rbmd::Real* charge, const rbmd::Real* p_sample_x,
    const rbmd::Real* p_sample_y, const rbmd::Real* p_sample_z,
    const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* density_real, rbmd::Real* density_imag) {
    // unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    //
    // CHECK_KERNEL(ComputePnumberChargeStructureFactor<<<blocks_per_grid,
    //                                                    BLOCK_SIZE, 0, 0>>>(
    //     box, num_atoms, p_number, charge, p_sample_x, p_sample_y, p_sample_z, px,
    //     py, pz, density_real, density_imag));

  //Opt::
  const unsigned int blocks_per_grid = p_number;
  CHECK_KERNEL(ComputePnumberChargeStructureFactorOpt<<<blocks_per_grid,
                                                     BLOCK_SIZE, 0, 0>>>(
      box, num_atoms, p_number, charge, p_sample_x, p_sample_y, p_sample_z, px,
      py, pz, density_real, density_imag));
   }

// RBE: RBE Force
void ComputeRBEForceOp<device::DEVICE_GPU>::operator()(
     Box box, const rbmd::Id num_atoms, const rbmd::Id p_number,
    const rbmd::Real alpha, const rbmd::Real qqr2e,
    const rbmd::Real* real_array, const rbmd::Real* imag_array,
    const rbmd::Real* charge, const rbmd::Real* p_sample_x,
    const rbmd::Real* p_sample_y, const rbmd::Real* p_sample_z,
    const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,rbmd::Real* flat_virial) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeRBEForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box, num_atoms, p_number, alpha, qqr2e, real_array, imag_array, charge,
        p_sample_x, p_sample_y, p_sample_z, px, py, pz, fx, fy, fz,flat_virial));
  }

  void ComputeRBEForceVirialOp<device::DEVICE_GPU>::operator()(
      Box box, const rbmd::Id p_number, const rbmd::Real alpha,
      const rbmd::Real qqrd2e, const rbmd::Real qsqsum, const rbmd::Real qsum,
      const rbmd::Real* real_array, const rbmd::Real* imag_array,
      const rbmd::Real* p_sample_x, const rbmd::Real* p_sample_y,
      const rbmd::Real* p_sample_z, rbmd::Real* virial_tensor,
      rbmd::Real* energy) {
    unsigned int blocks_per_grid = (p_number + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeRBEVirial0<<<1, BLOCK_SIZE, 0, 0>>>(
        box, p_number, alpha, qqrd2e, qsqsum, qsum, real_array, imag_array,
        p_sample_x, p_sample_y, p_sample_z, virial_tensor, energy));
  }

void ComputePnumberChargeStructureFactorSOGOp<device::DEVICE_GPU>::operator()(
     Box box, const rbmd::Id num_atoms, const rbmd::Id p_number,
    const rbmd::Real* charge, const rbmd::Real* p_sample_x,
    const rbmd::Real* p_sample_y, const rbmd::Real* p_sample_z,
    const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* density_real, rbmd::Real* density_imag) {

  //sog::
  const unsigned int blocks_per_grid = p_number;
  CHECK_KERNEL(ComputePnumberChargeStructureFactorSOG<<<blocks_per_grid,
                                                     BLOCK_SIZE, 0, 0>>>(
      box, num_atoms, p_number, charge, p_sample_x, p_sample_y, p_sample_z, px,
      py, pz, density_real, density_imag));
}

  }  // namespace op