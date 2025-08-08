#include "../common/rbmd_define.h"
#include "force/src/kspace_calculator_op/kspace_calculator_op.h"

#include "model/box.h"

static const double MY_PIS=1.77245385090551602729;
static const double MY_PI2=1.57079632679489661923;

namespace op {

//  Warp-level
__device__ __forceinline__ rbmd::Real warpReduceSum(rbmd::Real val) {
  // #pragma unroll
  for (int offset = WARP_SIZE / 2; offset > 0; offset /= 2) {
    val += __shfl_down_sync(0xFFFFFFFF, val, offset);
  }
  return val;
}

//  Block-level
__device__ __forceinline__ rbmd::Real blockReduceSum(rbmd::Real val) {

  __shared__ rbmd::Real shared[WARP_SIZE];

  int lane = threadIdx.x % WARP_SIZE;
  int wid = threadIdx.x / WARP_SIZE;

  // 1.
  val = warpReduceSum(val);

  // 2.
  if (lane == 0) {
    shared[wid] = val;
  }
  __syncthreads(); //

  // 3.
  val = (threadIdx.x < blockDim.x / WARP_SIZE) ? shared[lane] : 0.0;
  if (wid == 0) {
    val = warpReduceSum(val);
  }

  return val;
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

__global__ void ComputePnumberChargeStructureFactor1(
    Box   box, const rbmd::Id num_atoms, const rbmd::Id p_number,
    const rbmd::Real* __restrict__ charge, const rbmd::Real* __restrict__ p_sample_x,
    const rbmd::Real* __restrict__ p_sample_y, const rbmd::Real* __restrict__ p_sample_z,
    const rbmd::Real* __restrict__ px, const rbmd::Real* __restrict__ py, const rbmd::Real* __restrict__ pz,
    rbmd::Real* __restrict__ rhok_real_final, rbmd::Real* __restrict__ rhok_image_final) {

    // 1.
    const rbmd::Id P_index = blockIdx.x;
    if (P_index >= p_number) return;

    // 2.
    rbmd::Real k_x = p_sample_x[P_index] * 2 * M_PI / box._length[0];
    rbmd::Real k_y = p_sample_y[P_index] * 2 * M_PI / box._length[1];
    rbmd::Real k_z = p_sample_z[P_index] * 2 * M_PI / box._length[2];

    // 3.
    rbmd::Real local_real_sum = 0.0;
    rbmd::Real local_image_sum = 0.0;

    // 4. (Grid-Stride Loop)
    for (rbmd::Id N_index = threadIdx.x; N_index < num_atoms; N_index += blockDim.x) {
      //
      rbmd::Real chargei = charge[N_index];
      rbmd::Real p_x = px[N_index];
      rbmd::Real p_y = py[N_index];
      rbmd::Real p_z = pz[N_index];

      //
      rbmd::Real dot_product = k_x * p_x + k_y * p_y + k_z * p_z;

      //
      local_real_sum  += chargei * COS(dot_product);
      local_image_sum += chargei * SIN(dot_product);
    }

    // 5.
    rbmd::Real total_real  = blockReduceSum(local_real_sum);
    rbmd::Real total_image = blockReduceSum(local_image_sum);

    // 6.
    if (threadIdx.x == 0) {
      rhok_real_final[P_index] = total_real;
      rhok_image_final[P_index] = total_image;
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
      const rbmd::Id P,             // 随机采样波矢数量
      const rbmd::Real alpha,       // Ewald参数
      const rbmd::Real qqrd2e,      // 电荷转换系数
      const rbmd::Real qsqsum,      // 电荷平方和
      const rbmd::Real qsum,        // 电荷总和
      const rbmd::Real* real_array, // ρ(k)实部数组
      const rbmd::Real* imag_array, // ρ(k)虚部数组
      const rbmd::Real* p_sample_x, // 采样波矢x分量
      const rbmd::Real* p_sample_y, // 采样波矢y分量
      const rbmd::Real* p_sample_z, // 采样波矢z分量
      rbmd::Real* virial_tensor,     // 输出维里张量[6]
      rbmd::Real* energy)         // 输出能量
    {

    // 共享内存存放中间结果
    __shared__ rbmd::Real shared_virial[6];
    __shared__ rbmd::Real shared_energy;

    // 初始化共享内存
    if (threadIdx.x < 6) {
        shared_virial[threadIdx.x] = 0.0;
    }
    if (threadIdx.x == 0) {
      shared_energy = 0.0;
    }
    __syncthreads();

    // 计算重要性采样缩放因子S和体积V
    rbmd::Real S, V;
    ComputeS(box, alpha, S);
    V = box._length[0] * box._length[1] * box._length[2];

    // 计算能量和维里的公共系数
    const rbmd::Real energy_coef = (2 * M_PI * S) / (P * V);
    const rbmd::Real virial_coef = - (S / P) * M_PI / V;

    for (rbmd::Id i = threadIdx.x; i < P; i += blockDim.x){
        // 获取采样波矢
        const rbmd::Real m_x = __ldg(&p_sample_x[i]);
        const rbmd::Real m_y = __ldg(&p_sample_y[i]);
        const rbmd::Real m_z = __ldg(&p_sample_z[i]);
        //printf("test---  %f %f  %f ", m_x ,m_y,m_z);

        // 计算波矢k的分量 (2πm/L)
        const rbmd::Real k_x = 2 * M_PI * m_x / box._length[0];
        const rbmd::Real k_y = 2 * M_PI * m_y / box._length[1];
        const rbmd::Real k_z = 2 * M_PI * m_z / box._length[2];

        // 计算波矢模平方|k|²
        const rbmd::Real k2 = k_x * k_x + k_y * k_y + k_z * k_z;

        // 获取结构因子
        const rbmd::Real rho_real = __ldg(&real_array[i]);
        const rbmd::Real rho_imag = __ldg(&imag_array[i]);

        // 计算|ρ(k)|²
        const rbmd::Real rho2 = rho_real * rho_real + rho_imag * rho_imag;

      // ========== 能量计算 ==========
        const rbmd::Real energy_contribution = energy_coef * rho2 / k2;
        atomicAdd(&shared_energy, energy_contribution);

        //printf("test---  %f %f  %f %f", k2 ,rho2, S,V);
        const rbmd::Real base = virial_coef * rho2 / k2;


        // 括号内部系数 (1/4α² + 1/|k|²)
        const rbmd::Real inner_coef = (1.0 / (4 * alpha * alpha)) + (1.0 / k2);

        // 计算张量元素 (δ_βγ - 2k_βk_γ inner_coef)
        const rbmd::Real delta_term = 1.0; // δ_βγ项
        const rbmd::Real kterm_xx = 2 * k_x * k_x * inner_coef;
        const rbmd::Real kterm_yy = 2 * k_y * k_y * inner_coef;
        const rbmd::Real kterm_zz = 2 * k_z * k_z * inner_coef;
        const rbmd::Real kterm_xy = 2 * k_x * k_y * inner_coef;
        const rbmd::Real kterm_xz = 2 * k_x * k_z * inner_coef;
        const rbmd::Real kterm_yz = 2 * k_y * k_z * inner_coef;

        // 原子操作累加到共享内存
        atomicAdd(&shared_virial[0], base * (delta_term - kterm_xx)); // xx
        atomicAdd(&shared_virial[1], base * (delta_term - kterm_yy)); // yy
        atomicAdd(&shared_virial[2], base * (delta_term - kterm_zz)); // zz
        atomicAdd(&shared_virial[3], base * (-kterm_xy));              // xy
        atomicAdd(&shared_virial[4], base * (-kterm_xz));              // xz
        atomicAdd(&shared_virial[5], base * (-kterm_yz));              // yz
    }

    __syncthreads();

    // 第一个线程将结果复制到全局内存
    if (threadIdx.x == 0) {
      // 1. 处理能量结果
      shared_energy -= SQRT(alpha) * qsqsum / MY_PIS + MY_PI2 * qsum * qsum / (alpha * V);
      *energy = shared_energy * qqrd2e;

      // 2. 处理维里张量
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
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputePnumberChargeStructureFactor<<<blocks_per_grid,
                                                       BLOCK_SIZE, 0, 0>>>(
        box, num_atoms, p_number, charge, p_sample_x, p_sample_y, p_sample_z, px,
        py, pz, density_real, density_imag));

    //  p_number
    // const unsigned int blocks_per_grid = p_number;
    // CHECK_KERNEL(ComputePnumberChargeStructureFactor1<<<blocks_per_grid,
    //                                                    BLOCK_SIZE, 0, 0>>>(
    //     box, num_atoms, p_number, charge, p_sample_x, p_sample_y, p_sample_z, px,
    //     py, pz, density_real, density_imag));
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

    CHECK_KERNEL(ComputeRBEVirial<<<1, BLOCK_SIZE, 0, 0>>>(
        box, p_number, alpha, qqrd2e, qsqsum, qsum, real_array, imag_array,
        p_sample_x, p_sample_y, p_sample_z, virial_tensor, energy));
  }

  }  // namespace op