#include "ljforce_op/ljforce_op.h"
#include "base/rocm.h"

namespace op
{

#define THREADS_PER_BLOCK 256


	//force kernel 
	__global__
    void ComputeLJForce(
			const Box& box,
		    const rbmd::Id& num_atoms,
		    const rbmd::Id* atoms_type,
		    const rbmd::Real* eps,
		    const rbmd::Real* sigma,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* force_x,
			rbmd::Real* force_y,
			rbmd::Real* force_z)
	{
		rbmd::Real3 sum_f{ 0,0,0 };

		int tid1 = threadIdx.x + blockIdx.x * blockDim.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1];
			rbmd::Real eps_i = eps[typei];
			rbmd::Real sigma_i = sigma[typei];

			rbmd::Real x1 = px[tid1];
			rbmd::Real y1 = py[tid1];
			rbmd::Real z1 = pz[tid1];


			for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
			{
				rbmd::Id tid2 = id_verletlist[j];
				rbmd::Id typej = atoms_type[tid2];
				rbmd::Real eps_j = eps[typej];
				rbmd::Real sigma_j = sigma[typej];

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real x2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;

				//MinMirror(box, px12, py12, pz12);
				sum_f += LJForce(cut_off, px12, py12, pz12, eps_i, eps_j, sigma_i, sigma_j);
			}

			// save force
			atomicAdd(&force_x[tid1], sum_f[0]);
			atomicAdd(&force_y[tid1], sum_f[1]);
			atomicAdd(&force_z[tid1], sum_f[2]);
		}
	}


	//LJForce
	__device__
		rbmd::Real3 LJForce(
			const rbmd::Real cut_off,
			const rbmd::Real px12,
			const rbmd::Real py12,
			const rbmd::Real pz12,
			const rbmd::Real eps_i,
			const rbmd::Real eps_j,
			const rbmd::Real sigma_i,
			const rbmd::Real sigma_j)
	{
		rbmd::Real3 f{ 0, 0, 0 };

		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > small_value)
		{
			rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

			rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
			rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;
			rbmd::Real fx = -24 * Sqrt(eps_i * eps_j) * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * px12;
			rbmd::Real fy = -24 * Sqrt(eps_i * eps_j) * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * py12;
			rbmd::Real fz = -24 * Sqrt(eps_i * eps_j) * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * pz12;
			f{ fx,fy,fz };
		}
		return f;
	}

	//LJVirial
	__device__
		rbmd::Real6 LJVirial(
			const rbmd::Real cut_off,
			const rbmd::Real px12,
			const rbmd::Real py12,
			const rbmd::Real pz12,
			const rbmd::Real eps_i,
			const rbmd::Real eps_j,
			const rbmd::Real sigma_i,
			const rbmd::Real sigma_j)
	{
		rbmd::Real6 virial{ 0, 0, 0, 0, 0, 0 };

		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > small_value)
		{
			rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

			rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
			rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;
			rbmd::Real eps_ij = Sqrt(eps_i * eps_j);

			rbmd::Real fx = 0.5 * 24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * px12;
			rbmd::Real fy = 0.5 * 24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * py12;
			rbmd::Real fz = 0.5 * 24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * pz12;

			//compute virial
			virial[0] = px12 * fx; //xx
			virial[1] = py12 * fy; //yy
			virial[2] = pz12 * fz; //zz
			virial[3] = px12 * fy; //xy
			virial[4] = px12 * fz; //xz
			virial[5] = py12 * fz; //yz
		}
		return virial;
	}


	//
	struct LJforceOp<device::DEVICE_GPU>;
	{
		void operator()(
			Box& box,
			rbmd::Id& num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* force_x,
			rbmd::Real* force_y,
			rbmd::Real* force_z)
		{


			int block = (N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
			hipLaunchKernelGGL(HIP_KERNEL_NAME(ComputeLJForce<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK),
				box,num_atoms, atoms_type,sigma, eps, start_id, end_id, id_verletlist,px, py, pz, force_x, force_y, force_z);


			hipErrorCheck(hipGetLastError());
			hipErrorCheck(hipDeviceSynchronize());
		}
	};
	template struct LJforceOp<device::DEVICE_GPU>;

}

