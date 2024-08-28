#include "ljforce_op/ljforce_op.h"
#include "base/rocm.h"

namespace op
{

#define THREADS_PER_BLOCK 256


	//force kernel 
    //LJForce
	__global__
    void ComputeLJForce(
			const Box& box,
		    const rbmd::Id& num_atoms,
		    const rbmd::Id* atoms_type,
		    const rbmd::Id* molecular_type,
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
		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
	    rbmd::Real sum_fz = 0;

		int tid1 = threadIdx.x + blockIdx.x * blockDim.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1];
			rbmd::Id molecular_id_i=  molecular_type[tid1];
			rbmd::Real eps_i = eps[typei];
			rbmd::Real sigma_i = sigma[typei];

			rbmd::Real x1 = px[tid1];
			rbmd::Real y1 = py[tid1];
			rbmd::Real z1 = pz[tid1];


			for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
			{
				rbmd::Id tid2 = id_verletlist[j];

				rbmd::Id typej = atoms_type[tid2];
				rbmd::Id molecular_id_j = molecular_type[tid2];
				rbmd::Real eps_j = eps[typej];
				rbmd::Real sigma_j = sigma[typej];

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real x2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;

				if (molecular_id_i == molecular_id_j)
					return;

				//MinMirror(box, px12, py12, pz12);
				rbmd::Real f;
				f = LJForce(cut_off, px12, py12, pz12, eps_i, eps_j, sigma_i, sigma_j);
				sum_fx += f * px12;
				sum_fy += f * py12;
				sum_fz += f * pz12;

			}

			// save force
			//atomicAdd(&force_x[tid1],sum_fx);
			//atomicAdd(&force_y[tid1],sum_fy);
			//atomicAdd(&force_z[tid1],sum_fz);
			force_x[tid1] += sum_fx;
			force_y[tid1] += sum_fy;
			force_z[tid1] += sum_fz;

		}
	}


	__device__
		rbmd::Real LJForce(
			const rbmd::Real cut_off,
			const rbmd::Real px12,
			const rbmd::Real py12,
			const rbmd::Real pz12,
			const rbmd::Real eps_i,
			const rbmd::Real eps_j,
			const rbmd::Real sigma_i,
			const rbmd::Real sigma_j)
	{
		rbmd::Real f = 0;
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > small_value)
		{
			rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

			rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
			rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;
			 f = -24 * Sqrt(eps_i * eps_j) * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 ;
		}
		return f;
	}


	//LJVirial
	__global__
		void ComputeJVirial(
			const Box& box,
			const rbmd::Id& num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Id* molecular_type,
			const rbmd::Real* eps,
			const rbmd::Real* sigma,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* virial_xx,
			rbmd::Real* virial_yy,
			rbmd::Real* virial_zz,
			rbmd::Real* virial_xy,
			rbmd::Real* virial_xz,
			rbmd::Real* virial_yz)
	{

		rbmd::Real sum_virial_xx = 0;
		rbmd::Real sum_virial_yy= 0;
		rbmd::Real sum_virial_zz = 0;
		rbmd::Real sum_virial_xy = 0;
		rbmd::Real sum_virial_xz = 0;
		rbmd::Real sum_virial_yz = 0;

		int tid1 = threadIdx.x + blockIdx.x * blockDim.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1];
			rbmd::Id molecular_id_i = molecular_type[tid1];
			rbmd::Real eps_i = eps[typei];
			rbmd::Real sigma_i = sigma[typei];

			rbmd::Real x1 = px[tid1];
			rbmd::Real y1 = py[tid1];
			rbmd::Real z1 = pz[tid1];


			for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
			{
				rbmd::Id tid2 = id_verletlist[j];

				rbmd::Id typej = atoms_type[tid2];
				rbmd::Id molecular_id_j = molecular_type[tid2];
				rbmd::Real eps_j = eps[typej];
				rbmd::Real sigma_j = sigma[typej];

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real x2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;

				if (molecular_id_i == molecular_id_j)
					return;

				//MinMirror(box, px12, py12, pz12);
				rbmd::Real Virial_f;
				Virial_f = LJVirial(cut_off, px12, py12, pz12, eps_i, eps_j, sigma_i, sigma_j);

				rbmd::Real Virial_fx = Virial_f * px12;
				rbmd::Real Virial_fy = Virial_f * py12;
				rbmd::Real Virial_fz = Virial_f * pz12;

				sum_virial_xx +=  px12 * Virial_fx;
				sum_virial_yy +=  py12 * Virial_fy;
				sum_virial_zz +=  pz12 * Virial_fz;
				sum_virial_xy +=  px12 * Virial_fy;
				sum_virial_xz +=  px12 * Virial_fz;
				sum_virial_yz +=  py12 * Virial_fz;
			}

			// save virial
			virial_xx[tid1] += sum_virial_xx;
			virial_yy[tid1] += sum_virial_yy;
			virial_zz[tid1] += sum_virial_zz;
			virial_xy[tid1] += sum_virial_xy;
			virial_xz[tid1] += sum_virial_xz;
			virial_yz[tid1] += sum_virial_yz;

		}
	}
	__device__
		rbmd::Real LJVirial(
			const rbmd::Real cut_off,
			const rbmd::Real px12,
			const rbmd::Real py12,
			const rbmd::Real pz12,
			const rbmd::Real eps_i,
			const rbmd::Real eps_j,
			const rbmd::Real sigma_i,
			const rbmd::Real sigma_j)
	{
		rbmd::Real virial_f = 0;;

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

			 virial_f = 0.5 * 24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;

	}
		return virial_f;
	}


	//
	struct LJforceOp<device::DEVICE_GPU>;
	{
		void operator()(
			Box& box,
			rbmd::Id& num_atoms,
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
			rbmd::Real* force_x,
			rbmd::Real* force_y,
			rbmd::Real* force_z)
		{


			int block = (N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
			hipLaunchKernelGGL(HIP_KERNEL_NAME(ComputeLJForce<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK),
				box,num_atoms, atoms_type, molecular_type ,sigma, eps, start_id, end_id, id_verletlist,px, py, pz, force_x, force_y, force_z);


			hipErrorCheck(hipGetLastError());
			hipErrorCheck(hipDeviceSynchronize());
		}
	};
	template struct LJforceOp<device::DEVICE_GPU>;

}

