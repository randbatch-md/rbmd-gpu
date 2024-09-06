#include "ljforce_op/ljforce_op.h"
#include <cmath>
#include "box.h"


namespace op
{

#define THREADS_PER_BLOCK 256


	//force kernel 
    //LJForce
	__global__
    void ComputeLJForce(
		    Box* box,
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
			rbmd::Real* force_z,
		    rbmd::Real* evdwl)
	{
		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
	    rbmd::Real sum_fz = 0;
		rbmd::Real sum_eij = 0;

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

				//mix
				rbmd::Real eps_ij = sqrt(eps_i * eps_j);
				rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real z2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;

				//if (molecular_id_i == molecular_id_j)
					//continue; 

				//MinImageDistance(box, px12, py12, pz12);
				rbmd::Real f_ij,e_ij;

				LJForce(cut_off, px12, py12, pz12, eps_ij, sigma_ij,f_ij,e_ij);
				sum_fx += f_ij * px12;
				sum_fy += f_ij * py12;
				sum_fz += f_ij * pz12;

				sum_eij += e_ij;
			}

			// 使用 atomicAdd 保存力和势能，避免数据竞争
			//atomicAdd(&force_x[tid1],sum_fx);
			//atomicAdd(&force_y[tid1],sum_fy);
			//atomicAdd(&force_z[tid1],sum_fz);
			//atomicAdd(&evdwl[tid1], sum_eij);

			force_x[tid1] += sum_fx;
			force_y[tid1] += sum_fy;
			force_z[tid1] += sum_fz;

		    evdwl[tid1] += sum_eij;
		}
	}


	__device__
		LJForce(
			const rbmd::Real cut_off,
			const rbmd::Real px12,
			const rbmd::Real py12,
			const rbmd::Real pz12,
			const rbmd::Real eps_ij,
			const rbmd::Real sigma_ij,
			rbmd::Real& f_ij,
			rbmd::Real& e_ij)
	{
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > small_value)
		{
			rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
			rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;


			 f_ij = -24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 ;
			 e_ij = 0.5 * (4 * eps_ij * (sigmaij_6 / dis_6 - 1) * (sigmaij_6 / dis_6));
		}
	}

	__device__
		rbmd::Real CoulForce(
			const rbmd::Real cut_off,
			const rbmd::Real charge_pi，
			const rbmd::Real charge_pj，
			const rbmd::Real px12,
			const rbmd::Real py12,
			const rbmd::Real pz12)
	{
		rbmd::Real MY_pi = 3.14159265358979323846; //pi
		rbmd::Real MY_pis = 1.77245385090551602729; // sqrt(pi)

		rbmd::Real CoulForce = 0;
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real dis = sqrt(dis_2);
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > small_value)
		{
			rbmd::Real erfcx = sqrt(_alpha) * dis;
			rbmd::Real expx = -_alpha * dis * dis;
			rbmd::Real Gnearvalue = (1.0 - erf(erfcx)) / (dis * dis) +
				2 * sqrt(_alpha) * exp(expx) / (MY_pis * dis);

			CoulForce = -charge_pi * charge_pj * Gnearvalue / dis;


		}
		return CoulForce;
	}


	//LJVirial
	__global__
		void ComputeJVirial(
			Box* box,
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

				//mix
				rbmd::Real eps_ij = sqrt(eps_i * eps_j);
				rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real x2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;

				if (molecular_id_i == molecular_id_j)
					continue;

				MinMirror(box, px12, py12, pz12);
				rbmd::Real Virial_f;
				Virial_f = LJVirial(cut_off, px12, py12, pz12, eps_ij, sigma_ij);

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
			const rbmd::Real eps_ij,
			const rbmd::Real sigma_ij)
	{
		rbmd::Real virial_f = 0;;

		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > small_value)
		{
			rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
			rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

			 virial_f = 0.5 * 24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;

	}
		return virial_f;
	}


	//
	struct LJforceOp<device::DEVICE_GPU>;
	{
		void operator()(
			Box* box,
			const rbmd::Id& num_atoms,
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
			rbmd::Real* force_z
			rbmd::Real* evdwl)
		{


			int block = (N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
			hipLaunchKernelGGL(HIP_KERNEL_NAME(ComputeLJForce<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK),
				box,num_atoms, atoms_type, molecular_type ,sigma, eps, start_id, end_id, id_verletlist,px, py, pz, force_x, force_y, force_z,evdwl);


			hipErrorCheck(hipGetLastError());
			hipErrorCheck(hipDeviceSynchronize());
		}
	};
	template struct LJforceOp<device::DEVICE_GPU>;

}

