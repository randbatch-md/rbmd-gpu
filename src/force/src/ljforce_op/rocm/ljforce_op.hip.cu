#include "ljforce_op.h"
#include "model/box.h"
#include <hip/hip_runtime.h>
#include "../common/rbmd_define.h"

namespace op
{

	//force kernel 
    //LJForce
    __device__ void lj126(
		    rbmd::Real cut_off,
		    rbmd::Real px12,
		    rbmd::Real py12,
		    rbmd::Real pz12,
		    rbmd::Real eps_ij,
		    rbmd::Real sigma_ij,
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

	__device__ void lj126_rs(
		rbmd::Real rs,
		rbmd::Real px12,
		rbmd::Real py12,
		rbmd::Real pz12,
		rbmd::Real eps_ij,
		rbmd::Real sigma_ij,
		rbmd::Real& fs_ij)
	{
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real rs_2 = rs * rs;

		if (dis_2 < rs_2 && dis_2 > small_value)
		{
			rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
			rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

			fs_ij = -24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;
		}
	}

	__device__ void lj126_rcs(
		rbmd::Real rc,
		rbmd::Real rs,
		rbmd::Id pice_num,
		rbmd::Real px12,
		rbmd::Real py12,
		rbmd::Real pz12,
		rbmd::Real eps_ij,
		rbmd::Real sigma_ij,
		rbmd::Real& fcs_ij)
	{
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real rc_2 = rc * rc;
		const rbmd::Real rs_2 = rs * rs;

		if (dis_2 < rc_2 && dis_2 > rs_2)
		{
			rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
			rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

			fcs_ij = pice_num  * (-24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2);
		}
	}

	__global__ void ComputeLJForce(
						Box* box,
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
			            rbmd::Real* force_x,
			            rbmd::Real* force_y,
			            rbmd::Real* force_z,
			            rbmd::Real* evdwl,
					    rbmd::Real* total_evdwl)
	{
		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
	    rbmd::Real sum_fz = 0;
		rbmd::Real sum_eij = 0;

		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
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

				MinImageDistance(box, px12, py12, pz12);

				rbmd::Real f_ij;
				rbmd::Real e_ij;

				lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij,f_ij,e_ij);
				sum_fx += f_ij * px12;
				sum_fy += f_ij * py12;
				sum_fz += f_ij * pz12;
				sum_eij += e_ij;
			}

			force_x[tid1] = sum_fx;
			force_y[tid1] = sum_fy;
			force_z[tid1] = sum_fz;

		    evdwl[tid1] += sum_eij;
			atomicAdd(total_evdwl, sum_eij);
			//printf("--------test---evdwl[tid1]:%f---\n",evdwl[tid1]);
		}
		
	}

	__global__ void ComputeLJRBLForce(
		Box* box,
		const rbmd::Real rs,
		const rbmd::Real rc,
		const rbmd::Id pice_num,
		const rbmd::Id num_atoms,
		const rbmd::Id* atoms_type,
		const rbmd::Id* molecular_type,
		const rbmd::Real* sigma,
		const rbmd::Real* eps,
		const rbmd::Id* start_id,
		const rbmd::Id* end_id,
		const rbmd::Id* id_verletlist,
		const rbmd::Id* id_verletlist_rcs,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* corr_force_x,
		rbmd::Real* corr_force_y,
		rbmd::Real* corr_force_z,
		rbmd::Real* corr_value_x,
		rbmd::Real* corr_value_y,
		rbmd::Real* corr_value_z)
	{
		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
		rbmd::Real sum_fz = 0;
		rbmd::Real sum_eij = 0;

		rbmd::Real sum_fsx = 0;
		rbmd::Real sum_fsy = 0;
		rbmd::Real sum_fsz = 0;

		rbmd::Real sum_fcsx = 0;
		rbmd::Real sum_fcsy = 0;
		rbmd::Real sum_fcsz = 0;

		rbmd::Real  corr_fx_temp = 0;
		rbmd::Real  corr_fy_temp = 0;
		rbmd::Real  corr_fz_temp = 0;


		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1];
			rbmd::Id molecular_id_i = molecular_type[tid1];
			rbmd::Real eps_i = eps[typei];
			rbmd::Real sigma_i = sigma[typei];

			rbmd::Real x1 = px[tid1];
			rbmd::Real y1 = py[tid1];
			rbmd::Real z1 = pz[tid1];

			rbmd::Real fs_ij, fcs_ij;


			//rs 
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

				MinImageDistance(box, px12, py12, pz12);

				lj126_rs(rs, px12, py12, pz12, eps_ij, sigma_ij, fs_ij);			// the force of rs 

				sum_fsx += fs_ij * px12;
				sum_fsy += fs_ij * py12;
				sum_fsz += fs_ij * pz12;
			}

			//rcs
			for (int j = id_verletlist_rcs[tid1* pice_num]; j < pice_num; ++j)
			{

				rbmd::Id tid2 = id_verletlist_rcs[j];
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

				MinImageDistance(box, px12, py12, pz12);

				lj126_rcs(rc, rs, pice_num,px12, py12, pz12, eps_ij, sigma_ij, fcs_ij); //the force  rs to rc

				sum_fcsx += fcs_ij * px12;
				sum_fcsy += fcs_ij * py12;
				sum_fcsz += fcs_ij * pz12;
			}

			//sum
			sum_fx = sum_fsx + sum_fcsx;
			sum_fy = sum_fsy + sum_fcsy;
			sum_fz = sum_fsz + sum_fcsz;
			
			corr_force_x[tid1] = sum_fx;
			corr_force_y[tid1] = sum_fy;
			corr_force_z[tid1] = sum_fz;

			corr_fx_temp += corr_force_x[tid1];
			corr_fy_temp += corr_force_y[tid1];
			corr_fz_temp += corr_force_z[tid1];

			corr_force_x[tid1] = corr_force_x[tid1] - corr_fx_temp;
			corr_force_y[tid1] = corr_force_y[tid1] - corr_fy_temp;
			corr_force_z[tid1] = corr_force_z[tid1] - corr_fz_temp;


			//atomicAdd(corr_value_x, corr_fx_temp);
			//atomicAdd(corr_value_y, corr_fy_temp);
			//atomicAdd(corr_value_z, corr_fz_temp);
			//printf("--------test---evdwl[tid1]:%f---\n",evdwl[tid1]);
		}

	}

	__device__
		rbmd::Real CoulForce(
			 rbmd::Real cut_off,
			 rbmd::Real alpha,
			 rbmd::Real charge_pi,
			 rbmd::Real charge_pj,
			 rbmd::Real px12,
			 rbmd::Real py12,
			 rbmd::Real pz12)
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
			rbmd::Real erfcx = sqrt(alpha) * dis;
			rbmd::Real expx = -alpha * dis * dis;
			rbmd::Real Gnearvalue = (1.0 - erf(erfcx)) / (dis * dis) +
				2 * sqrt(alpha) * exp(expx) / (MY_pis * dis);

			CoulForce = -charge_pi * charge_pj * Gnearvalue / dis;


		}
		return CoulForce;
	}

	__device__ void EwaldForce(
		Box* box,
		const rbmd::Real alpha
		float3 M,
		float2 rhok_ri,
		const rbmd::Real* charge,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real ewald_force_x,
		rbmd::Real ewald_force_y
		rbmd::Real ewald_force_z)
	{
		rbmd::Real ewald_force;
		rbmd::Real volume = box->_length[0] * box->_length[1] * box->_length[2];

		float3 k{ 0,0,0 };
		K[0] = 2 * M_PI * M[0] / box->_length[0];
		K[1] = 2 * M_PI * M[1] / box->_length[1];
		K[2] = 2 * M_PI * M[2] / box->_length[2];

		rbmd::Real range_K_2 = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];
		rbmd::Real dot_product = K[0] * px + K[1] * py + K[2] * pz;

		rbmd::Real factor_a = -4 * M_PI * charge;
		rbmd::Real factor_b = exp(-range_K_2 / (4 * alpha));
		rbmd::Real factor_c = cos(dot_product) * rhok_ri[1];
		rbmd::Real factor_d = sin(dot_product) * rhok_ri[0];

		ewald_force = factor_a / (volume * range_K_2) * factor_b * (factor_c - factor_d);
		ewald_force_x = ewald_force * K[0];
		ewald_force_y = ewald_force * K[1];
		ewald_force_z = ewald_force * K[2];
	}


	__global__ void ComputeChargeStructureFactorComponent(
		const rbmd::Id num_atoms,
		float3 K,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		const rbmd::Real* charge,
		rbmd::Real* density_real,
		rbmd::Real* density_imag)
	{
		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Real local_charge = charge[tid1];
			rbmd::Real dot_product = K.x * px[tid1] + K.y * py[tid1] + K.z * pz[tid1];

			density_real[tid1] = local_charge * cos(dot_product);
			density_imag[tid1] = local_charge * sin(dot_product);
		}
	}

	__global__ void ComputeEwaldForce(
		Box* box,
		const rbmd::Id num_atoms,
	    const rbmd::Id  Kmax,
		const rbmd::Real alpha,
	    const float2  whole_rhok,
		const rbmd::Real* charge,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* ewald_force_x,
		rbmd::Real* ewald_force_y,
		rbmd::Real* ewald_force_z)
	{
		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;

		if (tid1 < num_atoms)
		{
			rbmd::Real p_x = px[tid1];
			rbmd::Real p_y = py[tid1];
			rbmd::Real p_z = pz[tid1];
			charge_i = charge[tid1];
			rbmd::Real force_x, force_y, force_z;

			for (rbmd::Id i = -Kmax; i <= Kmax; ++i)
			{
				for (rbmd::Id j = -Kmax; j <= Kmax; ++j)
				{
					for (rbmd::Id k = -Kmax; k <= Kmax; ++k)
					{
						if (!(i == 0 && j == 0 && k == 0)) 
						{
							float3 M = { (rbmd::Real)i, (rbmd::Real)j, (rbmd::Real)k };
							rbmd::Id indexEwald = (i + Kmax) * (2 * Kmax + 1) * 
												  (2 * Kmax + 1) + (j + Kmax) *
												  (2 * Kmax + 1) + (k + Kmax);
							 float2 rhok_i = whole_rhok[indexEwald];

							 EwaldForce(box, alpha,M, rhok_i, charge_i, p_x, p_y, p_z,
										force_x, force_y, force_z);
						}
					}
				}
			}
			ewald_force_x[tid1] = force_x;
			ewald_force_y[tid1] = force_y;
			ewald_force_z[tid1] = force_z;
	    }
	}


	//LJVirial
	//__global__ void ComputeJVirial(
	//		Box* box,
	//	    const rbmd::Real cut_off,
	//		const rbmd::Id& num_atoms,
	//		const rbmd::Id* atoms_type,
	//		const rbmd::Id* molecular_type,
	//		const rbmd::Real* eps,
	//		const rbmd::Real* sigma,
	//		const rbmd::Id* start_id,
	//		const rbmd::Id* end_id,
	//		const rbmd::Id* id_verletlist,
	//		const rbmd::Real* px,
	//		const rbmd::Real* py,
	//		const rbmd::Real* pz,
	//		rbmd::Real* virial_xx,
	//		rbmd::Real* virial_yy,
	//		rbmd::Real* virial_zz,
	//		rbmd::Real* virial_xy,
	//		rbmd::Real* virial_xz,
	//		rbmd::Real* virial_yz)
	//{

	//	rbmd::Real sum_virial_xx = 0;
	//	rbmd::Real sum_virial_yy= 0;
	//	rbmd::Real sum_virial_zz = 0;
	//	rbmd::Real sum_virial_xy = 0;
	//	rbmd::Real sum_virial_xz = 0;
	//	rbmd::Real sum_virial_yz = 0;

	//	int tid1 = threadIdx.x + blockIdx.x * blockDim.x;
	//	if (tid1 < num_atoms)
	//	{
	//		rbmd::Id typei = atoms_type[tid1];
	//		rbmd::Id molecular_id_i = molecular_type[tid1];
	//		rbmd::Real eps_i = eps[typei];
	//		rbmd::Real sigma_i = sigma[typei];

	//		rbmd::Real x1 = px[tid1];
	//		rbmd::Real y1 = py[tid1];
	//		rbmd::Real z1 = pz[tid1];


	//		for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
	//		{
	//			rbmd::Id tid2 = id_verletlist[j];

	//			rbmd::Id typej = atoms_type[tid2];
	//			rbmd::Id molecular_id_j = molecular_type[tid2];
	//			rbmd::Real eps_j = eps[typej];
	//			rbmd::Real sigma_j = sigma[typej];

	//			//mix
	//			rbmd::Real eps_ij = sqrt(eps_i * eps_j);
	//			rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

	//			rbmd::Real x2 = px[tid2];
	//			rbmd::Real y2 = py[tid2];
	//			rbmd::Real z2 = pz[tid2];
	//			rbmd::Real px12 = x2 - x1;
	//			rbmd::Real py12 = y2 - y1;
	//			rbmd::Real pz12 = z2 - z1;

	//			if (molecular_id_i == molecular_id_j)
	//				continue;

	//			//MinMirror(box, px12, py12, pz12);
	//			rbmd::Real Virial_f;
	//			Virial_f = LJVirial(cut_off, px12, py12, pz12, eps_ij, sigma_ij);

	//			rbmd::Real Virial_fx = Virial_f * px12;
	//			rbmd::Real Virial_fy = Virial_f * py12;
	//			rbmd::Real Virial_fz = Virial_f * pz12;

	//			sum_virial_xx +=  px12 * Virial_fx;
	//			sum_virial_yy +=  py12 * Virial_fy;
	//			sum_virial_zz +=  pz12 * Virial_fz;
	//			sum_virial_xy +=  px12 * Virial_fy;
	//			sum_virial_xz +=  px12 * Virial_fz;
	//			sum_virial_yz +=  py12 * Virial_fz;
	//		}

	//		// save virial
	//		virial_xx[tid1] += sum_virial_xx;
	//		virial_yy[tid1] += sum_virial_yy;
	//		virial_zz[tid1] += sum_virial_zz;
	//		virial_xy[tid1] += sum_virial_xy;
	//		virial_xz[tid1] += sum_virial_xz;
	//		virial_yz[tid1] += sum_virial_yz;

	//	}
	//}

	//__device__ rbmd::Real LJVirial(
	//		const rbmd::Real cut_off,
	//		const rbmd::Real px12,
	//		const rbmd::Real py12,
	//		const rbmd::Real pz12,
	//		const rbmd::Real eps_ij,
	//		const rbmd::Real sigma_ij)
	//{
	//	rbmd::Real virial_f = 0;;

	//	const rbmd::Real  small_value = 0.0001;
	//	const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
	//	const rbmd::Real cut_off_2 = cut_off * cut_off;

	//	if (dis_2 < cut_off_2 && dis_2 > small_value)
	//	{
	//		rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
	//		rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
	//		rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

	//		 virial_f = 0.5 * 24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;

	//}
	//	return virial_f;
	//}


	//
	void LJForceOp<device::DEVICE_GPU>::operator()(
						Box* box,
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
			            rbmd::Real* force_x,
			            rbmd::Real* force_y,
			            rbmd::Real* force_z,
			            rbmd::Real* evdwl,
						rbmd::Real* total_evdwl)
		{
		    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		    CHECK_KERNEL(ComputeLJForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (box, cut_off, num_atoms, atoms_type, molecular_type,
				sigma, eps, start_id, end_id, id_verletlist, px, py, pz, force_x, force_y, force_z, evdwl,total_evdwl));
		}

	//RBL
	void LJRBLForceOp<device::DEVICE_GPU>::operator()(
		Box* box,
		const rbmd::Real rs,
		const rbmd::Real rc,
		const rbmd::Id pice_num,
		const rbmd::Id num_atoms,
		const rbmd::Id* atoms_type,
		const rbmd::Id* molecular_type,
		const rbmd::Real* sigma,
		const rbmd::Real* eps,
		const rbmd::Id* start_id,
		const rbmd::Id* end_id,
		const rbmd::Id* id_verletlist,
		const rbmd::Id* id_verletlist_rcs,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* corr_force_x,
		rbmd::Real* corr_force_y,
		rbmd::Real* corr_force_z,
		rbmd::Real* corr_value_x,
		rbmd::Real* corr_value_y,
		rbmd::Real* corr_value_z )
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		CHECK_KERNEL(ComputeLJRBLForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (box, rs, rc, pice_num, num_atoms, atoms_type, molecular_type,
			sigma, eps, start_id, end_id, id_verletlist, id_verletlist_rcs, px, py, pz, corr_force_x, corr_force_y, corr_force_z, corr_value_x, corr_value_y, corr_value_z));
	}

	void ComputeChargeStructureFactorComponentOp<device::DEVICE_GPU>::operator()(
		const rbmd::Id num_atoms,
		const float3& k,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		const rbmd::Real* charge,
		rbmd::Real* density_real,
		rbmd::Real* density_imag)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
		
		CHECK_KERNEL(ComputeChargeStructureFactorComponent <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
			(num_atoms, k, px, py ,pz, charge, density_real, density_imag));

	}

	void ComputeEwaldForceOp<device::DEVICE_GPU>::operator()(
		Box* box,
		const rbmd::Id num_atoms,
		const rbmd::Id  Kmax,
		const rbmd::Real alpha,
		const float2  whole_rhok,
		const rbmd::Real* charge,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* ewald_force_x,
		rbmd::Real* ewald_force_y,
		rbmd::Real* ewald_force_z)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		CHECK_KERNEL(ComputeEwaldForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
			(box,num_atoms, Kmax, alpha, whole_rhok, charge,
		     px, py, pz, ewald_force_x, ewald_force_y, ewald_force_z));

	}

}

