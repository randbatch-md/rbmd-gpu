
#include "../common/rbmd_define.h"
#include "tersoff_op.h"
#include "tersoff.h"
#include "model/box.h"

rbmd::Id maxshort = 10;

namespace op{
//---------Math---------//
  __device__ void scale3(const rbmd::Real s, rbmd::Real *v)
  {
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
  }

  __device__ void scale3(const rbmd::Real s, const rbmd::Real *v, rbmd::Real *ans)
  {
    ans[0] = s * v[0];
    ans[1] = s * v[1];
    ans[2] = s * v[2];
  }

  __device__ rbmd::Real dot3(const rbmd::Real *v1, const rbmd::Real *v2)
  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  }

  __device__ void scaleadd3(const rbmd::Real s, const rbmd::Real *v1, const rbmd::Real *v2, rbmd::Real *ans)
  {
    ans[0] = s * v1[0] + v2[0];
    ans[1] = s * v1[1] + v2[1];
    ans[2] = s * v1[2] + v2[2];
  }

  __device__ void  add3(const rbmd::Real *v1, const rbmd::Real *v2, rbmd::Real *ans)
  {
    ans[0] = v1[0] + v2[0];
    ans[1] = v1[1] + v2[1];
    ans[2] = v1[2] + v2[2];
  }

//---------device---------//
  // inlined functions for efficiency
  __device__ rbmd::Real device_gijk(const rbmd::Real costheta, const TersoffParams *const  TersoffParams)
    {
      const rbmd::Real ters_c = TersoffParams->c * TersoffParams->c;
      const rbmd::Real ters_d = TersoffParams->d * TersoffParams->d;
      const rbmd::Real hcth = TersoffParams->costheta0 - costheta;

      return TersoffParams->gamma * (1.0 + ters_c / ters_d - ters_c / (ters_d + hcth * hcth));
    }

  __device__ rbmd::Real device_gijk_d(const rbmd::Real costheta, const TersoffParams *const TersoffParams)
    {
      const rbmd::Real ters_c = TersoffParams->c * TersoffParams->c;
      const rbmd::Real ters_d = TersoffParams->d * TersoffParams->d;
      const rbmd::Real hcth = TersoffParams->costheta0 - costheta;
      const rbmd::Real numerator = -2.0 * ters_c * hcth;
      const rbmd::Real denominator = 1.0 / (ters_d + hcth * hcth);
      return TersoffParams->gamma * numerator * denominator * denominator;
    }

  //fc
  __device__ rbmd::Real device_fc(rbmd::Real r, TersoffParams *TersoffParams)
  {
    rbmd::Real ters_R = TersoffParams->R;
    rbmd::Real ters_D = TersoffParams->D;

    if (r < ters_R-ters_D) return 1.0;
    if (r > ters_R+ters_D) return 0.0;
    return 0.5*(1.0 - SIN(M_PI_2*(r - ters_R)/ters_D));
  }

  __device__ rbmd::Real device_fc_d(rbmd::Real r, TersoffParams* TersoffParams)
    {
      rbmd::Real ters_R = TersoffParams->R;
      rbmd::Real ters_D = TersoffParams->D;

      if (r < ters_R-ters_D) return 0.0;
      if (r > ters_R+ters_D) return 0.0;
      return -(M_PI_4/ters_D) * COS(M_PI_2*(r - ters_R)/ters_D);
    }

  //  //fa
  __device__ rbmd::Real device_fa(rbmd::Real r, TersoffParams* TersoffParams)
    {
      if (r > TersoffParams->R + TersoffParams->D) return 0.0;
      return -TersoffParams->B * exp(-TersoffParams->lambda2 * r) * device_fc(r,TersoffParams);
    }

  __device__ rbmd::Real device_fa_d(rbmd::Real r, TersoffParams* TersoffParams)
    {
      if (r > TersoffParams->R + TersoffParams->D) return 0.0;
      return TersoffParams->B * exp(-TersoffParams->lambda2 * r) *
        (TersoffParams->lambda2 * device_fc(r,TersoffParams) - device_fc_d(r,TersoffParams));
    }


// //bij
__device__ rbmd::Real device_bij(rbmd::Real zeta, TersoffParams* TersoffParams)
  {
    rbmd::Real tmp = TersoffParams->beta * zeta;
    if (tmp > TersoffParams->c1) return 1.0/sqrt(tmp);
    if (tmp > TersoffParams->c2)
      return (1.0 - pow(tmp,-TersoffParams->n) / (2.0*TersoffParams->n))/sqrt(tmp);
    if (tmp < TersoffParams->c4) return 1.0;
    if (tmp < TersoffParams->c3)
      return 1.0 - pow(tmp,TersoffParams->n)/(2.0*TersoffParams->n);
    return pow(1.0 + pow(tmp,TersoffParams->n), -1.0/(2.0*TersoffParams->n));
  }

//
__device__ rbmd::Real device_bij_d(rbmd::Real zeta, TersoffParams* TersoffParams)
  {
    rbmd::Real tmp = TersoffParams->beta * zeta;
    if (tmp > TersoffParams->c1) return TersoffParams->beta * -0.5*pow(tmp,-1.5);
    if (tmp > TersoffParams->c2)
      return TersoffParams->beta * (-0.5*pow(tmp,-1.5) *
                            // error in negligible 2nd term fixed 9/30/2015
                            // (1.0 - 0.5*(1.0 +  1.0/(2.0*TersoffParams->powern)) *
                            (1.0 - (1.0 +  1.0/(2.0*TersoffParams->n)) *
                             pow(tmp,-TersoffParams->n)));
    if (tmp < TersoffParams->c4) return 0.0;
    if (tmp < TersoffParams->c3)
      return -0.5*TersoffParams->beta * pow(tmp,TersoffParams->n-1.0);

    rbmd::Real tmp_n = pow(tmp,TersoffParams->n);
    return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*TersoffParams->n)))*tmp_n / zeta;
  }

  // // zeta
  __device__ rbmd::Real zeta(TersoffParams* TersoffParams, rbmd::Real rsqij, rbmd::Real rsqik,
                           rbmd::Real *rij_hat, rbmd::Real *rik_hat)
    {
      rbmd::Real rij,rik,costheta,arg,ex_delr;

      rij = SQRT(rsqij);
      rik = SQRT(rsqik);
      costheta = dot3(rij_hat,rik_hat);

      if (TersoffParams->m_int == 3) arg = POW(TersoffParams->lambda3 * (rij-rik),3.0);
      else arg = TersoffParams->lambda3 * (rij-rik);

      if (arg > 69.0776) ex_delr = 1.e30;
      else if (arg < -69.0776) ex_delr = 0.0;
      else ex_delr = EXP(arg);

      return device_fc(rik,TersoffParams) * device_gijk(costheta,TersoffParams) * ex_delr;
    }

 __device__ void force_zeta(TersoffParams* TersoffParams, rbmd::Real rsq, rbmd::Real zeta_ij,
                             rbmd::Real &fforce, rbmd::Real &prefactor,rbmd::Real &eng)
  {
    rbmd::Real r,fa,fa_d,bij;

    r = SQRT(rsq);
    fa = device_fa(r,TersoffParams);
    fa_d = device_fa_d(r,TersoffParams);
    bij = device_bij(zeta_ij,TersoffParams);
    fforce = 0.5*bij*fa_d;
    prefactor = -0.5*fa * device_bij_d(zeta_ij,TersoffParams);
    eng = 0.5*bij*fa;
  }


__device__ void repulsive(TersoffParams* TersoffParams,rbmd::Real rsq, rbmd::Real &fforce,rbmd::Real &eng)
  {
    rbmd::Real r,tmp_fc,tmp_fc_d,tmp_exp;

    r = SQRT(rsq);
    tmp_fc = device_fc(r,TersoffParams);
    tmp_fc_d = device_fc_d(r,TersoffParams);
    tmp_exp = EXP(-TersoffParams->lambda1 * r);
    fforce = -TersoffParams->A * tmp_exp * (tmp_fc_d - tmp_fc*TersoffParams->lambda1) / r;
    eng = tmp_fc * TersoffParams->A * tmp_exp;
  }

__device__ void costheta_d(rbmd::Real *rij_hat, rbmd::Real rijinv,
                             rbmd::Real *rik_hat, rbmd::Real rikinv,
                             rbmd::Real *dri, rbmd::Real *drj, rbmd::Real *drk)
  {
    // first element is devative wrt Ri, second wrt Rj, third wrt Rk

    rbmd::Real cos_theta = dot3(rij_hat,rik_hat);

    scaleadd3(-cos_theta,rij_hat,rik_hat,drj);
    scale3(rijinv,drj);
    scaleadd3(-cos_theta,rik_hat,rij_hat,drk);
    scale3(rikinv,drk,drk);
    add3(drj,drk,dri);
    scale3(-1.0,dri);
  }

__device__ void zetaterm_d(rbmd::Real prefactor,rbmd::Real *rij_hat, rbmd::Real rij,
  rbmd::Real rijinv,rbmd::Real *rik_hat, rbmd::Real rik, rbmd::Real rikinv,
  rbmd::Real *dri, rbmd::Real *drj, rbmd::Real *drk,TersoffParams *TersoffParams)
   {
     rbmd::Real gijk,gijk_d,ex_delr,ex_delr_d,fc_v,dfc,cos_theta,tmp;
     rbmd::Real dcosdri[3],dcosdrj[3],dcosdrk[3];

     fc_v = device_fc(rik,TersoffParams);
     dfc = device_fc_d(rik,TersoffParams);
     if (TersoffParams->m_int == 3) tmp = POW(TersoffParams->lambda3,3.0) * (rij-rik);
     else tmp = TersoffParams->lambda3 * (rij-rik);

     if (tmp > 69.0776) ex_delr = 1.e30;
     else if (tmp < -69.0776) ex_delr = 0.0;
     else ex_delr = EXP(tmp);

     if (TersoffParams->m_int == 3)
       ex_delr_d = 3.0*POW((TersoffParams->lambda3),3.0) * POW((rij-rik),2.0)*ex_delr;
     else ex_delr_d = TersoffParams->lambda3 * ex_delr;

     cos_theta = dot3(rij_hat,rik_hat);
     gijk = device_gijk(cos_theta,TersoffParams);
     gijk_d = device_gijk_d(cos_theta,TersoffParams);
     costheta_d(rij_hat,rijinv,rik_hat,rikinv,dcosdri,dcosdrj,dcosdrk);

     // compute the derivative wrt Ri
     // dri = -dfc*gijk*ex_delr*rik_hat;
     // dri += fc*gijk_d*ex_delr*dcosdri;
     // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

     scale3(-dfc*gijk*ex_delr,rik_hat,dri);
     scaleadd3(fc_v*gijk_d*ex_delr,dcosdri,dri,dri);
     scaleadd3(fc_v*gijk*ex_delr_d,rik_hat,dri,dri);
     scaleadd3(-fc_v*gijk*ex_delr_d,rij_hat,dri,dri);
     scale3(prefactor,dri);

     // compute the derivative wrt Rj
     // drj = fc*gijk_d*ex_delr*dcosdrj;
     // drj += fc*gijk*ex_delr_d*rij_hat;

     scale3(fc_v*gijk_d*ex_delr,dcosdrj,drj);
     scaleadd3(fc_v*gijk*ex_delr_d,rij_hat,drj,drj);
     scale3(prefactor,drj);

     // compute the derivative wrt Rk
     // drk = dfc*gijk*ex_delr*rik_hat;
     // drk += fc*gijk_d*ex_delr*dcosdrk;
     // drk += -fc*gijk*ex_delr_d*rik_hat;

     scale3(dfc*gijk*ex_delr,rik_hat,drk);
     scaleadd3(fc_v*gijk_d*ex_delr,dcosdrk,drk,drk);
     scaleadd3(-fc_v*gijk*ex_delr_d,rik_hat,drk,drk);
     scale3(prefactor,drk);
   }

__device__ void attractive(TersoffParams* TersoffParams, rbmd::Real prefactor,
                             rbmd::Real rsqij, rbmd::Real rsqik,
                             rbmd::Real *rij_hat, rbmd::Real *rik_hat,
                             rbmd::Real *fi, rbmd::Real *fj, rbmd::Real *fk)
  {
    rbmd::Real rij,rijinv,rik,rikinv;

    rij = SQRT(rsqij);
    rik = SQRT(rsqik);

    // correct 1/r for shift in rsq
    int shift_flag = 0;
    rbmd::Real shift = 0.5;
    if (shift_flag == 1) {
      rijinv = 1.0/(rij - shift);
      rikinv = 1.0/(rik - shift);
    } else {
      rijinv = 1.0/rij;
      rikinv = 1.0/rik;
    }

    zetaterm_d(prefactor,rij_hat,rij,rijinv,rik_hat,rik,rikinv,fi,fj,fk,TersoffParams);
  }




////////////////////////////
__global__ void ComputeTerSoff(
    Box box,TersoffParams* params, ShiftFlag shift,const rbmd::Real cutmax,
    const rbmd::Id num_atoms,const rbmd::Id nelements,const rbmd::Id* atom_id_to_idx,
    const rbmd::Id* atoms_id, const rbmd::Id* atoms_type,const rbmd::Id* map,
    const rbmd::Id* elem3param,
    const rbmd::Id* start_id, const rbmd::Id* end_id,
    const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial,rbmd::Real* energy)
{
    __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
    temp_storage_elj;
    rbmd::Real sum_elj = 0;

    rbmd::Real forceshiftfac;
  rbmd::Real zeta_ij;
  rbmd::Real fforce, prefactor, fpair,eng;
  rbmd::Real fxtmp,fytmp, fztmp;
  rbmd::Real eng_repul,eng_zeta;
  rbmd::Real fi[3], fj[3], fk[3];
  rbmd::Real r1_hat[3],r2_hat[3];
  for (int i = 0; i < 3; ++i) {
      fi[i] = fj[i]= fk[i] =0.0;
      r1_hat[i] = r2_hat[i]=0.0;
  }
  const rbmd::Real cutshortsq = cutmax*cutmax;

  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_atoms)
  {
    rbmd::Id type1 = atoms_type[tid1];
    rbmd::Id itag = atoms_id[tid1];
    rbmd::Id itype = map[type1+1];

    rbmd::Real x1 = px[tid1];
    rbmd::Real y1 = py[tid1];
    rbmd::Real z1 = pz[tid1];
    fxtmp = fytmp = fztmp = 0.0;
    eng_repul= 0.0;
    eng_zeta= 0.0;
    // two-body interactions,
    for (int j1 = start_id[tid1]; j1 < end_id[tid1]; ++j1)
    {
      rbmd::Id tid2 = id_verletlist[j1];
      //rbmd::Id type2 = atoms_type[tid2];
      rbmd::Id jtag = atoms_id[tid2];
      rbmd::Real x2 = px[tid2];
      rbmd::Real y2 = py[tid2];
      rbmd::Real z2 = pz[tid2];

      //
      if (tid1 > tid2) {
        if ((tid1 + tid2) % 2 == 0) continue;
      }
      else if (tid1 < tid2) {
        if ((tid1 + tid2) % 2 == 1) continue;
      }
      else {
        if (z2 < z1) continue;
        if (z2 == z1 && y2 < y1) continue;
        if (z2 == z1 && y2 == y1 && x2 < x1) continue;
      }


    //
      rbmd::Real x12 = x1 - x2;
      rbmd::Real y12 = y1 - y2;
      rbmd::Real z12 = z1 - z2;
      MinImageDistance(box, x12, y12, z12);
      rbmd::Real r12_2 = x12 * x12 + y12 * y12 + z12 * z12;

      // shift rsq and store correction for force
      if (shift.shift_flag) {
        rbmd::Real  rsqtmp = r12_2 + shift.shift_value*shift.shift_value
        + 2*sqrt(r12_2)*shift.shift_value;
        forceshiftfac = SQRT(rsqtmp/r12_2);
        r12_2 = rsqtmp;
      }

      // //
      // if (tid1 > tid2) {
      //   if ((tid1 + tid2) % 2 == 0) continue;
      // }
      // else if (tid1 < tid2) {
      //   if ((tid1 + tid2) % 2 == 1) continue;
      // }
      // else {
      //   if (z2 < z1) continue;
      //   if (z2 == z1 && y2 < y1) continue;
      //   if (z2 == z1 && y2 == y1 && x2 < x1) continue;
      // }

      rbmd::Id type2 = atoms_type[tid2];
      rbmd::Id jtype = map[type2+1];
      rbmd::Id iparam_ij = elem3param[itype * nelements * nelements + jtype * nelements + jtype];

      if (r12_2 >= params[iparam_ij].cutsq)
        continue;

      repulsive(&params[iparam_ij],r12_2,fpair,eng);
      // correct force for shift in rsq
      if (shift.shift_flag) fpair *= forceshiftfac;
      //printf("fpair eng %f %f\n",fpair,eng);
      fxtmp += x12*fpair;
      fytmp += y12*fpair;
      fztmp += z12*fpair;
      eng_repul += eng;

      atomicAdd(&fx[tid2], -x12 * fpair);
      atomicAdd(&fy[tid2], -y12 * fpair);
      atomicAdd(&fz[tid2], -z12 * fpair);
    }

    //
    rbmd::Real fjxtmp,fjytmp, fjztmp;
    for (int jj = start_id[tid1]; jj < end_id[tid1]; ++jj)
    {
      rbmd::Id tid_j = id_verletlist[jj];
      rbmd::Id jtype = map[atoms_type[tid_j]+1];
      rbmd::Id jtag = atoms_id[tid_j];
      rbmd::Real x21 = px[tid_j] - x1;
      rbmd::Real y21 = py[tid_j] - y1;
      rbmd::Real z21 = pz[tid_j] - z1;

      MinImageDistance(box, x21, y21, z21);
      rbmd::Real rsq1 = x21*x21 + y21*y21 + z21*z21;
      if (shift.shift_flag)
        rsq1 += shift.shift_value*shift.shift_value + 2*SQRT(rsq1)*shift.shift_value;
      rbmd::Id iparam_ij = elem3param[itype * nelements * nelements + jtype * nelements + jtype];

      if (rsq1 >= params[iparam_ij].cutsq)
        continue;

      const rbmd::Real r1inv = 1 / SQRT(rsq1);
      r1_hat[0] = x21 * r1inv;
      r1_hat[1] = y21 * r1inv;
      r1_hat[2] = z21 * r1inv;

      fjxtmp =fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;

      // Loop over k neighbors
      for (int kk = start_id[tid1]; kk < end_id[tid1]; ++kk)
      {
        if (kk == jj)
          continue;

        rbmd::Id tid_k = id_verletlist[kk];
        rbmd::Id ktag = atoms_id[tid_k];

        rbmd::Id ktype = map[atoms_type[tid_k]+1];
        rbmd::Id iparam_ijk = elem3param[itype * nelements * nelements + jtype * nelements + ktype];

        rbmd::Real x31 = px[tid_k] - x1;
        rbmd::Real y31 = py[tid_k] - y1;
        rbmd::Real z31  = pz[tid_k] - z1;

        MinImageDistance(box, x31, y31, z31);
        rbmd::Real rsq2 = x31*x31 + y31*y31 + z31*z31;
        if (shift.shift_flag)
          rsq2 += shift.shift_value*shift.shift_value + 2*SQRT(rsq2)*shift.shift_value;

        if (rsq2 >= params[iparam_ij].cutsq)
          continue;

        rbmd::Real r2inv = 1/ SQRT(rsq2);
        r2_hat[0] = x31 * r2inv;
        r2_hat[1] = y31 * r2inv;
        r2_hat[2] = z31 * r2inv;

        zeta_ij += zeta(&params[iparam_ijk], rsq1, rsq2, r1_hat, r2_hat);
      }

      // zeta
      force_zeta(&params[iparam_ij], rsq1, zeta_ij, fforce, prefactor, eng);
      fpair = fforce * r1inv;

      fxtmp += x21 * fpair;
      fytmp += y21 * fpair;
      fztmp += z21 * fpair;
      fjxtmp -= x21 * fpair;
      fjytmp -= y21 * fpair;
      fjztmp -= z21 * fpair;
      eng_zeta += eng;

      // attractive
      for (int kk = start_id[tid1]; kk < end_id[tid1]; ++kk)
      {
        if (kk == jj)
          continue;

        rbmd::Id tid_k = id_verletlist[kk];
        rbmd::Id ktag = atoms_id[tid_k];

        rbmd::Id ktype = map[atoms_type[tid_k]+1];
        rbmd::Id iparam_ijk = elem3param[itype * nelements * nelements + jtype * nelements + ktype];

        rbmd::Real x31 = px[tid_k] - x1;
        rbmd::Real y31 = py[tid_k] - y1;
        rbmd::Real z31  = pz[tid_k] - z1;
        MinImageDistance(box, x31, y31, z31);
        rbmd::Real rsq2 = x31*x31 + y31*y31 + z31*z31;
        if (shift.shift_flag)
          rsq2 += shift.shift_value*shift.shift_value + 2*SQRT(rsq2)*shift.shift_value;

        if (rsq2 >= params[iparam_ij].cutsq)
          continue;

        rbmd::Real r2inv = 1 / SQRT(rsq2);
        r2_hat[0] = x31 * r2inv;
        r2_hat[1] = y31 * r2inv;
        r2_hat[2] = z31 * r2inv;

        attractive(&params[iparam_ijk], prefactor, rsq1, rsq2, r1_hat, r2_hat, fi, fj, fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];

        atomicAdd(&fx[tid_k], fk[0]);
        atomicAdd(&fy[tid_k], fk[1]);
        atomicAdd(&fz[tid_k], fk[2]);
      }
      atomicAdd(&fx[tid_j], fjxtmp);
      atomicAdd(&fy[tid_j], fjytmp);
      atomicAdd(&fz[tid_j], fjztmp);
    }
    //
    atomicAdd(&fx[tid1], fxtmp);
    atomicAdd(&fy[tid1], fytmp);
    atomicAdd(&fz[tid1], fztmp);
    sum_elj = eng_repul + eng_zeta;

  }
    rbmd::Real block_sum_elj =
    BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage_elj)
        .Sum(sum_elj);

    if (threadIdx.x == 0) {
      atomicAdd(energy, block_sum_elj);
    }
}

__global__ void ComputeTwoBodyTerSoff(
    Box box, TersoffParams* params, const rbmd::Real cutmax,
    const rbmd::Id num_atoms, const rbmd::Id nelements, const rbmd::Id* atoms_id,
    const rbmd::Id* atoms_type, const rbmd::Id* map, const rbmd::Id* elem3param,
    const rbmd::Id* start_id, const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
    const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz)
{
    rbmd::Real fxtmp,fytmp, fztmp;

  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_atoms) {
    rbmd::Id itag = atoms_id[tid1];
    rbmd::Id type1 = atoms_type[tid1];
    rbmd::Id itype = map[type1+1];
    rbmd::Real x1 = px[tid1];
    rbmd::Real y1 = py[tid1];
    rbmd::Real z1 = pz[tid1];

     fxtmp = fytmp = fztmp = 0.0;

    for (int j1 = start_id[tid1]; j1 < end_id[tid1]; ++j1) {
      rbmd::Id tid2 = id_verletlist[j1];
      rbmd::Id jtag = atoms_id[tid2];
      rbmd::Id type2 = atoms_type[tid2];
      rbmd::Real x2 = px[tid2];
      rbmd::Real y2 = py[tid2];
      rbmd::Real z2 = pz[tid2];

      // Distance calculation with PBC
      rbmd::Real x12 = x1 - x2;
      rbmd::Real y12 = y1 - y2;
      rbmd::Real z12 = z1 - z2;

      MinImageDistance(box, x12, y12, z12);
      rbmd::Real r12_2 = x12 * x12 + y12 * y12 + z12 * z12;

      // Unique pair check
      if (itag > jtag) {
        if ((itag + jtag) % 2 == 0) continue;
      }
      else if (itag < jtag) {
        if ((itag + jtag) % 2 == 1) continue;
      }
      else {
        if (z2 < z1) continue;
        if (z2 == z1 && y2 < y1) continue;
        if (z2 == z1 && y2 == y1 && x2 < x1) continue;
      }

      rbmd::Id jtype = map[type2+1];
      rbmd::Id iparam_ij = elem3param[itype * nelements * nelements + jtype * nelements + jtype];

      if (r12_2 >= params[iparam_ij].cutsq) continue;

      // Calculate repulsive force
      rbmd::Real fpair, eng;
      repulsive(&params[iparam_ij], r12_2, fpair, eng);

      // Update forces with atomic operations
      rbmd::Real fx_ij = x12 * fpair;
      rbmd::Real fy_ij = y12 * fpair;
      rbmd::Real fz_ij = z12 * fpair;

      fxtmp += fx_ij;
      fytmp += fy_ij;
      fztmp += fz_ij;

      // atomicAdd(&fx[tid1], fx_ij);
      // atomicAdd(&fy[tid1], fy_ij);
      // atomicAdd(&fz[tid1], fz_ij);

      // atomicAdd(&fx[tid2], -fx_ij);
      // atomicAdd(&fy[tid2], -fy_ij);
      // atomicAdd(&fz[tid2], -fz_ij);
    }
    fx[tid1]=  fxtmp;
    fy[tid1] = fytmp;
    fz[tid1] = fztmp;
    // atomicAdd(&fx[tid1], fxtmp);
    // atomicAdd(&fy[tid1], fytmp);
    // atomicAdd(&fz[tid1], fztmp);
  }

}

__global__ void ComputeThreeBodyTerSoff(
    Box box, TersoffParams* params, const rbmd::Real cutmax,
    const rbmd::Id num_atoms, const rbmd::Id nelements, const rbmd::Id* atoms_id,
    const rbmd::Id* atoms_type, const rbmd::Id* map, const rbmd::Id* elem3param,
    const rbmd::Id* start_id, const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
    const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz)
{
    rbmd::Real fxtmp_z,fytmp_z, fztmp_z;
    rbmd::Real fxtmp_a,fytmp_a, fztmp_a;
  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_atoms) {
    fxtmp_z = fytmp_z = fztmp_z = 0.0;
    fxtmp_a = fytmp_a = fztmp_a = 0.0;

    rbmd::Id itag = atoms_id[tid1];
    rbmd::Id type1 = atoms_type[tid1];
    rbmd::Id itype = map[type1+1];
    rbmd::Real x1 = px[tid1];
    rbmd::Real y1 = py[tid1];
    rbmd::Real z1 = pz[tid1];

    rbmd::Real fi[3] = {0.0, 0.0, 0.0};
    rbmd::Real r1_hat[3], r2_hat[3];

    rbmd::Real fjxtmp_z,fjytmp_z, fjztmp_z;
    rbmd::Real fjxtmp_a,fjytmp_a, fjztmp_a;
    for (int jj = start_id[tid1]; jj < end_id[tid1]; ++jj) {
      fjxtmp_z = fjytmp_z = fjztmp_z = 0.0;
      fjxtmp_a = fjytmp_a = fjztmp_a = 0.0;

      rbmd::Id tid_j = id_verletlist[jj];
      rbmd::Id jtype = map[atoms_type[tid_j] + 1];
      rbmd::Real x21 = px[tid_j] - x1;
      rbmd::Real y21 = py[tid_j] - y1;
      rbmd::Real z21 = pz[tid_j] - z1;

      MinImageDistance(box, x21, y21, z21);
      rbmd::Real rsq1 = x21*x21 + y21*y21 + z21*z21;

      rbmd::Id iparam_ij = elem3param[itype * nelements * nelements + jtype * nelements + jtype];
      if (rsq1 >= params[iparam_ij].cutsq) continue;

      const rbmd::Real r1inv = 1 / SQRT(rsq1);
      r1_hat[0] = x21 * r1inv;
      r1_hat[1] = y21 * r1inv;
      r1_hat[2] = z21 * r1inv;

      rbmd::Real zeta_ij = 0.0;
      rbmd::Real fj_force[3] = {0.0, 0.0, 0.0};

      // zeta
      for (int kk = start_id[tid1]; kk < end_id[tid1]; ++kk) {
        if (kk == jj) continue;

        rbmd::Id tid_k = id_verletlist[kk];
        rbmd::Id ktype = map[atoms_type[tid_k] + 1];

        rbmd::Real x31 = px[tid_k] - x1;
        rbmd::Real y31 = py[tid_k] - y1;
        rbmd::Real z31  = pz[tid_k] - z1;

        MinImageDistance(box, x31, y31, z31);
        rbmd::Real rsq2 = x31*x31 + y31*y31 + z31*z31;

        rbmd::Id iparam_ijk = elem3param[itype * nelements * nelements + jtype * nelements + ktype];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;

        rbmd::Real r2inv = 1 / SQRT(rsq2);
        r2_hat[0] = x31 * r2inv;
        r2_hat[1] = y31 * r2inv;
        r2_hat[2] = z31 * r2inv;

        zeta_ij += zeta(&params[iparam_ijk], rsq1, rsq2, r1_hat, r2_hat);
      }

      //
      rbmd::Real fforce, prefactor, eng;
      force_zeta(&params[iparam_ij], rsq1, zeta_ij, fforce, prefactor, eng);

      rbmd::Real fpair = fforce * r1inv;
      rbmd::Real fx_ij = x21 * fpair;
      rbmd::Real fy_ij = y21 * fpair;
      rbmd::Real fz_ij = z21 * fpair;

      // fi[0] += fx_ij;
      // fi[1] += fy_ij;
      // fi[2] += fz_ij;
      fxtmp_z += fx_ij;
      fytmp_z += fy_ij;
      fztmp_z += fz_ij;
      fjxtmp_z -= fx_ij;
      fjytmp_z -= fy_ij;
      fjztmp_z -= fz_ij;


      // atomicAdd(&fx[tid_j], -fx_ij);
      // atomicAdd(&fy[tid_j], -fy_ij);
      // atomicAdd(&fz[tid_j], -fz_ij);

      //
      for (int kk = start_id[tid1]; kk < end_id[tid1]; ++kk) {
        if (kk == jj) continue;

        rbmd::Id tid_k = id_verletlist[kk];
        rbmd::Id ktype = map[atoms_type[tid_k] + 1];

        rbmd::Real x31 = px[tid_k] - x1;
        rbmd::Real y31 = py[tid_k] - y1;
        rbmd::Real z31  = pz[tid_k] - z1;

        MinImageDistance(box, x31, y31, z31);
        rbmd::Real rsq2 = x31*x31 + y31*y31 + z31*z31;

        rbmd::Id iparam_ijk = elem3param[itype * nelements * nelements + jtype * nelements + ktype];
        TersoffParams param_ijk = params[iparam_ijk];

        if (rsq2 >= param_ijk.cutsq) continue;

        rbmd::Real r2inv = 1 / SQRT(rsq2);
        r2_hat[0] = x31 * r2inv;
        r2_hat[1] = y31 * r2inv;
        r2_hat[2] = z31 * r2inv;

        rbmd::Real fi_tmp[3], fj_tmp[3], fk_tmp[3];
        attractive(&param_ijk, prefactor, rsq1, rsq2, r1_hat, r2_hat,
                   fi_tmp, fj_tmp, fk_tmp);

        // //
        // fi[0] += fi_tmp[0];
        // fi[1] += fi_tmp[1];
        // fi[2] += fi_tmp[2];
        //
        // //
        // atomicAdd(&fx[tid_j], fj_tmp[0]);
        // atomicAdd(&fy[tid_j], fj_tmp[1]);
        // atomicAdd(&fz[tid_j], fj_tmp[2]);

        fxtmp_a += fi_tmp[0];
        fytmp_a += fi_tmp[1];
        fztmp_a += fi_tmp[2];
        fjxtmp_a += fj_tmp[0];
        fjytmp_a += fj_tmp[1];
        fjztmp_a += fj_tmp[2];

        // atomicAdd(&fx[tid_j], fj_tmp[0]);
        // atomicAdd(&fy[tid_j], fj_tmp[1]);
        // atomicAdd(&fz[tid_j], fj_tmp[2]);

        //
        atomicAdd(&fx[tid_k], fk_tmp[0]);
        atomicAdd(&fy[tid_k], fk_tmp[1]);
        atomicAdd(&fz[tid_k], fk_tmp[2]);
      }
      auto fjxtmp = fjxtmp_z + fjxtmp_a;
      auto fjytmp = fjytmp_z + fjytmp_a;
      auto fjztmp = fjztmp_z + fjztmp_a;

      atomicAdd(&fx[tid_j], fjxtmp);
      atomicAdd(&fy[tid_j], fjytmp);
      atomicAdd(&fz[tid_j], fjztmp);
    }

    //
    auto fxtmp = fxtmp_z + fxtmp_a;
    auto fytmp = fytmp_z + fytmp_a;
    auto fztmp = fztmp_z + fztmp_a;
    atomicAdd(&fx[tid1], fxtmp);
    atomicAdd(&fy[tid1], fytmp);
    atomicAdd(&fz[tid1], fztmp);
  }
}


void TerSoff<device::DEVICE_GPU>::operator()(
    Box box, TersoffParams* params, ShiftFlag shift,const rbmd::Real cutmax,
    const rbmd::Id num_atoms, const rbmd::Id nelements, const rbmd::Id* atom_id_to_idx,
    const rbmd::Id* atoms_id, const rbmd::Id* atoms_type, const rbmd::Id* map,
    const rbmd::Id* elem3param,
    const rbmd::Id* start_id, const rbmd::Id* end_id,
    const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial, rbmd::Real* energy)
{
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_KERNEL(ComputeTerSoff<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        box,params,shift,cutmax,num_atoms,nelements, atom_id_to_idx,atoms_id, atoms_type,
        map,elem3param,
        start_id,end_id, id_verletlist,
        px, py, pz,fx, fy, fz,
        flat_virial,energy));

  // //
  // CHECK_KERNEL(ComputeTwoBodyTerSoff<<<blocks_per_grid, BLOCK_SIZE>>>(
  //     box, params, cutmax, num_atoms, nelements, atoms_id, atoms_type,
  //     map, elem3param, start_id, end_id, id_verletlist,
  //     px, py, pz, fx, fy, fz));
  //
  // cudaDeviceSynchronize();  //
  //
  // //
  // CHECK_KERNEL(ComputeThreeBodyTerSoff<<<blocks_per_grid, BLOCK_SIZE>>>(
  //     box, params, cutmax, num_atoms, nelements, atoms_id, atoms_type,
  //     map, elem3param, start_id, end_id, id_verletlist,
  //     px, py, pz, fx, fy, fz));

 }

void ThreeBodyTerSoff<device::DEVICE_GPU>::operator()(
    Box box, TersoffParams* params, const rbmd::Real cutmax,
    const rbmd::Id num_atoms, const rbmd::Id nelements, const rbmd::Id* atom_id_to_idx,
    const rbmd::Id* atoms_id, const rbmd::Id* atoms_type, const rbmd::Id* map,
    const rbmd::Id* elem3param,
    const rbmd::Id* start_id, const rbmd::Id* end_id,
    const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial, rbmd::Real* energy)
  {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

    //
    CHECK_KERNEL(ComputeThreeBodyTerSoff<<<blocks_per_grid, BLOCK_SIZE>>>(
        box, params, cutmax, num_atoms, nelements, atoms_id, atoms_type,
        map, elem3param, start_id, end_id, id_verletlist,
        px, py, pz, fx, fy, fz));

  }

}


