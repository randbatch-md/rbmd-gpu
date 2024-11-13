// #include <hip/hip_runtime.h>
//
// #include "../common/rbmd_define.h"
// #include "tersoff_op.h"
// #include "tersoff.h"
// #include "model/box.h"
//
// rbmd::Id maxshort = 10;
//
// namespace op{
// //---------Math---------//
//   void scale3(const rbmd::Real s, rbmd::Real *v)
//   {
//     v[0] *= s;
//     v[1] *= s;
//     v[2] *= s;
//   }
//
//   void scale3(const rbmd::Real s, const rbmd::Real *v, rbmd::Real *ans)
//   {
//     ans[0] = s * v[0];
//     ans[1] = s * v[1];
//     ans[2] = s * v[2];
//   }
//
//   rbmd::Real dot3(const rbmd::Real *v1, const rbmd::Real *v2)
//   {
//     return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
//   }
//
//   void scaleadd3(const rbmd::Real s, const rbmd::Real *v1, const rbmd::Real *v2, rbmd::Real *ans)
//   {
//     ans[0] = s * v1[0] + v2[0];
//     ans[1] = s * v1[1] + v2[1];
//     ans[2] = s * v1[2] + v2[2];
//   }
//
//   void  add3(const rbmd::Real *v1, const rbmd::Real *v2, rbmd::Real *ans)
//   {
//     ans[0] = v1[0] + v2[0];
//     ans[1] = v1[1] + v2[1];
//     ans[2] = v1[2] + v2[2];
//   }
//
// //---------device---------//
//   // inlined functions for efficiency
//   inline rbmd::Real gijk(const rbmd::Real costheta, const TersoffParams* TersoffParams)
//     {
//       const rbmd::Real ters_c = TersoffParams->c * TersoffParams->c;
//       const rbmd::Real ters_d = TersoffParams->d * TersoffParams->d;
//       const rbmd::Real hcth = TersoffParams->h - costheta;
//
//       return TersoffParams->gamma * (1.0 + ters_c / ters_d - ters_c / (ters_d + hcth * hcth));
//     }
//
//   inline rbmd::Real gijk_d(const rbmd::Real costheta, const TersoffParams*TersoffParams)
//     {
//       const rbmd::Real ters_c = TersoffParams->c * TersoffParams->c;
//       const rbmd::Real ters_d = TersoffParams->d * TersoffParams->d;
//       const rbmd::Real hcth = TersoffParams->h - costheta;
//       const rbmd::Real numerator = -2.0 * ters_c * hcth;
//       const rbmd::Real denominator = 1.0 / (ters_d + hcth * hcth);
//       return TersoffParams->gamma * numerator * denominator * denominator;
//     }
//
//   //fc
//   rbmd::Real fc(rbmd::Real r, TersoffParams *TersoffParams)
//   {
//     rbmd::Real ters_R = TersoffParams->bigr;
//     rbmd::Real ters_D = TersoffParams->bigd;
//
//     if (r < ters_R-ters_D) return 1.0;
//     if (r > ters_R+ters_D) return 0.0;
//     return 0.5*(1.0 - SIN(M_PI_2*(r - ters_R)/ters_D));
//   }
//
//   rbmd::Real fc_d(rbmd::Real r, TersoffParams* TersoffParams)
//     {
//       rbmd::Real ters_R = TersoffParams->bigr;
//       rbmd::Real ters_D = TersoffParams->bigd;
//
//       if (r < ters_R-ters_D) return 0.0;
//       if (r > ters_R+ters_D) return 0.0;
//       return -(M_PI_4/ters_D) * COS(M_PI_2*(r - ters_R)/ters_D);
//     }
//
//    //fa
//   rbmd::Real fa(rbmd::Real r, TersoffParams* TersoffParams)
//     {
//       if (r > TersoffParams->bigr + TersoffParams->bigd) return 0.0;
//       return -TersoffParams->bigb * exp(-TersoffParams->lam2 * r) * fc(r,TersoffParams);
//     }
//
//   rbmd::Real fa_d(rbmd::Real r, TersoffParams* TersoffParams)
//     {
//       if (r > TersoffParams->bigr + TersoffParams->bigd) return 0.0;
//       return TersoffParams->bigb * exp(-TersoffParams->lam2 * r) *
//         (TersoffParams->lam2 * fc(r,TersoffParams) - fc_d(r,TersoffParams));
//     }
//
//   // zeta
//   rbmd::Real zeta(TersoffParams* TersoffParams, rbmd::Real rsqij, rbmd::Real rsqik,
//                            rbmd::Real *rij_hat, rbmd::Real *rik_hat)
//     {
//       rbmd::Real rij,rik,costheta,arg,ex_delr;
//
//       rij = SQRT(rsqij);
//       rik = SQRT(rsqik);
//       costheta = dot3(rij_hat,rik_hat);
//
//       if (TersoffParams->powermint == 3) arg = POW(TersoffParams->lam3 * (rij-rik),3.0);
//       else arg = TersoffParams->lam3 * (rij-rik);
//
//       if (arg > 69.0776) ex_delr = 1.e30;
//       else if (arg < -69.0776) ex_delr = 0.0;
//       else ex_delr = EXP(arg);
//
//       return fc(rik,TersoffParams) * gijk(costheta,TersoffParams) * ex_delr;
//     }
//
//   //bij
//   rbmd::Real bij(rbmd::Real zeta, TersoffParams* TersoffParams)
//     {
//       rbmd::Real tmp = TersoffParams->beta * zeta;
//       if (tmp > TersoffParams->c1) return 1.0/sqrt(tmp);
//       if (tmp > TersoffParams->c2)
//         return (1.0 - pow(tmp,-TersoffParams->powern) / (2.0*TersoffParams->powern))/sqrt(tmp);
//       if (tmp < TersoffParams->c4) return 1.0;
//       if (tmp < TersoffParams->c3)
//         return 1.0 - pow(tmp,TersoffParams->powern)/(2.0*TersoffParams->powern);
//       return pow(1.0 + pow(tmp,TersoffParams->powern), -1.0/(2.0*TersoffParams->powern));
//     }
//
//   rbmd::Real bij_d(rbmd::Real zeta, TersoffParams* TersoffParams)
//     {
//       rbmd::Real tmp = TersoffParams->beta * zeta;
//       if (tmp > TersoffParams->c1) return TersoffParams->beta * -0.5*pow(tmp,-1.5);
//       if (tmp > TersoffParams->c2)
//         return TersoffParams->beta * (-0.5*pow(tmp,-1.5) *
//                               // error in negligible 2nd term fixed 9/30/2015
//                               // (1.0 - 0.5*(1.0 +  1.0/(2.0*TersoffParams->powern)) *
//                               (1.0 - (1.0 +  1.0/(2.0*TersoffParams->powern)) *
//                                pow(tmp,-TersoffParams->powern)));
//       if (tmp < TersoffParams->c4) return 0.0;
//       if (tmp < TersoffParams->c3)
//         return -0.5*TersoffParams->beta * pow(tmp,TersoffParams->powern-1.0);
//
//       rbmd::Real tmp_n = pow(tmp,TersoffParams->powern);
//       return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*TersoffParams->powern)))*tmp_n / zeta;
//     }
//
//
//
//
//
// void repulsive(TersoffParams* TersoffParams,rbmd::Real rsq, rbmd::Real &fforce,
//   rbmd::Id eflag,rbmd::Real &eng)
//   {
//     rbmd::Real r,tmp_fc,tmp_fc_d,tmp_exp;
//
//     r = SQRT(rsq);
//     tmp_fc = fc(r,TersoffParams);
//     tmp_fc_d = fc_d(r,TersoffParams);
//     tmp_exp = EXP(-TersoffParams->lam1 * r);
//     fforce = -TersoffParams->biga * tmp_exp * (tmp_fc_d - tmp_fc*TersoffParams->lam1) / r;
//     if (eflag) eng = tmp_fc * TersoffParams->biga * tmp_exp;
//   }
//
//
//
// void costheta_d(rbmd::Real *rij_hat, rbmd::Real rijinv,
//                              rbmd::Real *rik_hat, rbmd::Real rikinv,
//                              rbmd::Real *dri, rbmd::Real *drj, rbmd::Real *drk)
//   {
//     // first element is devative wrt Ri, second wrt Rj, third wrt Rk
//
//     rbmd::Real cos_theta = dot3(rij_hat,rik_hat);
//
//     scaleadd3(-cos_theta,rij_hat,rik_hat,drj);
//     scale3(rijinv,drj);
//     scaleadd3(-cos_theta,rik_hat,rij_hat,drk);
//     scale3(rikinv,drk,drk);
//     add3(drj,drk,dri);
//     scale3(-1.0,dri);
//   }
//
// void zetaterm_d(rbmd::Real prefactor,rbmd::Real *rij_hat, rbmd::Real rij,
//   rbmd::Real rijinv,rbmd::Real *rik_hat, rbmd::Real rik, rbmd::Real rikinv,
//   rbmd::Real *dri, rbmd::Real *drj, rbmd::Real *drk,TersoffParams *TersoffParams)
//   {
//     rbmd::Real gijk,gijk_d,ex_delr,ex_delr_d,fc_v,dfc,cos_theta,tmp;
//     rbmd::Real dcosdri[3],dcosdrj[3],dcosdrk[3];
//
//     fc_v = fc(rik,TersoffParams);
//     dfc = fc_d(rik,TersoffParams);
//     if (TersoffParams->powermint == 3) tmp = POW((TersoffParams->lam3,3.0) * (rij-rik));
//     else tmp = TersoffParams->lam3 * (rij-rik);
//
//     if (tmp > 69.0776) ex_delr = 1.e30;
//     else if (tmp < -69.0776) ex_delr = 0.0;
//     else ex_delr = EXP(tmp);
//
//     if (TersoffParams->powermint == 3)
//       ex_delr_d = 3.0*POW((TersoffParams->lam3),3.0) * POW((rij-rik),2.0)*ex_delr;
//     else ex_delr_d = TersoffParams->lam3 * ex_delr;
//
//     cos_theta = dot3(rij_hat,rik_hat);
//     gijk = gijk(cos_theta,TersoffParams);
//     gijk_d = gijk_d(cos_theta,TersoffParams);
//     costheta_d(rij_hat,rijinv,rik_hat,rikinv,dcosdri,dcosdrj,dcosdrk);
//
//     // compute the derivative wrt Ri
//     // dri = -dfc*gijk*ex_delr*rik_hat;
//     // dri += fc*gijk_d*ex_delr*dcosdri;
//     // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);
//
//     scale3(-dfc*gijk*ex_delr,rik_hat,dri);
//     scaleadd3(fc_v*gijk_d*ex_delr,dcosdri,dri,dri);
//     scaleadd3(fc_v*gijk*ex_delr_d,rik_hat,dri,dri);
//     scaleadd3(-fc_v*gijk*ex_delr_d,rij_hat,dri,dri);
//     scale3(prefactor,dri);
//
//     // compute the derivative wrt Rj
//     // drj = fc*gijk_d*ex_delr*dcosdrj;
//     // drj += fc*gijk*ex_delr_d*rij_hat;
//
//     scale3(fc_v*gijk_d*ex_delr,dcosdrj,drj);
//     scaleadd3(fc_v*gijk*ex_delr_d,rij_hat,drj,drj);
//     scale3(prefactor,drj);
//
//     // compute the derivative wrt Rk
//     // drk = dfc*gijk*ex_delr*rik_hat;
//     // drk += fc*gijk_d*ex_delr*dcosdrk;
//     // drk += -fc*gijk*ex_delr_d*rik_hat;
//
//     scale3(dfc*gijk*ex_delr,rik_hat,drk);
//     scaleadd3(fc_v*gijk_d*ex_delr,dcosdrk,drk,drk);
//     scaleadd3(-fc_v*gijk*ex_delr_d,rik_hat,drk,drk);
//     scale3(prefactor,drk);
//   }
//
//
//
//
// __global__ void ComputeTerSoff(
//     Box* box, const rbmd::Real cut_off,
//     const rbmd::Id num_atoms,const rbmd::Id nelements,
//     const rbmd::Id* atoms_type,const rbmd::Id* map,
//     const rbmd::Id* tag,const rbmd::Id* elem3param ,
//     const rbmd::Id* start_id, const rbmd::Id* end_id,
//     const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
//     const rbmd::Real* pz, rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
//     rbmd::Real*neighshort)
// {
//   const rbmd::Real cut_2 = cut_off*cut_off;
//   unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
//   if (tid1 < num_atoms)
//   {
//     rbmd::Id type1 = atoms_type[tid1];
//     rbmd::Id itag = tag[type1];
//     rbmd::Id itype = map[type1];
//
//     rbmd::Real x1 = px[tid1];
//     rbmd::Real y1 = py[tid1];
//     rbmd::Real z1 = pz[tid1];
//     rbmd::Id numshort = 0;
//
//
//     // two-body interactions,
//     for (int j1 = start_id[tid1]; j1 < end_id[tid1]; ++j1)
//     {
//       rbmd::Id tid2 = id_verletlist[j1];
//       rbmd::Id type2 = atoms_type[tid2];
//       rbmd::Real x2 = px[tid2];
//       rbmd::Real y2 = py[tid2];
//       rbmd::Real z2 = pz[tid2];
//
//       rbmd::Real x12 = x2 - x1;
//       rbmd::Real y12 = y2 - y1;
//       rbmd::Real z12 = z2 - z1;
//
//       MinImageDistance(box, x12, y12, z12);
//       rbmd::Real r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
//       if (r12_2 < cut_2)
//       {
//         neighshort[numshort++] = tid2;
//         if (numshort >= maxshort) {
//           maxshort += maxshort/2;
//         }
//       }
//
//       // rbmd::Id jtag = tag[tid2];
//       // if (itag > jtag){
//       //   if ((itag+jtag) % 2 == 0) continue;
//       // } else if (itag < jtag) {
//       //   if ((itag+jtag) % 2 == 1) continue;
//       // }else{
//       //   if (z2  < z1 ) continue;
//       //   if (z2 == z1 && y2 < y1 ) continue;
//       //   if (z2  == z1 && y2 == y1  && x2  < x1) continue;
//       // }
//
//       rbmd::Id jtype = map[type2];
//       rbmd::Id iparam_ij = elem3param[itype * nelements * nelements
//                                        + jtype * nelements + jtype];
//       if (r12_2 >= params[iparam_ij].cutsq) continue;
//
//       rbmd::Real fpair,evdwl;
//       repulsive(&params[iparam_ij],r12_2,fpair,evdwl);
//
//     }
//     // three-body interactions
//
//
//   }
// }
//
// }
//
//
