#pragma once

#include <cmath>
#include <limits>
#include <algorithm> // For std::swap and std::fabs
#include "../common/rbmd_define.h"

namespace MathLib {

static constexpr rbmd::Id nmaxfactorial = 167;

static const double nfac_table[] = {
  1,
  1,
  2,
  6,
  24,
  120,
  720,
  5040,
  40320,
  362880,
  3628800,
  39916800,
  479001600,
  6227020800,
  87178291200,
  1307674368000,
  20922789888000,
  355687428096000,
  6.402373705728e+15,
  1.21645100408832e+17,
  2.43290200817664e+18,
  5.10909421717094e+19,
  1.12400072777761e+21,
  2.5852016738885e+22,
  6.20448401733239e+23,
  1.5511210043331e+25,
  4.03291461126606e+26,
  1.08888694504184e+28,
  3.04888344611714e+29,
  8.8417619937397e+30,
  2.65252859812191e+32,
  8.22283865417792e+33,
  2.63130836933694e+35,
  8.68331761881189e+36,
  2.95232799039604e+38,
  1.03331479663861e+40,
  3.71993326789901e+41,
  1.37637530912263e+43,
  5.23022617466601e+44,
  2.03978820811974e+46,
  8.15915283247898e+47,
  3.34525266131638e+49,
  1.40500611775288e+51,
  6.04152630633738e+52,
  2.65827157478845e+54,
  1.1962222086548e+56,
  5.50262215981209e+57,
  2.58623241511168e+59,
  1.24139155925361e+61,
  6.08281864034268e+62,
  3.04140932017134e+64,
  1.55111875328738e+66,
  8.06581751709439e+67,
  4.27488328406003e+69,
  2.30843697339241e+71,
  1.26964033536583e+73,
  7.10998587804863e+74,
  4.05269195048772e+76,
  2.35056133128288e+78,
  1.3868311854569e+80,
  8.32098711274139e+81,
  5.07580213877225e+83,
  3.14699732603879e+85,
  1.98260831540444e+87,
  1.26886932185884e+89,
  8.24765059208247e+90,
  5.44344939077443e+92,
  3.64711109181887e+94,
  2.48003554243683e+96,
  1.71122452428141e+98,
  1.19785716699699e+100,
  8.50478588567862e+101,
  6.12344583768861e+103,
  4.47011546151268e+105,
  3.30788544151939e+107,
  2.48091408113954e+109,
  1.88549470166605e+111,
  1.45183092028286e+113,
  1.13242811782063e+115,
  8.94618213078297e+116,
  7.15694570462638e+118,
  5.79712602074737e+120,
  4.75364333701284e+122,
  3.94552396972066e+124,
  3.31424013456535e+126,
  2.81710411438055e+128,
  2.42270953836727e+130,
  2.10775729837953e+132,
  1.85482642257398e+134,
  1.65079551609085e+136,
  1.48571596448176e+138,
  1.3520015276784e+140,
  1.24384140546413e+142,
  1.15677250708164e+144,
  1.08736615665674e+146,
  1.03299784882391e+148,
  9.91677934870949e+149,
  9.61927596824821e+151,
  9.42689044888324e+153,
  9.33262154439441e+155,
  9.33262154439441e+157,
  9.42594775983835e+159,
  9.61446671503512e+161,
  9.90290071648618e+163,
  1.02990167451456e+166,
  1.08139675824029e+168,
  1.14628056373471e+170,
  1.22652020319614e+172,
  1.32464181945183e+174,
  1.44385958320249e+176,
  1.58824554152274e+178,
  1.76295255109024e+180,
  1.97450685722107e+182,
  2.23119274865981e+184,
  2.54355973347219e+186,
  2.92509369349301e+188,
  3.3931086844519e+190,
  3.96993716080872e+192,
  4.68452584975429e+194,
  5.5745857612076e+196,
  6.68950291344912e+198,
  8.09429852527344e+200,
  9.8750442008336e+202,
  1.21463043670253e+205,
  1.50614174151114e+207,
  1.88267717688893e+209,
  2.37217324288005e+211,
  3.01266001845766e+213,
  3.8562048236258e+215,
  4.97450422247729e+217,
  6.46685548922047e+219,
  8.47158069087882e+221,
  1.118248651196e+224,
  1.48727070609069e+226,
  1.99294274616152e+228,
  2.69047270731805e+230,
  3.65904288195255e+232,
  5.01288874827499e+234,
  6.91778647261949e+236,
  9.61572319694109e+238,
  1.34620124757175e+241,
  1.89814375907617e+243,
  2.69536413788816e+245,
  3.85437071718007e+247,
  5.5502938327393e+249,
  8.04792605747199e+251,
  1.17499720439091e+254,
  1.72724589045464e+256,
  2.55632391787286e+258,
  3.80892263763057e+260,
  5.71338395644585e+262,
  8.62720977423323e+264,
  1.31133588568345e+267,
  2.00634390509568e+269,
  3.08976961384735e+271,
  4.78914290146339e+273,
  7.47106292628289e+275,
  1.17295687942641e+278,
  1.85327186949373e+280,
  2.94670227249504e+282,
  4.71472363599206e+284,
  7.59070505394721e+286,
  1.22969421873945e+289,
  2.0044015765453e+291,
  3.28721858553429e+293,
  5.42391066613159e+295,
  9.00369170577843e+297,
  1.503616514865e+300, // nmaxfactorial = 167
};

/**
 * @brief Calculate the dot product of two 3D vectors.
 */
inline rbmd::Real dot3(const rbmd::Real* a, const rbmd::Real* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * @brief Calculate the cross product of two 3D vectors c = a x b.
 */
inline void cross3(const rbmd::Real* a, const rbmd::Real* b, rbmd::Real* c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/**
 * @brief Reverse the 3D vector: v = -v
 */
inline void negate3(rbmd::Real* v) {
    v[0] = -v[0];
    v[1] = -v[1];
    v[2] = -v[2];
}


/**
 * @brief Calculate the eigenvalues and eigenvectors of a 3x3 symmetric matrix using the Jacobi rotation method.
 * @param a input 3x3 symmetric matrix (the content will be modified)
 * @param d  output array of eigenvalues (with a size of 3)
 * @param v  output eigenvectors matrix (3x3), and the eigenvectors are column vectors.
 * @return 0 indicates successful convergence
 */
inline int jacobi3(rbmd::Real a[3][3], rbmd::Real d[3], rbmd::Real v[3][3]) {
    constexpr int max_rotations = 50;

    // Initialize the feature vector matrix v to be the identity matrix
    v[0][0] = v[1][1] = v[2][2] = 1.0;
    v[0][1] = v[0][2] = v[1][0] = v[1][2] = v[2][0] = v[2][1] = 0.0;

    // Initialize the characteristic value d and the intermediate vector b
    d[0] = a[0][0];
    d[1] = a[1][1];
    d[2] = a[2][2];
    rbmd::Real b[3] = {d[0], d[1], d[2]};
    rbmd::Real z[3] = {0.0, 0.0, 0.0};

    for (int iter = 0; iter < max_rotations; ++iter) {
        // Check for convergence: If the sum of all off-diagonal elements is sufficiently small
        rbmd::Real sum_off_diag = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
        if (sum_off_diag < 1.0e-12) {
            return 0; // Successful Contraction
        }

        // Select the element with the largest absolute value that is not on the diagonal for rotation.
        int p = 0, q = 1;
        if (std::fabs(a[0][2]) > std::fabs(a[p][q])) { p = 0; q = 2; }
        if (std::fabs(a[1][2]) > std::fabs(a[p][q])) { p = 1; q = 2; }

        rbmd::Real t, tau, c, s;
        if (std::fabs(a[p][q]) < 1.0e-12) {
            t = 0.0;
        } else {
            tau = (d[q] - d[p]) / (2.0 * a[p][q]);
            t = 1.0 / (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
            if (tau < 0.0) t = -t;
        }

        c = 1.0 / std::sqrt(1.0 + t * t);
        s = t * c;

        // Update matrix a, eigenvalue d, and eigenvector v
        rbmd::Real h = t * a[p][q];
        z[p] -= h;
        z[q] += h;
        d[p] -= h;
        d[q] += h;
        a[p][q] = 0.0;

        for (int j = 0; j < p; ++j) {
            rbmd::Real gj = a[j][p];
            rbmd::Real hj = a[j][q];
            a[j][p] = c * gj - s * hj;
            a[j][q] = s * gj + c * hj;
        }
        for (int j = p + 1; j < q; ++j) {
            rbmd::Real gpj = a[p][j];
            rbmd::Real hj = a[j][q];
            a[p][j] = c * gpj - s * hj;
            a[j][q] = s * gpj + c * hj;
        }
        for (int j = q + 1; j < 3; ++j) {
            rbmd::Real gpj = a[p][j];
            rbmd::Real gqj = a[q][j];
            a[p][j] = c * gpj - s * gqj;
            a[q][j] = s * gpj + c * gqj;
        }
        for (int j = 0; j < 3; ++j) {
            rbmd::Real vjp = v[j][p];
            rbmd::Real vjq = v[j][q];
            v[j][p] = c * vjp - s * vjq;
            v[j][q] = s * vjp + c * vjq;
        }
    }
    return 1; // Reaching the maximum number of iterations and still not converging
}


/**
 * @brief Calculate the angular velocity based on the angular momentum and the diagonalized moment of inertia tensor.
 * @param angmom Angular momentum L (space frame)
 * @param ex The first main axis (feature vector)
 * @param ey The second main axis (feature vector)
 * @param ez The third main axis (feature vector)
 * @param idiag Principal moment of inertia (eigenvalue)
 * @param w The output angular velocity ω (space frame)
 */
inline void angmom_to_omega(const rbmd::Real* angmom,
                            const rbmd::Real* ex, const rbmd::Real* ey, const rbmd::Real* ez,
                            const rbmd::Real* idiag, rbmd::Real* w) {
    rbmd::Real wbody[3];

    // Project the angular momentum in the spatial coordinate system
    // onto the object coordinate system (the principal axis coordinate system)
    // w_body = I_body^-1 * L_body = (1/I_diag) * (P^T * L_space)
    if (std::fabs(idiag[0]) < 1.0e-12) {
        wbody[0] = 0.0;
    } else {
        wbody[0] = dot3(angmom, ex) / idiag[0];
    }

    if (std::fabs(idiag[1]) < 1.0e-12) {
        wbody[1] = 0.0;
    } else {
        wbody[1] = dot3(angmom, ey) / idiag[1];
    }

    if (std::fabs(idiag[2]) < 1.0e-12) {
        wbody[2] = 0.0;
    } else {
        wbody[2] = dot3(angmom, ez) / idiag[2];
    }

    // Convert the angular velocity in the object coordinate system back to the spatial coordinate system
    // w_space = P * w_body
    w[0] = wbody[0] * ex[0] + wbody[1] * ey[0] + wbody[2] * ez[0];
    w[1] = wbody[0] * ex[1] + wbody[1] * ey[1] + wbody[2] * ez[1];
    w[2] = wbody[0] * ex[2] + wbody[1] * ey[2] + wbody[2] * ez[2];
}

/* ----------------------------------------------------------------------
   factorial n vial lookup from precomputed table
------------------------------------------------------------------------- */

inline  rbmd::Real factorial(const rbmd::Id n)
{
  if (n < 0 || n > nmaxfactorial)
    return std::numeric_limits<double>::quiet_NaN();

  return nfac_table[n];
}

} // namespace MathLib

