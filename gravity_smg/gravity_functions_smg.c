/** @file gravity_functions_smg.c Documented gravity_functions_smg module
 *
 * Emilio Bellini, Ignacy Sawicki, Miguel Zumalacarregui, 2024
 *
 * This module contains all the complicated expressions
 * used in hi_class. In particular, they are casted in
 * different functions, depending on their type:
 * - background (first Es and then Ps and Rs) as a functions of the Gs
 * - alphas as a functions of the Gs
 * - As as a functions of the alphas
 * - Bs as a functions of the Gs
 * - Cs as a functions of the As
 *
 * The purpose of this module is twofold:
 * - isolate long expressions in one single file
 * - have a file that can be optimized for specific theories.
 *   Indeed, the problem of numerical errors can be alleviated
 *   when analytic cancellations greatly simplify long expressions
 */

#include "gravity_functions_smg.h"


/**
* Get gravity functions Es from Gs. These are the functions
* entering the Friedmann time-time constraint, which is
* cubic in H for beyond Horndeski theories. In background_smg
* these functions are used to solve algebrically for H. The form
* is the following:
*
* E0 + E1 H + E3 H^3 = E2 H^2
*
* NOTE: H is NOT the curly H!!
* NOTE: added rho_tot and p_tot do not contain _smg
*
* @param pba           Input: pointer to background structure
* @param a             Input: scale factor
* @param pvecback_B    Input: vector containing all {B} quantities
* @param pvecback      Output: vector of background quantities (assumed to be already allocated)
* @param pgf           Input: pointer to G_functions_and_derivs structure
* @return the error status
*/
int gravity_functions_Es_from_Gs_smg(
                                     struct background *pba,
                                     double a,
                                     double * pvecback_B,
                                     double * pvecback,
                                     struct G_functions_and_derivs *pgf
                                     ) {

  /* Define local variables */
  double phi = pvecback_B[pba->index_bi_phi_smg];
  double phi_prime = pvecback_B[pba->index_bi_phi_prime_smg];
  double X = 0.5*pow(phi_prime/a,2);
  double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];

  double G2=pgf->G2;
  double G2_X=pgf->G2_X, G2_XX=pgf->G2_XX, G2_XXX=pgf->G2_XXX;
  double G2_phi=pgf->G2_phi, G2_Xphi=pgf->G2_Xphi, G2_XXphi=pgf->G2_XXphi;
  double G2_phiphi=pgf->G2_phiphi, G2_Xphiphi=pgf->G2_Xphiphi;
  double G3_X=pgf->G3_X, G3_XX=pgf->G3_XX, G3_XXX=pgf->G3_XXX;
  double G3_phi=pgf->G3_phi, G3_Xphi=pgf->G3_Xphi, G3_XXphi=pgf->G3_XXphi;
  double G3_phiphi=pgf->G3_phiphi, G3_Xphiphi=pgf->G3_Xphiphi;
  double G3_phiphiphi=pgf->G3_phiphiphi;
  double G4=pgf->G4, DG4=pgf->DG4;
  double G4_X=pgf->G4_X, G4_XX=pgf->G4_XX, G4_XXX=pgf->G4_XXX, G4_XXXX=pgf->G4_XXXX;
  double G4_phi=pgf->G4_phi, G4_Xphi=pgf->G4_Xphi, G4_XXphi=pgf->G4_XXphi, G4_XXXphi=pgf->G4_XXXphi;
  double G4_phiphi=pgf->G4_phiphi, G4_Xphiphi=pgf->G4_Xphiphi, G4_XXphiphi=pgf->G4_XXphiphi;
  double G4_phiphiphi=pgf->G4_phiphiphi, G4_Xphiphiphi=pgf->G4_Xphiphiphi;
  double G5_X=pgf->G5_X, G5_XX=pgf->G5_XX, G5_XXX=pgf->G5_XXX, G5_XXXX=pgf->G5_XXXX;
  double G5_phi=pgf->G5_phi, G5_Xphi=pgf->G5_Xphi, G5_XXphi=pgf->G5_XXphi, G5_XXXphi=pgf->G5_XXXphi;
  double G5_phiphi=pgf->G5_phiphi, G5_Xphiphi=pgf->G5_Xphiphi, G5_XXphiphi=pgf->G5_XXphiphi;
  double G5_phiphiphi=pgf->G5_phiphiphi, G5_Xphiphiphi=pgf->G5_Xphiphiphi;
  double F4=pgf->F4;
  double F4_X=pgf->F4_X, F4_XX=pgf->F4_XX, F4_XXX=pgf->F4_XXX;
  double F4_phi=pgf->F4_phi, F4_Xphi=pgf->F4_Xphi, F4_XXphi=pgf->F4_XXphi;
  double F4_phiphi=pgf->F4_phiphi, F4_Xphiphi=pgf->F4_Xphiphi;
  double F5=pgf->F5;
  double F5_X=pgf->F5_X, F5_XX=pgf->F5_XX, F5_XXX=pgf->F5_XXX;
  double F5_phi=pgf->F5_phi, F5_Xphi=pgf->F5_Xphi, F5_XXphi=pgf->F5_XXphi;
  double F5_phiphi=pgf->F5_phiphi, F5_Xphiphi=pgf->F5_Xphiphi;

  pvecback[pba->index_bg_E0_smg] =
  1./3.*(
    + 3.*rho_tot - G2 + 2.*X*(G2_X - G3_phi)
  );

  pvecback[pba->index_bg_E1_smg] =
  - 2.*(
    + G4_phi - X*(G3_X - 2.*G4_Xphi)
  )*phi_prime/a;

  pvecback[pba->index_bg_E2_smg] =
  2.*(
    1./2. + DG4 - X*(4.*G4_X - 3.*G5_phi)
    - 2.*pow(X,2)*(2.*G4_XX - G5_Xphi + 10.*F4) - 8.*pow(X,3)*F4_X
  );

  pvecback[pba->index_bg_E3_smg] =
  2./3.*X*(
    + 5.*G5_X + 2.*X*(G5_XX - 42.*F5) - 24.*pow(X,2)*F5_X
  )*phi_prime/a;

  return _SUCCESS_;
}


/**
* Get gravity functions Ps and Rs from Gs. These are the functions
* entering the Friedmann space-space equation and the scalar field
* one. In background_smg these functions are used to solve algebrically
* for H' and phi''. The form is the following:
*
* P0 phi'' + P1 H' + P2 = 0
* R0 phi'' + R1 H' + R2 = 0
*
* NOTE: H is NOT the curly H!!
* NOTE: added rho_tot and p_tot do not contain _smg
*
* @param pba           Input: pointer to background structure
* @param a             Input: scale factor
* @param pvecback_B    Input: vector containing all {B} quantities
* @param pvecback      Output: vector of background quantities (assumed to be already allocated)
* @param pgf           Input: pointer to G_functions_and_derivs structure
* @return the error status
*/
int gravity_functions_Ps_and_Rs_from_Gs_smg(
                                            struct background *pba,
                                            double a,
                                            double * pvecback_B,
                                            double * pvecback,
                                            struct G_functions_and_derivs *pgf
                                            ) {

  /* Define local variables */
  double H = pvecback[pba->index_bg_H];
  double phi = pvecback_B[pba->index_bi_phi_smg];
  double phi_prime = pvecback_B[pba->index_bi_phi_prime_smg];
  double X = 0.5*pow(phi_prime/a,2);
  double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];

  double G2=pgf->G2;
  double G2_X=pgf->G2_X, G2_XX=pgf->G2_XX, G2_XXX=pgf->G2_XXX;
  double G2_phi=pgf->G2_phi, G2_Xphi=pgf->G2_Xphi, G2_XXphi=pgf->G2_XXphi;
  double G2_phiphi=pgf->G2_phiphi, G2_Xphiphi=pgf->G2_Xphiphi;
  double G3_X=pgf->G3_X, G3_XX=pgf->G3_XX, G3_XXX=pgf->G3_XXX;
  double G3_phi=pgf->G3_phi, G3_Xphi=pgf->G3_Xphi, G3_XXphi=pgf->G3_XXphi;
  double G3_phiphi=pgf->G3_phiphi, G3_Xphiphi=pgf->G3_Xphiphi;
  double G3_phiphiphi=pgf->G3_phiphiphi;
  double G4=pgf->G4, DG4=pgf->DG4;
  double G4_X=pgf->G4_X, G4_XX=pgf->G4_XX, G4_XXX=pgf->G4_XXX, G4_XXXX=pgf->G4_XXXX;
  double G4_phi=pgf->G4_phi, G4_Xphi=pgf->G4_Xphi, G4_XXphi=pgf->G4_XXphi, G4_XXXphi=pgf->G4_XXXphi;
  double G4_phiphi=pgf->G4_phiphi, G4_Xphiphi=pgf->G4_Xphiphi, G4_XXphiphi=pgf->G4_XXphiphi;
  double G4_phiphiphi=pgf->G4_phiphiphi, G4_Xphiphiphi=pgf->G4_Xphiphiphi;
  double G5_X=pgf->G5_X, G5_XX=pgf->G5_XX, G5_XXX=pgf->G5_XXX, G5_XXXX=pgf->G5_XXXX;
  double G5_phi=pgf->G5_phi, G5_Xphi=pgf->G5_Xphi, G5_XXphi=pgf->G5_XXphi, G5_XXXphi=pgf->G5_XXXphi;
  double G5_phiphi=pgf->G5_phiphi, G5_Xphiphi=pgf->G5_Xphiphi, G5_XXphiphi=pgf->G5_XXphiphi;
  double G5_phiphiphi=pgf->G5_phiphiphi, G5_Xphiphiphi=pgf->G5_Xphiphiphi;
  double F4=pgf->F4;
  double F4_X=pgf->F4_X, F4_XX=pgf->F4_XX, F4_XXX=pgf->F4_XXX;
  double F4_phi=pgf->F4_phi, F4_Xphi=pgf->F4_Xphi, F4_XXphi=pgf->F4_XXphi;
  double F4_phiphi=pgf->F4_phiphi, F4_Xphiphi=pgf->F4_Xphiphi;
  double F5=pgf->F5;
  double F5_X=pgf->F5_X, F5_XX=pgf->F5_XX, F5_XXX=pgf->F5_XXX;
  double F5_phi=pgf->F5_phi, F5_Xphi=pgf->F5_Xphi, F5_XXphi=pgf->F5_XXphi;
  double F5_phiphi=pgf->F5_phiphi, F5_Xphiphi=pgf->F5_Xphiphi;

  pvecback[pba->index_bg_P0_smg] =
  - 2./3.*(
    + G4_phi - X*(G3_X - 2.*G4_Xphi)
    - 2.*(
      + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
    )*H*phi_prime/a
    - X*(
      + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
    )*pow(H,2)
  )/a;

  pvecback[pba->index_bg_P1_smg] =
  - 2./3.*(
    + 1. + 2.*DG4 - 2.*X*(2.*G4_X - G5_phi)
    - 8.*F4*pow(X,2) - 2.*X*(G5_X - 12.*X*F5 )*H*phi_prime/a
  );

  pvecback[pba->index_bg_P2_smg] =
  - 1./3.*a*(
    + 3.*(rho_tot + p_tot)
    + 2.*X*(G2_X - 2.*G3_phi + 2.*G4_phiphi)
    - 4.*(
      + G4_phi - X*(2.*G3_X - 6.*G4_Xphi + G5_phiphi) + 4.*pow(X,2)*F4_phi
    )*H*phi_prime/a
    + 4.*X*(
      + 5.*(G4_X - G5_phi) + 2.*X*(5.*G4_XX - 3.*G5_Xphi + 20.*F4)
      + 4.*pow(X,2)*(5.*F4_X + 3.*F5_phi)
    )*pow(H,2)
    + 4.*X*(
      + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
    )*pow(H,3)*phi_prime/a
  );

  pvecback[pba->index_bg_R0_smg] =
  1./3.*(
    + (
      + G2_X - 2.*G3_phi + 2.*X*(G2_XX - G3_Xphi)
    )/a/H
    + 6.*(
      + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
    )*pow(a,-2)*phi_prime
    + 2.*(
      + 3.*G5_X + X*(7.*G5_XX - 120.*F5)
      + 2.*pow(X,2)*(G5_XXX - 66.*F5_X) - 24.*pow(X,3)*F5_XX
    )*pow(a,-2)*pow(H,2)*phi_prime
    + 6.*(
      + G4_X - G5_phi + X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
      + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X) + 8.*pow(X,3)*F4_XX
    )*H/a
  );

  pvecback[pba->index_bg_R1_smg] =
  2.*(
    - (G4_phi - X*(G3_X - 2.*G4_Xphi))/H
    + 2.*(
      + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
    )*phi_prime/a
    + X*(
      + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
    )*H
  );

  pvecback[pba->index_bg_R2_smg] =
  1./3.*(
    - (G2_phi - 2.*X*(G2_Xphi - G3_phiphi))*a/H
    + 2.*(
      + G2_X - 2.*G3_phi - X*(G2_XX - 4.*G3_Xphi + 6.*G4_Xphiphi)
    )*phi_prime
    - 2.*(
      + 6.*G4_phi - 3.*X*(G3_X - G5_phiphi)
      + 6.*pow(X,2)*(G3_XX - 4.*G4_XXphi + G5_Xphiphi - 6.*F4_phi)
      - 24.*pow(X,3)*F4_Xphi
    )*a*H
    + 4.*(
      + 3.*(G4_X - G5_phi) - X*(3.*G4_XX - 4.*G5_Xphi)
      - 2.*pow(X,2)*(3.*G4_XXX - 2.*G5_XXphi + 18.*F4_X + 12.*F5_phi)
      - 12.*pow(X,3)*(F4_XX + F5_Xphi)
    )*pow(H,2)*phi_prime
    + 2.*X*(
      + 3.*G5_X - 4.*X*(2.*G5_XX - 15.*F5)
      - 4.*pow(X,2)*(G5_XXX - 48.*F5_X) + 48.*pow(X,3)*F5_XX
    )*a*pow(H,3)
  );

  return _SUCCESS_;
}


/**
* Get building blocks from Gs. These are the basic functions
* used to build an effective theory for dark energy. They are:
* - density and pressure of the background _smg field
* - an effective Planck mass and the alphas (K, B, M, T, H)
* - current and shift of the background _smg field (these are not
*   passed to the perturbation module, as they are dependent on
*   the others)
*
* NOTE: H is NOT the curly H!!
* NOTE: added rho_tot and p_tot do not contain _smg
*
* @param pba           Input: pointer to background structure
* @param a             Input: scale factor
* @param pvecback_B    Input: vector containing all {B} quantities
* @param pvecback      Output: vector of background quantities (assumed to be already allocated)
* @param pgf           Input: pointer to G_functions_and_derivs structure
* @return the error status
*/
int gravity_functions_building_blocks_from_Gs_smg(
                                                  struct background *pba,
                                                  double a,
                                                  double * pvecback_B,
                                                  double * pvecback,
                                                  struct G_functions_and_derivs *pgf
                                                  ) {

  /* Define local variables */
  double H = pvecback[pba->index_bg_H];
  double phi = pvecback_B[pba->index_bi_phi_smg];
  double phi_prime = pvecback_B[pba->index_bi_phi_prime_smg];
  double X = 0.5*pow(phi_prime/a,2);

  double G2=pgf->G2;
  double G2_X=pgf->G2_X, G2_XX=pgf->G2_XX, G2_XXX=pgf->G2_XXX;
  double G2_phi=pgf->G2_phi, G2_Xphi=pgf->G2_Xphi, G2_XXphi=pgf->G2_XXphi;
  double G2_phiphi=pgf->G2_phiphi, G2_Xphiphi=pgf->G2_Xphiphi;
  double G3_X=pgf->G3_X, G3_XX=pgf->G3_XX, G3_XXX=pgf->G3_XXX;
  double G3_phi=pgf->G3_phi, G3_Xphi=pgf->G3_Xphi, G3_XXphi=pgf->G3_XXphi;
  double G3_phiphi=pgf->G3_phiphi, G3_Xphiphi=pgf->G3_Xphiphi;
  double G3_phiphiphi=pgf->G3_phiphiphi;
  double G4=pgf->G4, DG4=pgf->DG4;
  double G4_X=pgf->G4_X, G4_XX=pgf->G4_XX, G4_XXX=pgf->G4_XXX, G4_XXXX=pgf->G4_XXXX;
  double G4_phi=pgf->G4_phi, G4_Xphi=pgf->G4_Xphi, G4_XXphi=pgf->G4_XXphi, G4_XXXphi=pgf->G4_XXXphi;
  double G4_phiphi=pgf->G4_phiphi, G4_Xphiphi=pgf->G4_Xphiphi, G4_XXphiphi=pgf->G4_XXphiphi;
  double G4_phiphiphi=pgf->G4_phiphiphi, G4_Xphiphiphi=pgf->G4_Xphiphiphi;
  double G5_X=pgf->G5_X, G5_XX=pgf->G5_XX, G5_XXX=pgf->G5_XXX, G5_XXXX=pgf->G5_XXXX;
  double G5_phi=pgf->G5_phi, G5_Xphi=pgf->G5_Xphi, G5_XXphi=pgf->G5_XXphi, G5_XXXphi=pgf->G5_XXXphi;
  double G5_phiphi=pgf->G5_phiphi, G5_Xphiphi=pgf->G5_Xphiphi, G5_XXphiphi=pgf->G5_XXphiphi;
  double G5_phiphiphi=pgf->G5_phiphiphi, G5_Xphiphiphi=pgf->G5_Xphiphiphi;
  double F4=pgf->F4;
  double F4_X=pgf->F4_X, F4_XX=pgf->F4_XX, F4_XXX=pgf->F4_XXX;
  double F4_phi=pgf->F4_phi, F4_Xphi=pgf->F4_Xphi, F4_XXphi=pgf->F4_XXphi;
  double F4_phiphi=pgf->F4_phiphi, F4_Xphiphi=pgf->F4_Xphiphi;
  double F5=pgf->F5;
  double F5_X=pgf->F5_X, F5_XX=pgf->F5_XX, F5_XXX=pgf->F5_XXX;
  double F5_phi=pgf->F5_phi, F5_Xphi=pgf->F5_Xphi, F5_XXphi=pgf->F5_XXphi;
  double F5_phiphi=pgf->F5_phiphi, F5_Xphiphi=pgf->F5_Xphiphi;

  /* Computing density, pressure (also rewritten as shift and current) */

  /* Energy density of the field */
  pvecback[pba->index_bg_rho_smg] =
    - (G2 - 2.*X*(G2_X - G3_phi))/3.
    - 2.*(G4_phi - X*(G3_X - 2.*G4_Xphi))*H*phi_prime/a
    - 2.*(
      + DG4 - X*(4.*G4_X - 3.*G5_phi) - 2.*pow(X,2)*(2.*G4_XX - G5_Xphi + 10.*F4)
      - 8.*pow(X,3)*F4_X
    )*pow(H,2)
    + 2./3.*X*(
      + 5.*G5_X + 2.*X*(G5_XX - 42.*F5) - 24.*pow(X,2)*F5_X
    )*pow(H,3)*phi_prime/a;

  /* Pressure of the field */
  pvecback[pba->index_bg_p_smg] =
  (
    + G2 - 2.*X*(G3_phi - 2.*G4_phiphi)
    + 2.*(
      + 3.*DG4 - X*(2.*G4_X + G5_phi)
      + 2.*pow(X,2)*(4.*G4_XX - 3.*G5_Xphi + 10.*F4)
      + 8.*pow(X,3)*(2.*F4_X + 3.*F5_phi)
    )*pow(H,2)
    + 2.*(
      + G4_phi + X*(G3_X - 6.*G4_Xphi + 2.*G5_phiphi) - 8.*pow(X,2)*F4_phi
    )*H*phi_prime/a
    + 2.*X*(
      + G5_X + 2.*X*(G5_XX - 18.*F5) - 24.*pow(X,2)*F5_X
    )*pow(H,3)*phi_prime/a
    + 4.*(
      + DG4 - X*(2.*G4_X - G5_phi) - 4.*pow(X,2)*F4
      - X*(G5_X - 12.*X*F5)*H*phi_prime/a
    )*pvecback[pba->index_bg_H_prime]/a
    + 2.*(
      + G4_phi - X*(G3_X - 2.*G4_Xphi)
      - X*(3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X)*pow(H,2)
      - 2.*(
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
      )*H*phi_prime/a
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/3.;

  /* Current of the field */
  pvecback[pba->index_bg_current_smg] =
    + (G2_X - G3_phi)*phi_prime/a + 6.*X*(G3_X - 2.*G4_Xphi)*H
    + 3.*(
      + 2.*G4_X - G5_phi + 2.*X*(2.*G4_XX - G5_Xphi + 8.*F4)
      + 8.*pow(X,2)*F4_X
    )*pow(H,2)*phi_prime/a
    + 2.*X*(
      + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
    )*pow(H,3);

  /* Shift of the field */
  pvecback[pba->index_bg_shift_smg] =
    + G2_phi + 2.*G3_phi*H*phi_prime/a
    + 12.*(G4_phi + 2.*pow(X,2)*F4_phi)*pow(H,2)
    + 2.*(
      + 3.*G5_phi - 2.*X*G5_Xphi - 12.*pow(X,2)*F5_phi
    )*pow(H,3)*phi_prime/a
    + 6.*(G4_phi + G5_phi*H*phi_prime/a)*pvecback[pba->index_bg_H_prime]/a
    + (
      + G3_phi + 6.*G4_Xphi*H*phi_prime/a
      + 3.*(G5_phi + 2.*X*G5_Xphi)*pow(H,2)
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg];


  /* Computing alphas at the end (alpha_T, alpha_M depend on phi'') */

  /* Planck mass */
  pvecback[pba->index_bg_delta_M2_smg] =
    + 2.*(DG4 - X*(2.*G4_X - G5_phi) - 4.*F4*pow(X,2))
    - 2.*X*(G5_X - 12.*X*F5)*H*phi_prime/a;

  pvecback[pba->index_bg_M2_smg] = 1. + pvecback[pba->index_bg_delta_M2_smg];

  /* alpha_K kineticity */
  pvecback[pba->index_bg_kineticity_over_phiphi_smg] =
  (
    + G2_X - 2.*G3_phi + 2.*X*(G2_XX - G3_Xphi)
    + 6.*(
      + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
    )*H*phi_prime/a
    + 2.*(
      + 3.*G5_X + X*(7.*G5_XX - 120.*F5)
      + 2.*pow(X,2)*(G5_XXX - 66.*F5_X) - 24.*pow(X,3)*F5_XX
    )*pow(H,3)*phi_prime/a
    + 6.*(
      + G4_X - G5_phi + X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
      + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X) + 8.*pow(X,3)*F4_XX
    )*pow(H,2)
  )/pvecback[pba->index_bg_M2_smg];

  pvecback[pba->index_bg_kineticity_smg] = 2.*pow(H,-2)*X*pvecback[pba->index_bg_kineticity_over_phiphi_smg];

  /* alpha_B braiding */
  pvecback[pba->index_bg_braiding_over_phi_smg] =
  2.*(
    - G4_phi + X*(G3_X - 2.*G4_Xphi)
    + X*(
      + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
    )*pow(H,2)
    + 2.*(
      + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
    )*H*phi_prime/a
  )/pvecback[pba->index_bg_M2_smg];

  pvecback[pba->index_bg_braiding_smg] = phi_prime/a/H*pvecback[pba->index_bg_braiding_over_phi_smg];

  /* alpha_T: tensor speed excess */
  pvecback[pba->index_bg_tensor_excess_smg] =
  2.*X*(
    + 2*(G4_X - G5_phi + 2.*X*F4)
    + 2*(G5_X - 6.*X*F5)*H*phi_prime/a
    - G5_X*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* alpha_M: Planck mass running */
  pvecback[pba->index_bg_M2_running_smg] =
  2.*(
    + 2.*X*(
      + G4_X - G5_phi + 2.*X*(G4_XX - G5_Xphi + 4.*F4)
      + 4.*pow(X,2)*(F4_X + 3.*F5_phi)
    )
    + (G4_phi - X*(2.*G4_Xphi - G5_phiphi) - 4.*pow(X,2)*F4_phi)*phi_prime/a/H
    + X*(3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X)*H*phi_prime/a
    - X*(G5_X - 12.*X*F5)*pvecback[pba->index_bg_H_prime]*pow(a,-2)*phi_prime/H
    + (
      - 3.*X*G5_X
      + 2.*pow(X,2)*(30.*F5 - G5_XX)
      + 24.*pow(X,3)*F5_X
      - (
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
      )*phi_prime/a/H
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* alpha_H: Beyond Horndeski */
  pvecback[pba->index_bg_beyond_horndeski_over_phi_smg] =
  4.*X*(
    F4*phi_prime/a - 6.*X*H*F5
  )*H/pvecback[pba->index_bg_M2_smg];

  pvecback[pba->index_bg_beyond_horndeski_smg] = pvecback[pba->index_bg_beyond_horndeski_over_phi_smg]*phi_prime/a/H;

  return _SUCCESS_;
}


/**
* Get gravity functions As from alphas.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Input/Output: vector of background quantities
* @param pvecback_derivs      Input: vector of derivatives
* @return the error status
*/
int gravity_functions_As_from_alphas_smg(
                                         struct background *pba,
                                         double * pvecback,
                                         double * pvecback_derivs
                                         ) {

  // basic background quantities
  double a = pvecback[pba->index_bg_a];
  double H = pvecback[pba->index_bg_H];
  double H_p = pvecback[pba->index_bg_H_prime];
  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];
  double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
  double p_smg = pvecback[pba->index_bg_p_smg];
  double rho_smg = pvecback[pba->index_bg_rho_smg];
  // factor to convert lna derivatives to tau derivatives
  double factor = a*H;

  // alphas
  double M2 = pvecback[pba->index_bg_M2_smg];
  double DelM2 = pvecback[pba->index_bg_delta_M2_smg];
  double kin = pvecback[pba->index_bg_kineticity_smg];
  double bra = pvecback[pba->index_bg_braiding_smg];
  double run = pvecback[pba->index_bg_M2_running_smg];
  double ten = pvecback[pba->index_bg_tensor_excess_smg];
  double beh = pvecback[pba->index_bg_beyond_horndeski_smg];
  double dM2 = pvecback[pba->index_bg_delta_M2_smg];

  // need to update the time derivatives of the interesting functions
  double kin_p = factor*pvecback_derivs[pba->index_bg_kineticity_smg];
  double bra_p = factor*pvecback_derivs[pba->index_bg_braiding_smg];
  double run_p = factor*pvecback_derivs[pba->index_bg_M2_running_smg];
  double ten_p = factor*pvecback_derivs[pba->index_bg_tensor_excess_smg];
  double beh_p = factor*pvecback_derivs[pba->index_bg_beyond_horndeski_smg];
  double p_tot_p = factor*pvecback_derivs[pba->index_bg_p_tot_wo_smg];
  double p_smg_p = factor*pvecback_derivs[pba->index_bg_p_smg];

  // kinetic term D
  pvecback[pba->index_bg_kinetic_D_smg] = kin + 3./2.*pow(bra,2);

  // A0
	pvecback[pba->index_bg_A0_smg] =
	1./2.*(
    + bra - 3.*(rho_smg + p_smg + (rho_tot + p_tot)*DelM2/M2)*pow(H,-2)
	);

  // A1
	pvecback[pba->index_bg_A1_smg] =
	+ (1. + ten)*kin
	- 3.*(beh*(1. + run) + run - ten + beh_p/a/H)*bra;

  // A2
	pvecback[pba->index_bg_A2_smg] =
	- (kin + 3./2.*pow(bra,2))*(2. + run)
	- 9./4.*bra*(
	  + (2. - bra)*(rho_smg+p_smg)
	  + (2.*DelM2/M2 - bra)*(rho_tot+p_tot)
	)*pow(H,-2)
	- 3./2.*bra*bra_p/a/H;

  // A3
	pvecback[pba->index_bg_A3_smg] = bra*beh;

  // A4
  pvecback[pba->index_bg_A4_smg] =
  9./2.*kin*(
    + (2. - bra)*(rho_smg+p_smg)
    + (2.*DelM2/M2 - bra)*(rho_tot+p_tot)
  )*pow(H,-2)
  + 3.*(bra*kin_p - kin*bra_p)/a/H;

  // A5
	pvecback[pba->index_bg_A5_smg] = - beh*kin;

  // A6
	pvecback[pba->index_bg_A6_smg] =
	+ 9./4.*(
	  + (2.*kin + 9.*bra)*(2.*DelM2/M2 - bra)
	  + 4.*(kin + 3./2.*pow(bra,2))*run
	)
	+ 9.*(kin + 9./2.*bra)*rho_smg*pow(H,-2)/M2
	+ 9./2.*(
	  + (kin + 9.*bra)*(2.*DelM2/M2 - bra)
	  + 2.*(kin + 3./2.*pow(bra,2))*run
	)*pow(H,-2)*p_tot
	+ 81./4.*bra*(
	  + 2.*rho_smg*(p_tot + p_smg)/M2
	  - 2.*(1./M2 - 2. + bra)*p_tot*p_smg
	  + (2. - bra)*pow(p_smg,2)
	  + (2.*DelM2 - bra*M2)*pow(p_tot,2)/M2
	)*pow(H,-4)
	+ 9./2.*(
	  - 9.*bra*(1./M2 - 2. + bra)
	  + kin*(2. - bra)
	  + 2.*(kin + 3./2.*pow(bra,2))*run
	)*pow(H,-2)*p_smg
	+ 3.*(
	  + bra*kin_p
	  - (kin - 9./2.*bra - 9./2.*bra*pow(H,-2)*(p_tot + p_smg))*bra_p
	)/a/H
	+ 9.*(
	  + (kin*DelM2/M2 + 3./2.*pow(bra,2))*p_tot_p
	  + (kin + 3./2.*pow(bra,2))*p_smg_p
	)*pow(H,-3)/a;

  // A7
	pvecback[pba->index_bg_A7_smg] =
	- 2.*kin*beh
	+ 3.*bra*(bra + 2.*beh)*(1. + run)
	+ 2.*(kin + 3.*bra)*(run - ten)
	+ 9./2.*bra*(
	  + (2. - bra - 2.*beh)*(rho_smg + p_smg)
	  + (2.*DelM2/M2 - bra - 2.*beh)*(rho_tot + p_tot)
	)*pow(H,-2)
	+ 3.*bra*(bra_p + 2.*beh_p)/a/H;

  // A8
	pvecback[pba->index_bg_A8_smg] = run - ten - beh;

  // A9
	pvecback[pba->index_bg_A9_smg] =
	+ 3./4.*(
	  + (2. - bra)*(rho_smg + p_smg)
	  + (2.*DelM2/M2 - bra)*(rho_tot + p_tot)
	)*pow(H,-2)
	+ 1./2.*bra_p/a/H;

  // A10
	pvecback[pba->index_bg_A10_smg] =
	bra + 2.*run - (2. - bra)*ten + 2.*(1. + run)*beh + 2.*beh_p/a/H;

  // A11
	pvecback[pba->index_bg_A11_smg] =
	- (kin + 3./2.*pow(bra,2))*(4. + run)
	+ 3./4.*(
	  + (4.*kin + 6.*bra + 3.*pow(bra,2))*(rho_smg + p_smg)
	  + (4.*kin + 6.*bra*DelM2/M2 + 3.*pow(bra,2))*(rho_tot + p_tot)
	)*pow(H,-2)
	- (kin_p + 3./2.*bra*bra_p)/a/H;

  // A12
	pvecback[pba->index_bg_A12_smg] =
	+ kin/2. - 3.*bra*(3./2. - bra) - run*(kin + 3./2.*pow(bra,2))
	- 9./4.*(
		+ (6.*DelM2/M2 + bra*(2./M2 - 7. + 2.*bra))*pow(H,-2)*rho_tot
		+ (6.*DelM2/M2  - 2.*kin + bra*(2./M2 - 5. - 2.*bra))*pow(H,-2)*p_tot
		+ (6. - 5.*bra - 2.*pow(bra,2) - 2.*kin)*pow(H,-2)*p_smg
		- 6.*(1./M2 - 2. + bra)*pow(H,-4)*p_tot*p_smg
		+ 3.*(2. - bra)*pow(H,-4)*pow(p_smg,2)
		+ 3.*(2.*DelM2/M2 - bra)*(
			+ rho_tot*p_tot + pow(p_tot,2) + rho_tot*p_smg
		)*pow(H,-4)
		+ (2. - bra)*(
			+ 3. - 2.*bra + 3.*pow(H,-2)*(p_tot + p_smg)
		)*pow(H,-2)*rho_smg
		+ 2.*bra*pow(H,-3)*p_tot_p/a/M2
	)
	- (
	 	+ kin_p
	 	+ 3./2.*(3. + bra + 3.*pow(H,-2)*(p_tot + p_smg))*bra_p
	)/a/H;

  // A13
	pvecback[pba->index_bg_A13_smg] =
	- bra - 2.*run + (2. - bra)*ten - (2. + bra + 2.*run)*beh
	- 3./2.*(
	 	+ (2. - bra - 2.*beh)*(rho_smg + p_smg)*pow(H,-2)
	 	+ (2.*DelM2/M2 - bra - 2.*beh)*(rho_tot + p_tot)*pow(H,-2)
	)
	- (bra_p + 2.*beh_p)/a/H;

  // A14
	pvecback[pba->index_bg_A14_smg] = - (kin + 3.*bra)/2.;

  // A15
	pvecback[pba->index_bg_A15_smg] = - 1./2.*bra - beh;

  // A16
	pvecback[pba->index_bg_A16_smg] =
	- 1./2.*(kin + 3.*bra)
	+ 9./4.*(
	 + (2. - bra)*(rho_smg + p_smg)
	 + (2.*DelM2/M2 - bra)*(rho_tot + p_tot)
	)*pow(H,-2);


  // NOTE: this could be written in terms of the A's and B's variables.

  pvecback[pba->index_bg_lambda_1_smg] = (run + (-1.)*ten)*(-3.)*bra + (1. + ten)*kin;

  pvecback[pba->index_bg_lambda_2_smg] = (- 2.*dM2 + bra*M2)*(rho_tot + p_tot)*(-3.)/2.*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)/2.*pow(H,-2) + pow(H,-1)*bra_p*pow(a,-1);

	pvecback[pba->index_bg_lambda_3_smg] = (2. + run)*(-1.)/2.*pvecback[pba->index_bg_kinetic_D_smg] + (-3.)/4.*bra*pvecback[pba->index_bg_lambda_2_smg];

	pvecback[pba->index_bg_lambda_4_smg] = kin*pvecback[pba->index_bg_lambda_2_smg] + (2.*kin*bra_p + (-1.)*bra*kin_p)*(-1.)*pow(H,-1)*pow(a,-1);

	pvecback[pba->index_bg_lambda_5_smg] = (bra + 2.*run + (-2.)*ten + bra*ten)*3./2.*bra + (run + (-1.)*ten)*pvecback[pba->index_bg_kinetic_D_smg] + 3./2.*bra*pvecback[pba->index_bg_lambda_2_smg];

	pvecback[pba->index_bg_lambda_6_smg] = 3./2.*(((9./2.*bra + kin)*dM2*pow(M2,-1) + (-9.)/4.*pow(bra,2) - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*pow(rho_tot,2) + ((9.*bra + kin)*dM2*pow(M2,-1) + (-9.)/2.*pow(bra,2) - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*rho_tot*p_tot + 9./2.*bra*(dM2 - M2*bra/2.)*pow(M2,-1)*pow(p_tot,2) + (kin*dM2*pow(M2,-1) - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*(rho_tot + p_tot)*rho_smg + ((kin - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*rho_smg + ((9.*bra + kin)*(2. - bra)/2. + pvecback[pba->index_bg_kinetic_D_smg]*run - 9./2.*bra*pow(M2,-1))*rho_tot + 9.*bra*(1. - bra/2. - pow(M2,-1)/2.)*p_tot)*(rho_smg + p_smg) + 9./2.*bra*(1. - bra/2.)*pow(rho_smg + p_smg,2))*pow(H,-4) + (((9.*bra*(rho_tot + p_tot) - 2.*kin*(rho_tot + rho_smg)) + (rho_smg + p_smg)*9.*bra)*bra_p/2. + (rho_tot + rho_smg)*bra*kin_p + (2.*dM2*kin + 3.*pow(bra,2)*M2)*3./2.*pow(M2,-1)*p_tot_p + 3.*pvecback[pba->index_bg_kinetic_D_smg]*p_smg_p)*pow(H,-3)*pow(a,-1)/2.;

	pvecback[pba->index_bg_lambda_7_smg] = ((-2.) + bra)*(4. + run)*(-1.)/8.*pvecback[pba->index_bg_kinetic_D_smg] + ((-2.)*(2. + dM2) + bra*M2)*(rho_tot + p_tot)*3./16.*pow(H,-2)*pvecback[pba->index_bg_kinetic_D_smg]*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*3./16.*pow(H,-2)*pvecback[pba->index_bg_kinetic_D_smg] + (pvecback[pba->index_bg_kinetic_D_smg]*bra_p + ((-2.) + bra)*((-3.)*bra*bra_p + (-1.)*kin_p))*1./8.*pow(H,-1)*pow(a,-1);

	pvecback[pba->index_bg_lambda_8_smg] = ((-2.) + bra)*(4. + run)*1./8.*pvecback[pba->index_bg_kinetic_D_smg] + 3./8.*(rho_tot + p_tot)*(((-9.)*bra + (-2.)*pvecback[pba->index_bg_kinetic_D_smg]*(3. + 2.*dM2 - bra*M2))*(-1.)/2. + (-rho_tot*dM2 - (p_smg + rho_smg*M2))*9.*pow(H,-2)*pow(M2,-1))*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)/8.*pow(H,-2)*pvecback[pba->index_bg_kinetic_D_smg] + (-2.*dM2 + bra*M2)*(rho_tot + p_tot)*(p_tot + p_smg)*27./16.*pow(H,-4)*pow(M2,-2) + ((-9.)*(rho_tot + p_tot) + (-6.)*bra*pow(H,2)*M2 + 3.*pow(bra,2)*pow(H,2)*M2 + (-1.)*pow(H,2)*pvecback[pba->index_bg_kinetic_D_smg]*M2)*1./8.*pow(H,-3)*pow(M2,-1)*bra_p*pow(a,-1) + ((-2.) + bra)*1./8.*pow(H,-1)*kin_p*pow(a,-1) + ((-2.) + bra)*9./16.*bra*pow(H,-3)*pow(M2,-1)*p_tot_p*pow(a,-1);

	pvecback[pba->index_bg_lambda_9_smg] = ((-2.) + 3.*bra)*pvecback[pba->index_bg_kinetic_D_smg] + 2.*pvecback[pba->index_bg_lambda_3_smg] + (pvecback[pba->index_bg_kinetic_D_smg] + (-1.)*pvecback[pba->index_bg_lambda_2_smg])*(((-3.) + 2.*bra)*(-3.)/2. + (p_tot + p_smg)*9./2.*pow(H,-2)) + (3.*bra*bra_p + kin_p)*(-1.)*pow(H,-1)*pow(a,-1) + (-9.)/2.*bra*pow(H,-3)*pow(M2,-1)*p_tot_p*pow(a,-1);

	pvecback[pba->index_bg_lambda_10_smg] = (pvecback[pba->index_bg_kinetic_D_smg] + (-1.)*pvecback[pba->index_bg_lambda_3_smg])*(-2.) + (3.*bra*dM2 + kin*M2)*(rho_tot + p_tot)*3.*pow(H,-2)*pow(M2,-1) + (3.*bra + kin)*(rho_smg + p_smg)*3.*pow(H,-2) + (-1.)*pow(H,-1)*kin_p*pow(a,-1);

  pvecback[pba->index_bg_lambda_11_smg] = bra + 2.*run - (2.-bra)*ten;


  pvecback[pba->index_bg_cs2num_smg] = 
    + (2. - bra)*(bra + 2.*beh + 2.*run + 2.*beh*run - 2.*ten + bra*ten)/2.
    + 3./2.*(2. - bra)*(1. + beh)*pow(H,-2)*(rho_smg + p_smg)
    - 3./2.*(1. + beh)*(2. + 2.*beh - 2.*M2 + bra*M2)*pow(H,-2)/M2*(rho_tot + p_tot)
    + (1. + beh)*bra_p/a/H
    + (2. - bra)*beh_p/a/H;

  // NOTE: this is to regularize cs2 when both the numerator and denominator are
  // below numerical precision
	if (pvecback[pba->index_bg_cs2num_smg] == pvecback[pba->index_bg_kinetic_D_smg]) {
		pvecback[pba->index_bg_cs2_smg] = 1.;
	}
	else {
		pvecback[pba->index_bg_cs2_smg] = pvecback[pba->index_bg_cs2num_smg]/pvecback[pba->index_bg_kinetic_D_smg];
	}

  // NOTE: Geff and slip are calculated for Horndeski. Consider extending them to
  // beyond Horndeski
	double beta_1 = (run + (-1.)*ten)*2. + (1. + ten)*bra;
	double beta_2 = 2.*beta_1 + (2. + (-2.)*M2 + bra*M2)*(rho_tot + p_tot)*(-3.)*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)*pow(H,-2) + 2.*pow(H,-1)*bra_p*pow(a,-1);

	if (bra*beta_1 == 0.) {
		pvecback[pba->index_bg_G_eff_smg] = 1./M2;
	}
	else {
		pvecback[pba->index_bg_G_eff_smg] = (1. - bra*beta_1*pow(bra*beta_1 - beta_2,-1))/M2;
	}

  if (2.*(run - ten)*beta_1 + ten*beta_2 == 0.) {
		pvecback[pba->index_bg_slip_eff_smg] = 1.;
	}
	else {
		pvecback[pba->index_bg_slip_eff_smg] = 1. - (2.*(run - ten)*beta_1 + ten*beta_2)*pow((run - ten)*2.*beta_1 + (1. + ten)*beta_2,-1);
	}

  return _SUCCESS_;
}


/**
* Get gravity functions Bs from Gs.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Input/Output: vector of background quantities
* @param pvecback_derivs      Input: vector of derivatives
* @return the error status
*/
int gravity_functions_Bs_from_Gs_smg(
  struct background *pba,
  double a,
  double * pvecback_B,
  double * pvecback,
  struct G_functions_and_derivs *pgf
) {

  /* Define local variables */
  double H = pvecback[pba->index_bg_H];
  double phi = pvecback_B[pba->index_bi_phi_smg];
  double phi_prime = pvecback_B[pba->index_bi_phi_prime_smg];
  double X = 0.5*pow(phi_prime/a,2);

  double G2=pgf->G2;
  double G2_X=pgf->G2_X, G2_XX=pgf->G2_XX, G2_XXX=pgf->G2_XXX;
  double G2_phi=pgf->G2_phi, G2_Xphi=pgf->G2_Xphi, G2_XXphi=pgf->G2_XXphi;
  double G2_phiphi=pgf->G2_phiphi, G2_Xphiphi=pgf->G2_Xphiphi;
  double G3_X=pgf->G3_X, G3_XX=pgf->G3_XX, G3_XXX=pgf->G3_XXX;
  double G3_phi=pgf->G3_phi, G3_Xphi=pgf->G3_Xphi, G3_XXphi=pgf->G3_XXphi;
  double G3_phiphi=pgf->G3_phiphi, G3_Xphiphi=pgf->G3_Xphiphi;
  double G3_phiphiphi=pgf->G3_phiphiphi;
  double G4=pgf->G4, DG4=pgf->DG4;
  double G4_X=pgf->G4_X, G4_XX=pgf->G4_XX, G4_XXX=pgf->G4_XXX, G4_XXXX=pgf->G4_XXXX;
  double G4_phi=pgf->G4_phi, G4_Xphi=pgf->G4_Xphi, G4_XXphi=pgf->G4_XXphi, G4_XXXphi=pgf->G4_XXXphi;
  double G4_phiphi=pgf->G4_phiphi, G4_Xphiphi=pgf->G4_Xphiphi, G4_XXphiphi=pgf->G4_XXphiphi;
  double G4_phiphiphi=pgf->G4_phiphiphi, G4_Xphiphiphi=pgf->G4_Xphiphiphi;
  double G5_X=pgf->G5_X, G5_XX=pgf->G5_XX, G5_XXX=pgf->G5_XXX, G5_XXXX=pgf->G5_XXXX;
  double G5_phi=pgf->G5_phi, G5_Xphi=pgf->G5_Xphi, G5_XXphi=pgf->G5_XXphi, G5_XXXphi=pgf->G5_XXXphi;
  double G5_phiphi=pgf->G5_phiphi, G5_Xphiphi=pgf->G5_Xphiphi, G5_XXphiphi=pgf->G5_XXphiphi;
  double G5_phiphiphi=pgf->G5_phiphiphi, G5_Xphiphiphi=pgf->G5_Xphiphiphi;
  double F4=pgf->F4;
  double F4_X=pgf->F4_X, F4_XX=pgf->F4_XX, F4_XXX=pgf->F4_XXX;
  double F4_phi=pgf->F4_phi, F4_Xphi=pgf->F4_Xphi, F4_XXphi=pgf->F4_XXphi;
  double F4_phiphi=pgf->F4_phiphi, F4_Xphiphi=pgf->F4_Xphiphi;
  double F5=pgf->F5;
  double F5_X=pgf->F5_X, F5_XX=pgf->F5_XX, F5_XXX=pgf->F5_XXX;
  double F5_phi=pgf->F5_phi, F5_Xphi=pgf->F5_Xphi, F5_XXphi=pgf->F5_XXphi;
  double F5_phiphi=pgf->F5_phiphi, F5_Xphiphi=pgf->F5_Xphiphi;

  /* Computing B functions (intermediate step for perturbations) */

  /* B0_smg */
  pvecback[pba->index_bg_B0_smg] =
  (
    + G4_phi - X*(3.*G3_X - 10.*G4_Xphi + 2.*G5_phiphi) + 8.*pow(X,2)*F4_phi
    - 3.*(
      + G4_X - G5_phi + 2.*X*(G4_XX - 2./3.*G5_Xphi + 4.*F4)
      + 4.*pow(X,2)*(F4_X + F5_phi)
    )*H*phi_prime/a
    - X*(3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X)*pow(H,2)
    - (G2_X/2. - G3_phi + G4_phiphi)*phi_prime/a/H
  )/pvecback[pba->index_bg_M2_smg];

  /* B1_smg */
  pvecback[pba->index_bg_B1_smg] =
  6.*(
    + (
      + G4_phi + X*(3.*G3_X - 16.*G4_Xphi + 6.*G5_phiphi)
      + 2.*pow(X,2)*(G3_XX - 20.*F4_phi - 6.*G4_XXphi + 2.*G5_Xphiphi)
      - 16.*pow(X,3)*F4_Xphi
    )
    + X*(
      + 3.*G5_X + 12.*X*(G5_XX - 15.*F5)
      + 4.*pow(X,2)*(G5_XXX - 60.*F5_X) - 48.*pow(X,3)*F5_XX
    )*pow(H,2)
    + (
      + (G2_X - 2.*G3_phi + 4.*G4_phiphi)/2. - X*(G3_Xphi - 2.*G4_Xphiphi)
    )*phi_prime/a/H
    + (
      + G4_X - G5_phi + X*(40.*F4 + 14.*G4_XX + (-13.)*G5_Xphi)
      + 2.*pow(X,2)*(34.*F4_X + 4.*G4_XXX + 36.*F5_phi + (-3.)*G5_XXphi)
      + 8.*pow(X,3)*(2.*F4_XX + 3.*F5_Xphi)
    )*H*phi_prime/a
    - 2.*(
      + X*(
        + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
      )
      + (
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
      )*phi_prime/H/a
    )/a*pvecback[pba->index_bg_H_prime]
    + (
      - 2.*(G4_X - G5_phi)
      - 2.*X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
      - 4.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X)
      - 16.*pow(X,3)*F4_XX
      - (G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi))*phi_prime/a/H
      - (
        + 3.*G5_X - X*(120.*F5 - 7.*G5_XX)
        - 2.*pow(X,2)*(66.*F5_X - G5_XXX) - 24.*pow(X,3)*F5_XX
      )*H*phi_prime/a
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* B2_smg */
  pvecback[pba->index_bg_B2_smg] =
  6.*(
    + (G2_phi/2. - X*(G3_phiphi - 2.*G4_phiphiphi))*pow(H,-2)
    + 3.*G4_phi - X*(2.*G4_Xphi + G5_phiphi)
    + 2.*pow(X,2)*(4.*G4_XXphi - 3.*G5_Xphiphi + 10.*F4_phi)
    + 8.*pow(X,3)*(2.*F4_Xphi + 3.*F5_phiphi)
    + X*(
      + G5_Xphi + 2.*X*(G5_XXphi - 18.*F5_phi) - 24.*pow(X,2)*F5_Xphi
    )*H*phi_prime/a
    + (
      + G4_phiphi + X*(G3_Xphi - 6.*G4_Xphiphi + 2.*G5_phiphiphi)
      - 8.*pow(X,2)*F4_phiphi
    )*phi_prime/a/H
    - 2.*(
      + X*(G5_Xphi - 12.*X*F5_phi)*phi_prime/a
      - (G4_phi - X*(2.*G4_Xphi - G5_phiphi) - 4.*pow(X,2)*F4_phi)/H
    )*pvecback[pba->index_bg_H_prime]/a/H
    + (
      - X*(3.*G5_Xphi + 2.*X*(G5_XXphi - 30.*F5_phi) - 24.*pow(X,2)*F5_Xphi)
      + (G4_phiphi - X*(G3_Xphi - 2.*G4_Xphiphi))*pow(H,-2)
      - 2.*(
        + G4_Xphi - G5_phiphi + X*(2.*G4_XXphi - G5_Xphiphi + 8.*F4_phi)
        + 4.*pow(X,2)*F4_Xphi
      )*phi_prime/a/H
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* B3_smg */
  pvecback[pba->index_bg_B3_smg] =
  4.*(
    + G4_phi - X*(2.*G4_Xphi - G5_phiphi) - 4.*pow(X,2)*F4_phi
    + X*(
      + G5_X + 2.*X*(G5_XX - 18.*F5) - 24.*pow(X,2)*F5_X
    )*pow(H,2)
    + 2.*X*(
      + G4_XX - G5_Xphi + 2.*F4 + 2.*X*(F4_X + 3.*F5_phi)
    )*H*phi_prime/a
    - X*(G5_X - 12.*X*F5)*pvecback[pba->index_bg_H_prime]/a
    - (
      + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 6.*F4) + 4.*pow(X,2)*F4_X
      + (G5_X + X*(G5_XX - 24.*F5) - 12.*pow(X,2)*F5_X)*H*phi_prime/a
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* B4_smg */
  pvecback[pba->index_bg_B4_smg] =
  2.*(
    + G4_phi - X*(2.*G4_Xphi - G5_phiphi) - 4.*pow(X,2)*F4_phi
    + 2.*X*(
      + 2.*F4 + G4_XX - G5_Xphi + 2.*X*(F4_X + 3.*F5_phi)
    )*H*phi_prime/a
    + X*(
      + G5_X + 2.*X*(G5_XX - 18.*F5) - 24.*pow(X,2)*F5_X
    )*pow(H,2)
    - X*(G5_X - 12.*X*F5)*pvecback[pba->index_bg_H_prime]/a
    - (
      + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 6.*F4) + 4.*pow(X,2)*F4_X
      + (
        + G5_X - X*(24.*F5 - G5_XX) - 12.*pow(X,2)*F5_X
      )*H*phi_prime/a
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* B5_smg */
  pvecback[pba->index_bg_B5_smg] =
  (
    - 3.*G4_phi + X*(3.*G3_X - 4.*G4_Xphi - 2.*G5_phiphi)
    - 2.*pow(X,2)*(G3_XX - 6.*G4_XXphi + 2.*G5_Xphiphi - 12.*F4_phi) + 16.*pow(X,3)*F4_Xphi
    + (
      + G2_X - 2.*G3_phi + 2.*X*(G3_Xphi - 2.*G4_Xphiphi)
    )*phi_prime/a/H/2.
    + (
      + 5.*(G4_X - G5_phi) - X*(2.*G4_XX - 5.*G5_Xphi - 8.*F4)
      - 2.*pow(X,2)*( + 4.*G4_XXX - 3.*G5_XXphi + 22.*F4_X + 24.*F5_phi)
      - 8.*pow(X,3)*(2.*F4_XX + 3.*F5_Xphi)
    )*H*phi_prime/a
    + X*(
      + 3.*G5_X - 4.*X*(2.*G5_XX - 15.*F5)
      - 4.*pow(X,2)*(G5_XXX - 48.*F5_X) + 48.*pow(X,3)*F5_XX
    )*pow(H,2)
    + 2.*(
      + X*(3.*G5_X - 2.*X*(30.*F5 - G5_XX) - 24.*pow(X,2)*F5_X)
      + (
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
      )*phi_prime/a/H
    )*pvecback[pba->index_bg_H_prime]/a
    + (
      + 2.*(
        + G4_X - G5_phi + X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
        + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X) + 8.*pow(X,3)*F4_XX
      )
      + (
        + 3.*G5_X + X*(7.*G5_XX - 120.*F5)
        + 2.*pow(X,2)*(G5_XXX - 66.*F5_X) - 24.*pow(X,3)*F5_XX
      )*H*phi_prime/a
      + (
        + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
      )*phi_prime/a/H
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* B6_smg */
  pvecback[pba->index_bg_B6_smg] =
  4.*(
    + G4_phi - X*(2.*G4_Xphi - G5_phiphi)
    + 2.*X*(G4_XX - G5_Xphi)*H*phi_prime/a
    + X*(G5_X + 2.*X*G5_XX)*pow(H,2)
    - X*G5_X*pvecback[pba->index_bg_H_prime]/a
    - (
      + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi)
      + (G5_X + X*G5_XX)*H*phi_prime/a
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* B7_smg */
  pvecback[pba->index_bg_B7_smg] =
  - 2.*(
    + G2_X - 2.*G3_phi
    - X*(G2_XX - 8.*G3_Xphi + 18.*G4_Xphiphi)
    - 2.*pow(X,2)*(G2_XXX - 4.*G3_XXphi + 6.*G4_XXphiphi)
    + (
      + G2_Xphi - 2.*G3_phiphi + 2.*X*(G2_XXphi - G3_Xphiphi)
    )*phi_prime/a/H/2.
    + 6.*(
      + G4_X - G5_phi - X*(G4_XX - 2.*G5_Xphi)
      - 4./3.*pow(X,2)*(9.*G4_XXX - 7.*G5_XXphi + 45.*F4_X + 30.*F5_phi)
      - 4./3.*pow(X,3)*(3.*G4_XXXX - 2.*G5_XXXphi + 39.*F4_XX + 33.*F5_Xphi)
      - 8.*pow(X,4)*(F4_XXX + F5_XXphi)
    )*pow(H,2)
    + (
      + 3.*G5_X - X*(13.*G5_XX - 120.*F5) - 4.*pow(X,2)*(5.*G5_XXX - 159.*F5_X)
      - 4.*pow(X,3)*(G5_XXXX - 96.*F5_XX) + 48.*pow(X,4)*F5_XXX
    )*pow(H,3)*phi_prime/a
    + 3.*(
      + G3_X - 2.*G4_Xphi - G5_phiphi
      - X*(3.*G3_XX - 16.*G4_XXphi + 5.*G5_Xphiphi - 24.*F4_phi)
      - 2.*pow(X,2)*(G3_XXX - 4.*G4_XXXphi + G5_XXphiphi - 18.*F4_Xphi)
      + 8.*pow(X,3)*F4_XXphi
    )*H*phi_prime/a
    + 6.*(
      + G4_X - G5_phi + X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
      + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X) + 8.*pow(X,3)*F4_XX
      - (
        - 3.*G5_X - X*(7.*G5_XX - 120.*F5) - 2.*pow(X,2)*(G5_XXX - 66.*F5_X)
        + 24.*pow(X,3)*F5_XX
      )*H*phi_prime/a/2.
      + (
        + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
      )*phi_prime/a/H/2.
    )*pvecback[pba->index_bg_H_prime]/a
    + 3.*(
      + (G3_X - 3.*G4_Xphi)
      + X*(5.*G3_XX - 12.*G4_XXphi)
      + 2.*pow(X,2)*(G3_XXX - 2.*G4_XXXphi)
      + (
        + G5_X + 3.*X*(3.*G5_XX - 40.*F5) + 4.*pow(X,2)*(2.*G5_XXX - 75.*F5_X)
        + 4./3.*pow(X,3)*(G5_XXXX - 108.*F5_XX) - 16.*pow(X,4)*F5_XXX
      )*pow(H,2)
      + (
        + G2_XX/2. - 2./3.*G3_Xphi + X*(G2_XXX - G3_XXphi)/3.
      )*phi_prime/a/H
      + (
        + 3.*(3.*G4_XX - 2.*G5_Xphi + 8.*F4)
        + X*(16.*G4_XXX - 9.*G5_XXphi + 96.*F4_X)
        + 2.*pow(X,2)*(2.*G4_XXXX - G5_XXXphi + 30.*F4_XX)
        + 8.*pow(X,3)*F4_XXX
      )*H*phi_prime/a
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* B8_smg */
  pvecback[pba->index_bg_B8_smg] =
  - 2.*(
    + G2_X/2. - G3_phi + X*(G3_Xphi - 2.*G4_Xphiphi)
    + (
      + 3.*(G4_X - G5_phi) - X*(2.*G4_XX - 3.*G5_Xphi + 8.*F4)
      - 2.*pow(X,2)*(4.*G4_XXX - 3.*G5_XXphi + 26.*F4_X + 24.*F5_phi)
      - 8.*pow(X,3)*(2.*F4_XX + 3.*F5_Xphi)
    )*pow(H,2)
    + (
      + G3_X - 3.*G4_Xphi - X*(G3_XX - 6.*G4_XXphi + 2.*G5_Xphiphi - 12.*F4_phi)
      + 8.*pow(X,2)*F4_Xphi
    )*H*phi_prime/a
    + (
      + G5_X - 3.*X*(G5_XX - 20.*F5) - 2.*pow(X,2)*(G5_XXX - 54.*F5_X)
      + 24.*pow(X,3)*F5_XX
    )*pow(H,3)*phi_prime/a
    + 2.*(
      + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 10.*F4) + 4.*pow(X,2)*F4_X
      + (
        + G5_X + X*(G5_XX - 36.*F5) - 12.*pow(X,2)*F5_X
      )*H*phi_prime/a
    )*pvecback[pba->index_bg_H_prime]/a
    + (
      + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
      + (
        + G5_X + 5.*X*(G5_XX - 24.*F5) + 2.*pow(X,2)*(G5_XXX - 66.*F5_X)
        - 24.*pow(X,3)*F5_XX
      )*pow(H,2)
      + 2.*(
        + 3.*G4_XX - 2.*G5_Xphi + 12.*F4 + X*(2.*G4_XXX - G5_XXphi + 18.*F4_X)
        + 4.*pow(X,2)*F4_XX
      )*H*phi_prime/a
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* B9_smg */
  pvecback[pba->index_bg_B9_smg] =
  2.*(
    + 6.*G4_phiphi
    - 3.*X*(G3_Xphi - G5_phiphiphi)
    + 6.*pow(X,2)*(G3_XXphi - 4.*G4_XXphiphi + G5_Xphiphiphi - 6.*F4_phiphi)
    - 24.*pow(X,3)*F4_Xphiphi
    + (
      + G2_phiphi - (G2_Xphiphi - G3_phiphiphi)*2.*X
    )*pow(H,-2)/2.
    - (
      + G2_Xphi - 2.*G3_phiphi - X*(G2_XXphi - 4.*G3_Xphiphi + 6.*G4_Xphiphiphi)
    )*phi_prime/a/H
    - 2.*(
      + 3.*(G4_Xphi - G5_phiphi) - X*(3.*G4_XXphi - 4.*G5_Xphiphi)
      - 2.*pow(X,2)*(3.*G4_XXXphi - 2.*G5_XXphiphi + 18.*F4_Xphi + 12.*F5_phiphi)
      - 12.*pow(X,3)*(F4_XXphi + F5_Xphiphi)
    )*H*phi_prime/a
    - X*(
      + 3.*G5_Xphi - 4.*X*(2.*G5_XXphi - 15.*F5_phi)
      - 4.*pow(X,2)*(G5_XXXphi - 48.*F5_Xphi) + 48.*pow(X,3)*F5_XXphi
    )*pow(H,2)
    - 3.*(
      + X*(3.*G5_Xphi + 2.*X*(G5_XXphi - 30.*F5_phi) - 24.*pow(X,2)*F5_Xphi)
      - (G4_phiphi - X*(G3_Xphi - 2.*G4_Xphiphi))*pow(H,-2)
      + 2.*(
        + G4_Xphi - G5_phiphi + X*(8.*F4_phi + 2.*G4_XXphi - G5_Xphiphi)
        + 4.*pow(X,2)*F4_Xphi
      )*phi_prime/a/H
    )*pvecback[pba->index_bg_H_prime]/a
    - (
      + 3.*(G4_Xphi - G5_phiphi)
      + 3.*X*(8.*G4_XXphi - 5.*G5_Xphiphi + 24.*F4_phi)
      + 6.*pow(X,2)*(2.*G4_XXXphi - G5_XXphiphi + 18.*F4_Xphi)
      + 24.*pow(X,3)*F4_XXphi
      + (
        + G2_Xphi - 2.*G3_phiphi + 2.*X*(G2_XXphi - G3_Xphiphi)
      )*pow(H,-2)/2.
      + 3.*(
        + G3_Xphi - 3.*G4_Xphiphi + X*(G3_XXphi - 2.*G4_XXphiphi)
      )*phi_prime/a/H
      + (
        + 3.*G5_Xphi + X*(7.*G5_XXphi - 120.*F5_phi)
        - 2.*pow(X,2)*(66.*F5_Xphi - G5_XXXphi) - 24.*pow(X,3)*F5_XXphi
      )*H*phi_prime/a
    )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
  )/pvecback[pba->index_bg_M2_smg];

  /* B10_smg */
  pvecback[pba->index_bg_B10_smg] =
  (
    + 3.*(
      + G4_phi - X*(3.*G3_X - 8.*G4_Xphi) - 2.*pow(X,2)*(G3_XX - 2.*G4_XXphi)
    )
    - (
      + G2_X - 2.*G3_phi + 2.*X*(G2_XX - G3_Xphi)
    )*phi_prime/a/H/2.

    - X*(
      + 15.*G5_X + 20.*X*(G5_XX - 21.*F5) + 4.*pow(X,2)*(G5_XXX - 84.*F5_X)
      - 48.*pow(X,3)*F5_XX
    )*pow(H,2)
    - 3.*(
      + 3.*(G4_X - G5_phi) + X*(12.*G4_XX - 7.*G5_Xphi + 40.*F4)
      + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 22.*F4_X) + 8.*pow(X,3)*F4_XX
    )*H*phi_prime/a
  )/pvecback[pba->index_bg_M2_smg];

  /* B11_smg */
  pvecback[pba->index_bg_B11_smg] =
  (
    + G4_phi
    - X*(G3_X - 2.*G4_Xphi)
    - X*(
      + 3.*G5_X + 2.*X*(G5_XX - 42.*F5) - 24.*pow(X,2)*F5_X
    )*pow(H,2)
    - 2.*(
      + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 10.*F4) + 4.*pow(X,2)*F4_X
    )*H*phi_prime/a
  )/pvecback[pba->index_bg_M2_smg];

  /* B12_smg */
  pvecback[pba->index_bg_B12_smg] =
  (
    + 3.*G4_phi
    - 3.*X*(4.*G4_Xphi - 3.*G5_phiphi)
    - 6.*pow(X,2)*(10.*F4_phi + 2.*G4_XXphi - G5_Xphiphi)
    - 24.*pow(X,3)*F4_Xphi
    + (
      + G2_phi/2. - X*(G2_Xphi - G3_phiphi)
    )*pow(H,-2)
    - X*(
      + 5.*G5_Xphi + 2.*X*(G5_XXphi - 42.*F5_phi) - 24.*pow(X,2)*F5_Xphi
    )*H*phi_prime/a
    + 3.*(
      + G4_phiphi - X*(G3_Xphi - 2.*G4_Xphiphi)
    )*phi_prime/a/H
  )/pvecback[pba->index_bg_M2_smg];

  return _SUCCESS_;
}


/**
* Get gravity functions Cs from Bs.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Input/Output: vector of background quantities
* @return the error status
*/
int gravity_functions_Cs_from_Bs_smg(
                                     struct background *pba,
                                     double * pvecback,
                                     double * pvecback_derivs
                                     ) {

  // basic background quantities
  double a = pvecback[pba->index_bg_a];
  double H = pvecback[pba->index_bg_H];
  double H_p = pvecback[pba->index_bg_H_prime];
  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];
  double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
  double p_smg = pvecback[pba->index_bg_p_smg];
  double rho_smg = pvecback[pba->index_bg_rho_smg];
  // factor to convert lna derivatives to tau derivatives
  double factor = a*H;

  // alphas
  double M2 = pvecback[pba->index_bg_M2_smg];
  double DelM2 = pvecback[pba->index_bg_delta_M2_smg];
  double kin = pvecback[pba->index_bg_kineticity_smg];
  double bra = pvecback[pba->index_bg_braiding_smg];
  double run = pvecback[pba->index_bg_M2_running_smg];
  double ten = pvecback[pba->index_bg_tensor_excess_smg];
  double beh = pvecback[pba->index_bg_beyond_horndeski_smg];
  double dM2 = pvecback[pba->index_bg_delta_M2_smg];

  // need to update the time derivatives of the interesting functions
  double kin_p = factor*pvecback_derivs[pba->index_bg_kineticity_smg];
  double bra_p = factor*pvecback_derivs[pba->index_bg_braiding_smg];
  double run_p = factor*pvecback_derivs[pba->index_bg_M2_running_smg];
  double ten_p = factor*pvecback_derivs[pba->index_bg_tensor_excess_smg];
  double beh_p = factor*pvecback_derivs[pba->index_bg_beyond_horndeski_smg];
  double p_tot_p = factor*pvecback_derivs[pba->index_bg_p_tot_wo_smg];
  double p_smg_p = factor*pvecback_derivs[pba->index_bg_p_smg];

  // alphas over scalar field
  double kin_ss = pvecback[pba->index_bg_kineticity_over_phiphi_smg];
  double bra_s = pvecback[pba->index_bg_braiding_over_phi_smg];
  double beh_s = pvecback[pba->index_bg_beyond_horndeski_over_phi_smg];

  // Bs
  double B0 = pvecback[pba->index_bg_B0_smg];
  double B1 = pvecback[pba->index_bg_B1_smg];
  double B2 = pvecback[pba->index_bg_B2_smg];
  double B3 = pvecback[pba->index_bg_B3_smg];
  double B4 = pvecback[pba->index_bg_B4_smg];
  double B5 = pvecback[pba->index_bg_B5_smg];
  double B6 = pvecback[pba->index_bg_B6_smg];
  double B7 = pvecback[pba->index_bg_B7_smg];
  double B8 = pvecback[pba->index_bg_B8_smg];
  double B9 = pvecback[pba->index_bg_B9_smg];
  double B10 = pvecback[pba->index_bg_B10_smg];
  double B11 = pvecback[pba->index_bg_B11_smg];
  double B12 = pvecback[pba->index_bg_B12_smg];

  // kinetic term D over phiphi
  pvecback[pba->index_bg_kinetic_D_over_phiphi_smg] =
  + kin_ss + 3./2.*pow(bra_s,2);

  // C0
  pvecback[pba->index_bg_C0_smg] = B0;

  // C1
  pvecback[pba->index_bg_C1_smg] = kin_ss*(1. + ten) - 3./2.*bra_s*B6;

  // C2
  pvecback[pba->index_bg_C2_smg] = - kin_ss*(2. + run) - 3.*bra_s*B5;

  // C3
  pvecback[pba->index_bg_C3_smg] = bra_s*beh_s;

  // C4
  pvecback[pba->index_bg_C4_smg] = kin_ss*B1 - 3.*bra_s*B7;

  // C5
  pvecback[pba->index_bg_C5_smg] = - kin_ss*beh_s;

  // C6
  pvecback[pba->index_bg_C6_smg] = kin_ss*B2 - 3.*bra_s*B9;

  // C7
  pvecback[pba->index_bg_C7_smg] = kin_ss*B3 - 3.*bra_s*B8;

  // C8
  pvecback[pba->index_bg_C8_smg] = B4;

  // C9
  pvecback[pba->index_bg_C9_smg] = - bra_s - bra_s*run/2. + B5;

  // C10
  pvecback[pba->index_bg_C10_smg] = bra_s*(1. + ten) + B6;

  // C11
  pvecback[pba->index_bg_C11_smg] = bra_s*B1/2. + B7;

  // C12
  pvecback[pba->index_bg_C12_smg] = bra_s*B2/2. + B9;

  // C13
  pvecback[pba->index_bg_C13_smg] = bra_s*B3/2. + B8;

  // C14
  pvecback[pba->index_bg_C14_smg] = B10;

  // C15
  pvecback[pba->index_bg_C15_smg] = B11;

  // C16
  pvecback[pba->index_bg_C16_smg] = B12;

  return _SUCCESS_;
}
