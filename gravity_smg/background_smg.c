/** @file background_smg.c Documented background_smg module
 *
 * Emilio Bellini, Ignacy Sawicki, Miguel Zumalacarregui, 2024
 *
 * Additional functions for the background module.
 * It contains all the hi_class related functions (_smg)
 * that are used by background.c. In this way the main hi_class
 * modifications are stored here and the standard Class modules
 * remain cleaner.
 *
 * The following nomenclature has been adopted:
 *
 * -# all the functions end with "_smg" to make them easily
 *    recognizable
 * -# all the functions starting with "background" are
 *    directly called by background.c or the classy wrapper
 * -# all the functions that do not start with "background"
 *    are only used internally in background_smg.c
 */

#include "background_smg.h"


/**
* Fill the modified gravity part of background_functions.
* First all the Horndeski functions G_i(X,phi) are computed.
* A loop is used to allow different implementations without erasing previous ones
* Note that in CLASS units a canonical field has G2 = X ~ [Mpc]^-2
* This implies that phi_prime ~ [Mpc]^-1
* and phi is dimensionless (we can think of phi as given in units of the Planck mass
* - A default module to numerically compute the derivatives when no analytic functions are given should be added.
* Numerical derivatives may further serve as a consistency check.
*
* @param pba           Input: pointer to background structure
* @param a             Input: scale factor
* @param pvecback_B    Input: vector containing all {B} quantities
* @param pvecback      Output: vector of background quantities (assumed to be already allocated)
* @param ptr_rho_tot   Output: pointer to the total (without smg) energy density
* @param ptr_p_tot     Output: pointer to the total (without smg) pressure
* @param ptr_rho_de    Output: pointer to the dark energy density (smg)
* @return the error status
*/
int background_gravity_functions_smg(
				 														 struct background *pba,
				 													 	 double a,
				 													   double * pvecback_B,
				 													   double * pvecback,
         									 					 double * ptr_rho_tot,
				 													   double * ptr_p_tot,
				 													   double * ptr_rho_de
				 												 		 ){


	/** - save rho_tot and p_tot without smg (it is important that smg is the last thing to be computed) */
	pvecback[pba->index_bg_rho_tot_wo_smg] = *ptr_rho_tot;
	pvecback[pba->index_bg_p_tot_wo_smg] = *ptr_p_tot;

	// scalar field + curvature not yet implemented
	class_test(pba->K !=0 ,
	     pba->error_message,
	     "has_smg with curvature K = %e not yet implemented",pba->K);

  if (pba->field_evolution_smg == _TRUE_) {

		/* Define local variables */
		double x,f,df;
    int n, n_max=100;
		double H;

		/* Get phi and phi_prime from the integrator */
		pvecback[pba->index_bg_phi_smg] = pvecback_B[pba->index_bi_phi_smg];
    pvecback[pba->index_bg_phi_prime_smg] = pvecback_B[pba->index_bi_phi_prime_smg];

		/* declare G functions and set defaults */
		struct G_functions_and_derivs gf = DEFAULT_G_FUNCTIONS_AND_DERIVS;

		/* update G functions and derivatives */
 	  class_call(gravity_models_get_Gs_smg(pba, a, pvecback_B, &gf),
 	    pba->error_message,
 	    pba->error_message
 	  );

		/* get Es functions. Notation:
		* E0 + E1 H + E3 H^3 = E2 H^2
		*/
		class_call(gravity_functions_Es_from_Gs_smg(pba, a, pvecback_B, pvecback, &gf),
			pba->error_message,
			pba->error_message
		);

		double E0 = pvecback[pba->index_bg_E0_smg];
		double E1 = pvecback[pba->index_bg_E1_smg];
		double E2 = pvecback[pba->index_bg_E2_smg];
		double E3 = pvecback[pba->index_bg_E3_smg];

		class_test(E3*pow(E0,1./2.) > 1e-10 && pba->initial_conditions_set_smg == _FALSE_,
							 pba->error_message,
							 "E3=%e is large in Friedmann constraint when setting ICs ",  E3);

		/* get Hubble, either solving the cubic Friedmann
		* equation or from the integrator.
		*/
		if ((pba->hubble_evolution == _FALSE_) || (pba->initial_conditions_set_smg == _FALSE_)) {
      /* Use Newton's method to solve the cubic Friedmann equation */
      x = sqrt(E0);
      f = E3*x*x*x -E2*x*x + E1*x + E0;
      n = 0;
      while (fabs(f/E0)> 1e-8 && n < 100){
				f = E3*x*x*x - E2*x*x + E1*x + E0;
				df = 3.*E3*x*x - 2.*E2*x + E1;
				x -= f/df;
				n++;
      }
      H=x;
      if ((pba->background_verbose > 5) && (pba->initial_conditions_set_smg == _FALSE_ ))
				printf(" Initial H = %e, sqrt(rho) = %e, ratio = %e, n=%i \n", H, sqrt(pvecback[pba->index_bg_rho_tot_wo_smg]),sqrt(pvecback[pba->index_bg_rho_tot_wo_smg])/H,n);

 	    class_test(isnan(H),
 		       			 pba->error_message,
 	               "H=%e is not a number at a = %e. phi = %e, phi_prime = %e, E0=%e, E1=%e, E2=%e, E3=%e ",
 		       			 H,a,pvecback[pba->index_bg_phi_smg],pvecback[pba->index_bg_phi_prime_smg],E0,E1,E2,E3);

 		}
		else {
			H = exp(pvecback_B[pba->index_bi_logH]);
		}
		pvecback[pba->index_bg_H] = H;

		/* get Ps and Rs functions. Notation:
		* P0 phi'' + P1 H' + P2 = 0
		* R0 phi'' + R1 H' + R2 = 0
		*/
 	  class_call(gravity_functions_Ps_and_Rs_from_Gs_smg(pba, a, pvecback_B, pvecback, &gf),
 	    pba->error_message,
 	    pba->error_message
 	  );

		double P0 = pvecback[pba->index_bg_P0_smg];
    double P1 = pvecback[pba->index_bg_P1_smg];
    double P2 = pvecback[pba->index_bg_P2_smg];
    double R0 = pvecback[pba->index_bg_R0_smg];
    double R1 = pvecback[pba->index_bg_R1_smg];
    double R2 = pvecback[pba->index_bg_R2_smg];

    class_test_except((P1*R0 - P0*R1) == 0 ,
	       pba->error_message,
         free(pvecback);free(pvecback_B),
               "scalar field mixing with metric has degenerate denominator at a = %e, phi = %e, phi_prime = %e \n with P1 = %e, R0 =%e, P0=%e, R1=%e, \n H=%e, E0=%e, E1=%e, E2=%e, E3=%e",
	       a,pvecback[pba->index_bg_phi_smg],pvecback[pba->index_bg_phi_prime_smg], P1, R0, P0, R1,
	       pvecback[pba->index_bg_H],E0,E1,E2,E3);

    /* Friedmann space-space equation with friction added */
    pvecback[pba->index_bg_H_prime] = (R0*P2 - P0*R2)/(P0*R1 - P1*R0);
    /* choose sign for friction depending on the derivative */
    if ((2.*E2*H - E1 - 3*E3*H*H)>=0)
      pvecback[pba->index_bg_H_prime] += - a*pba->hubble_friction*(E2*H*H - E0 - E1*H - E3*H*H*H);
    else{
      pvecback[pba->index_bg_H_prime] += a*pba->hubble_friction*(E2*H*H - E0 - E1*H - E3*H*H*H);
    }

    /* Field equation */
    pvecback[pba->index_bg_phi_prime_prime_smg] = (P1*R2 - R1*P2)/(P0*R1 - P1*R0);

		/* get alphas, background density and pressure, shift and current. */
		class_call(gravity_functions_building_blocks_from_Gs_smg(pba, a, pvecback_B, pvecback, &gf),
 	    pba->error_message,
 	    pba->error_message
 	  );

		/* get Bs, intermediate step for perturbations. */
		class_call(gravity_functions_Bs_from_Gs_smg(pba, a, pvecback_B, pvecback, &gf),
 	    pba->error_message,
 	    pba->error_message
 	  );

	}
	// end of if pba->field_evolution_smg
  else{

		double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
	  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];

		/* get background parametrizations. */
		class_call(gravity_models_get_back_par_smg(pba, a, pvecback, pvecback_B),
 	    pba->error_message,
 	    pba->error_message
 	  );

		/* add _smg to rho_tot */
		rho_tot += pvecback[pba->index_bg_rho_smg];
    p_tot += pvecback[pba->index_bg_p_smg];

    pvecback[pba->index_bg_H] = sqrt(rho_tot-pba->K/a/a);
    /** - compute derivative of H with respect to conformal time */
    pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;

    //add friction term
    if (pba->hubble_evolution == _TRUE_ && pba->initial_conditions_set_smg == _TRUE_){
      pvecback[pba->index_bg_H] = exp(pvecback_B[pba->index_bi_logH]);
      /** - compute derivative of H with respect to conformal time */
      pvecback[pba->index_bg_H_prime] += - a*pba->hubble_friction*(pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H] - rho_tot - pba->K/a/a);
    }

    // Compute time derivative of rho_smg
    if (pba->rho_evolution_smg == _TRUE_){
        pvecback[pba->index_bg_rho_prime_smg] = -3.*a*pvecback[pba->index_bg_H]*(1.+pvecback[pba->index_bg_w_smg])*pvecback[pba->index_bg_rho_smg];
    }


		/* initialize the values to the defaults */
    pvecback[pba->index_bg_kineticity_smg] = 0;
    pvecback[pba->index_bg_braiding_smg] = 0.;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_M2_smg] = 1.;
    pvecback[pba->index_bg_delta_M2_smg] = 0.;
    pvecback[pba->index_bg_M2_running_smg] = 0.;

		/* get background parametrizations. */
		class_call(gravity_models_get_alphas_par_smg(pba, a, pvecback, pvecback_B),
 	    pba->error_message,
 	    pba->error_message
 	  );

	}
	//end of parameterized mode

  // add a value to the kineticity to avoid problems with perturbations in certain models.
  // NOTE: this needs to be done here to avoid interfering with the equations
  pvecback[pba->index_bg_kineticity_smg] += pba->kineticity_safe_smg;

  //Derivatives of the BS functions and others. Set to zero here and computed numerically once the background is integrated (needed so that debuggers don't complain).

  pvecback[pba->index_bg_kineticity_prime_smg] = 0.;
  pvecback[pba->index_bg_braiding_prime_smg] = 0.;
  pvecback[pba->index_bg_M2_running_prime_smg] = 0.;
  pvecback[pba->index_bg_tensor_excess_prime_smg] = 0.;
  pvecback[pba->index_bg_beyond_horndeski_prime_smg] = 0.;
  pvecback[pba->index_bg_H_prime_prime] = 0.;
  pvecback[pba->index_bg_p_tot_wo_prime_smg] = 0.;
  pvecback[pba->index_bg_p_prime_smg] = 0.;
  pvecback[pba->index_bg_cs2_smg] = 0.;
  pvecback[pba->index_bg_kinetic_D_smg] = 0.;
  pvecback[pba->index_bg_kinetic_D_prime_smg] = 0.;
	if (pba->field_evolution_smg == _TRUE_) {
	  pvecback[pba->index_bg_kinetic_D_over_phiphi_smg] = 0.;
	  pvecback[pba->index_bg_kinetic_D_over_phiphi_prime_smg] = 0.;
	}
  pvecback[pba->index_bg_cs2num_smg] = 0.;
  pvecback[pba->index_bg_cs2num_prime_smg] = 0.;
  pvecback[pba->index_bg_A0_smg] = 0.;
  pvecback[pba->index_bg_A1_smg] = 0.;
  pvecback[pba->index_bg_A2_smg] = 0.;
  pvecback[pba->index_bg_A3_smg] = 0.;
  pvecback[pba->index_bg_A4_smg] = 0.;
  pvecback[pba->index_bg_A5_smg] = 0.;
  pvecback[pba->index_bg_A6_smg] = 0.;
  pvecback[pba->index_bg_A7_smg] = 0.;
  pvecback[pba->index_bg_A8_smg] = 0.;
  pvecback[pba->index_bg_A9_smg] = 0.;
  pvecback[pba->index_bg_A10_smg] = 0.;
  pvecback[pba->index_bg_A11_smg] = 0.;
  pvecback[pba->index_bg_A12_smg] = 0.;
  pvecback[pba->index_bg_A13_smg] = 0.;
  pvecback[pba->index_bg_A14_smg] = 0.;
  pvecback[pba->index_bg_A15_smg] = 0.;
  pvecback[pba->index_bg_A16_smg] = 0.;
  pvecback[pba->index_bg_A9_prime_smg] = 0.;
  pvecback[pba->index_bg_A10_prime_smg] = 0.;
  pvecback[pba->index_bg_A12_prime_smg] = 0.;
  pvecback[pba->index_bg_A13_prime_smg] = 0.;
	if (pba->field_evolution_smg == _TRUE_) {
		pvecback[pba->index_bg_C0_smg] = 0.;
	  pvecback[pba->index_bg_C1_smg] = 0.;
	  pvecback[pba->index_bg_C2_smg] = 0.;
	  pvecback[pba->index_bg_C3_smg] = 0.;
	  pvecback[pba->index_bg_C4_smg] = 0.;
	  pvecback[pba->index_bg_C5_smg] = 0.;
	  pvecback[pba->index_bg_C6_smg] = 0.;
	  pvecback[pba->index_bg_C7_smg] = 0.;
	  pvecback[pba->index_bg_C8_smg] = 0.;
	  pvecback[pba->index_bg_C9_smg] = 0.;
	  pvecback[pba->index_bg_C10_smg] = 0.;
	  pvecback[pba->index_bg_C11_smg] = 0.;
	  pvecback[pba->index_bg_C12_smg] = 0.;
	  pvecback[pba->index_bg_C13_smg] = 0.;
	  pvecback[pba->index_bg_C14_smg] = 0.;
	  pvecback[pba->index_bg_C15_smg] = 0.;
	  pvecback[pba->index_bg_C16_smg] = 0.;
	  pvecback[pba->index_bg_C9_prime_smg] = 0.;
	  pvecback[pba->index_bg_C10_prime_smg] = 0.;
	  pvecback[pba->index_bg_C12_prime_smg] = 0.;
	  pvecback[pba->index_bg_C13_prime_smg] = 0.;
	}
  pvecback[pba->index_bg_lambda_1_smg] = 0.;
  pvecback[pba->index_bg_lambda_2_smg] = 0.;
  pvecback[pba->index_bg_lambda_3_smg] = 0.;
  pvecback[pba->index_bg_lambda_4_smg] = 0.;
  pvecback[pba->index_bg_lambda_5_smg] = 0.;
  pvecback[pba->index_bg_lambda_6_smg] = 0.;
  pvecback[pba->index_bg_lambda_7_smg] = 0.;
  pvecback[pba->index_bg_lambda_8_smg] = 0.;
  pvecback[pba->index_bg_lambda_9_smg] = 0.;
  pvecback[pba->index_bg_lambda_10_smg] = 0.;
  pvecback[pba->index_bg_lambda_11_smg] = 0.;
  pvecback[pba->index_bg_lambda_2_prime_smg] = 0.;
  pvecback[pba->index_bg_lambda_8_prime_smg] = 0.;
  pvecback[pba->index_bg_lambda_9_prime_smg] = 0.;
  pvecback[pba->index_bg_lambda_11_prime_smg] = 0.;
  pvecback[pba->index_bg_G_eff_smg] = 0.;
  pvecback[pba->index_bg_slip_eff_smg] = 0.;
  if (pba->field_evolution_smg == _FALSE_ && pba->M2_evolution_smg == _FALSE_){
    pvecback[pba->index_bg_M2_running_smg] = 0.;
  }

    /* Check required conditions for the gravity_models. */
  if ( (pba->skip_stability_tests_smg == _FALSE_) && (pba->parameters_tuned_smg == _TRUE_) && (pba->Omega_smg_debug == 0) ) {
		if (pba->is_quintessence_smg == _TRUE_){
		    /*  Check that w is not lower than w < -1 for quintessence */
		   class_test( (pvecback[pba->index_bg_p_smg]/pvecback[pba->index_bg_rho_smg] < -(1 + pba->quintessence_w_safe_smg)),
		       pba->error_message,
		       "Dark Energy equation of state at a = %e is w = %e, lower than w = -(1 + th), with threshold, th = %e. Quintessence models can only have  w > -1. Aborting.\n", a, pvecback[pba->index_bg_p_smg]/pvecback[pba->index_bg_rho_smg], pba->quintessence_w_safe_smg);

		   class_test(( (pvecback[pba->index_bg_rho_smg] / (pba->Omega0_smg*pow(pba->H0,2))  + 1e-3)  < 1),
		           pba->error_message,
		           "Dark Energy density, rho_smg = %e, at a = %e is lower than its given present value, Omega0_smg*H0^2= %e, since quintessence models have  w > -1, it will never be reached. Aborting.\n", pvecback[pba->index_bg_rho_smg], a, pba->Omega0_smg*pow(pba->H0,2) );
		}
  }

	*ptr_rho_tot += pvecback[pba->index_bg_rho_smg];
	*ptr_p_tot += pvecback[pba->index_bg_p_smg];
	//divide relativistic & nonrelativistic (not very meaningful for oscillatory models)

	//TODO: need to define meaningfully -> separate early universe (IC, BBN...) from late (Halofit...)
	//BUG: causes problem with halofit!, if not, causes bug with Brans-Dicke
	*ptr_rho_de += pvecback[pba->index_bg_rho_smg];

	/** - compute w_smg */
	if (pba->rho_evolution_smg == _FALSE_) {
		pvecback[pba->index_bg_w_smg] = pvecback[pba->index_bg_p_smg] / pvecback[pba->index_bg_rho_smg];
	}

  return _SUCCESS_;

}


/**
 * Free all memory space allocated by hi_class in background.
 *
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */
int background_free_smg(
			  								struct background *pba
											  ) {

	if (pba->parameters_smg != NULL)
	  free(pba->parameters_smg);
	//dealocate parameters_2_smg only for parameterizations
	if (pba->field_evolution_smg == _FALSE_ && pba->parameters_2_smg != NULL)
	  free(pba->parameters_2_smg);

	return _SUCCESS_;
}


/**
 * Define hi_class bg indices.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */
int background_define_indices_bg_smg(
				 														 struct background *pba,
				 													   int * index_bg
			 													 		 ) {

	class_define_index(pba->index_bg_phi_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_phi_prime_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_phi_prime_prime_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_rho_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_rho_prime_smg,pba->rho_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_w_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_current_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_shift_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_M2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_delta_M2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_kineticity_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_tensor_excess_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_M2_running_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_kineticity_over_phiphi_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_over_phi_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_over_phi_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_over_phi_prime_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_over_phi_prime_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_kineticity_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_M2_running_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_tensor_excess_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_prime_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_cs2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_cs2num_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_cs2num_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_kinetic_D_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_kinetic_D_prime_smg,_TRUE_,*index_bg,1);
	if (pba->field_evolution_smg == _TRUE_) {
		class_define_index(pba->index_bg_kinetic_D_over_phiphi_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_kinetic_D_over_phiphi_prime_smg,_TRUE_,*index_bg,1);
	}
	class_define_index(pba->index_bg_A0_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A1_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A3_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A4_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A5_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A6_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A7_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A8_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A9_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A10_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A11_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A12_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A13_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A14_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A15_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A16_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A9_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A10_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A12_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A13_prime_smg,_TRUE_,*index_bg,1);

	if (pba->field_evolution_smg == _TRUE_) {
		class_define_index(pba->index_bg_B0_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B1_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B2_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B3_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B4_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B5_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B6_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B7_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B8_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B9_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B10_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B11_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B12_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C0_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C1_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C2_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C3_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C4_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C5_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C6_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C7_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C8_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C9_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C10_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C11_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C12_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C13_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C14_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C15_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C16_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C9_prime_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C10_prime_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C12_prime_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C13_prime_smg,_TRUE_,*index_bg,1);
	}

	class_define_index(pba->index_bg_lambda_1_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_3_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_4_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_5_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_6_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_7_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_8_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_9_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_10_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_11_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_2_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_8_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_9_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_11_prime_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_E0_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_E1_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_E2_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_E3_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_P0_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_P1_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_P2_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_R0_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_R1_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_R2_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_rho_tot_wo_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_tot_wo_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_H_prime_prime,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_tot_wo_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_prime_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_G_eff_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_slip_eff_smg,_TRUE_,*index_bg,1);

  return _SUCCESS_;
}


/**
 * Define hi_class bi indices.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */
int background_define_indices_bi_smg(
				 														 struct background *pba,
				 													   int * index_bi
			 													     ) {

	/* -> scalar field and its derivative wrt conformal time (only if needs to evolve the field)
	* plus other parameters that might be integrated in certain parameterizations
	*/
	class_define_index(pba->index_bi_phi_smg, pba->field_evolution_smg,*index_bi,1);
	class_define_index(pba->index_bi_phi_prime_smg, pba->field_evolution_smg,*index_bi,1);

	//if model needs to integrate M2 from alpha_M, declare an index
	class_define_index(pba->index_bi_delta_M2_smg, pba->M2_evolution_smg,*index_bi,1);

	/* index for the smg energy density */
	class_define_index(pba->index_bi_rho_smg, pba->rho_evolution_smg,*index_bi,1);

  return _SUCCESS_;
}


/**
* In background_solve, after quantities {B} and {C} have been
* inegrated, this routine computes derived quantities, such as
* numerical derivatives of the alphas and functions that depend
* on combinations of alphas and derivatives. This is also the
* place where stability tests are done (in stability_tests_smg).
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
* @return the error status
*/
int background_solve_smg(
												 struct background *pba,
                         double * pvecback,
                         double * pvecback_integration
											   ) {

	/* define local variables */
	double a;

	/* necessary for calling array_interpolate(), but never used */
	int last_index;
	int i;

	/** - second loop over lines, overwrite derivatives that can't be analytically computed from background_functions
	 * Fill the derivatives of the Bellini-Sawicki functions in pvecback
	 * This is done just by overwriting the pvecback entries corresponding to the relevant indice
	 */

	double * pvecback_derivs;
	class_alloc(pvecback_derivs,pba->bg_size*sizeof(double),pba->error_message);

 	for (i=0; i < pba->bt_size; i++) {

		// write the derivatives in the structure
		class_call(array_derivate_spline(pba->loga_table, // x_array
				     pba->bt_size, // int n_lines
				     pba->background_table, // array
				     pba->d2background_dloga2_table, // double * array_splined
				     pba->bg_size, // n_columns
				     pba->loga_table[i], // double x -> tau
				     &last_index, // int* last_index // this is something for the interpolation to talk to each other when using a loop
				     pvecback_derivs, // double * result
				     pba->bg_size, //result_size, from 1 to n_columns
				     pba->error_message),
		pba->error_message,
		pba->error_message);

		//Need to update pvecback
		class_call(background_at_tau(pba,
				 pba->tau_table[i],
				 long_info,
				 inter_normal,
				 &last_index, //should be no problem to use the same one as for the derivatives
				 pvecback),
		pba->error_message,
		pba->error_message);

		a = pvecback[pba->index_bg_a];

 	  /* - indices for scalar field (modified gravity) */
 	  class_call(derivatives_alphas_smg(pba, pvecback, pvecback_derivs, i),
 	    pba->error_message,
 	    pba->error_message
 	  );


		class_call(gravity_functions_As_from_alphas_smg(pba, pvecback, pvecback_derivs),
							 pba->error_message,
							 pba->error_message);

		copy_to_background_table_smg(pba, i, pba->index_bg_kinetic_D_smg, pvecback[pba->index_bg_kinetic_D_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A0_smg, pvecback[pba->index_bg_A0_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A1_smg, pvecback[pba->index_bg_A1_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A2_smg, pvecback[pba->index_bg_A2_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A3_smg, pvecback[pba->index_bg_A3_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A4_smg, pvecback[pba->index_bg_A4_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A5_smg, pvecback[pba->index_bg_A5_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A6_smg, pvecback[pba->index_bg_A6_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A7_smg, pvecback[pba->index_bg_A7_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A8_smg, pvecback[pba->index_bg_A8_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A9_smg, pvecback[pba->index_bg_A9_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A10_smg, pvecback[pba->index_bg_A10_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A11_smg, pvecback[pba->index_bg_A11_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A12_smg, pvecback[pba->index_bg_A12_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A13_smg, pvecback[pba->index_bg_A13_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A14_smg, pvecback[pba->index_bg_A14_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A15_smg, pvecback[pba->index_bg_A15_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_A16_smg, pvecback[pba->index_bg_A16_smg]);

		if (pba->field_evolution_smg == _TRUE_) {

			class_call(gravity_functions_Cs_from_Bs_smg(pba, pvecback, pvecback_derivs),
								 pba->error_message,
								 pba->error_message);

			copy_to_background_table_smg(pba, i, pba->index_bg_kinetic_D_over_phiphi_smg, pvecback[pba->index_bg_kinetic_D_over_phiphi_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C0_smg, pvecback[pba->index_bg_C0_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C1_smg, pvecback[pba->index_bg_C1_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C2_smg, pvecback[pba->index_bg_C2_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C3_smg, pvecback[pba->index_bg_C3_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C4_smg, pvecback[pba->index_bg_C4_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C5_smg, pvecback[pba->index_bg_C5_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C6_smg, pvecback[pba->index_bg_C6_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C7_smg, pvecback[pba->index_bg_C7_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C8_smg, pvecback[pba->index_bg_C8_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C9_smg, pvecback[pba->index_bg_C9_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C10_smg, pvecback[pba->index_bg_C10_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C11_smg, pvecback[pba->index_bg_C11_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C12_smg, pvecback[pba->index_bg_C12_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C13_smg, pvecback[pba->index_bg_C13_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C14_smg, pvecback[pba->index_bg_C14_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C15_smg, pvecback[pba->index_bg_C15_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_C16_smg, pvecback[pba->index_bg_C16_smg]);
		}

		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_1_smg, pvecback[pba->index_bg_lambda_1_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_2_smg, pvecback[pba->index_bg_lambda_2_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_3_smg, pvecback[pba->index_bg_lambda_3_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_4_smg, pvecback[pba->index_bg_lambda_4_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_5_smg, pvecback[pba->index_bg_lambda_5_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_6_smg, pvecback[pba->index_bg_lambda_6_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_7_smg, pvecback[pba->index_bg_lambda_7_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_8_smg, pvecback[pba->index_bg_lambda_8_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_9_smg, pvecback[pba->index_bg_lambda_9_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_10_smg, pvecback[pba->index_bg_lambda_10_smg]);
	  copy_to_background_table_smg(pba, i, pba->index_bg_lambda_11_smg, pvecback[pba->index_bg_lambda_11_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_cs2num_smg, pvecback[pba->index_bg_cs2num_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_cs2_smg, pvecback[pba->index_bg_cs2_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_G_eff_smg, pvecback[pba->index_bg_G_eff_smg]);
		copy_to_background_table_smg(pba, i, pba->index_bg_slip_eff_smg, pvecback[pba->index_bg_slip_eff_smg]);

		/* Here we update the minimum values of the stability quantities
		* test will be performed based on the lowest values
		*/
		if (a > pba->a_min_stability_test_smg && pba->parameters_tuned_smg == _TRUE_){
			if (pvecback[pba->index_bg_kinetic_D_smg] < pba->min_D_smg){
				pba->min_D_smg = pvecback[pba->index_bg_kinetic_D_smg];
				pba->a_min_D_smg = a;
			}
			if (pvecback[pba->index_bg_cs2_smg] <= pba->min_cs2_smg){
				pba->min_cs2_smg = pvecback[pba->index_bg_cs2_smg];
				pba->a_min_cs2_smg = a;
			}
			if (pvecback[pba->index_bg_M2_smg] < pba->min_M2_smg){
				pba->min_M2_smg = pvecback[pba->index_bg_M2_smg];
				pba->a_min_M2_smg = a;
			}
			if (pvecback[pba->index_bg_tensor_excess_smg] + 1. < pba->min_ct2_smg){
				pba->min_ct2_smg = 1. + pvecback[pba->index_bg_tensor_excess_smg];
				pba->a_min_ct2_smg = a;
			}
			if (pvecback[pba->index_bg_braiding_smg] < pba->min_bra_smg){
				pba->min_bra_smg = pvecback[pba->index_bg_braiding_smg];
			}
			if (pvecback[pba->index_bg_braiding_smg] > pba->max_bra_smg){
				pba->max_bra_smg = pvecback[pba->index_bg_braiding_smg];
			}
		}

	}

	class_call_except(
	  stability_tests_smg(pba, pvecback, pvecback_integration),
	  pba->error_message,
	  pba->error_message,
	  free(pvecback_derivs);free(pvecback);free(pvecback_integration);background_free(pba);
	);

	 /* Yet another (third!) loop to make sure the background table makes sense
	 */
	for (i=0; i < pba->bt_size; i++) {

		//Need to update pvecback
		class_call(background_at_tau(pba,
				 pba->tau_table[i],
				 long_info,
				 inter_normal,
				 &last_index, //should be no problem to use the same one as for the derivatives
				 pvecback),
		pba->error_message,
		pba->error_message);

		// NOTE: note that the derivative is now calculated w.r.t. loga, while our _prime are w.r.t. tau
		double factor = pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H];
		double d_over_dtau;

    //write the derivatives in the structure
    class_call(array_derivate_spline(pba->loga_table, // x_array
			       pba->bt_size, // int n_lines
			       pba->background_table, // array
			       pba->d2background_dloga2_table, // double * array_splined
			       pba->bg_size, // n_columns
			       pba->loga_table[i], // double x -> tau
			       &last_index, // int* last_index // this is something for the interpolation to talk to each other when using a loop
			       pvecback_derivs, // double * result
			       pba->bg_size, //result_size, from 1 to n_columns
			       pba->error_message),
	  pba->error_message,
	  pba->error_message);

    //cs2num'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_cs2num_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_cs2num_prime_smg, d_over_dtau);

    //D'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_kinetic_D_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_kinetic_D_prime_smg, d_over_dtau);

		if (pba->field_evolution_smg == _TRUE_) {
			//D_over_phiphi'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_kinetic_D_over_phiphi_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_kinetic_D_over_phiphi_prime_smg, d_over_dtau);
		}

    //A9'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_A9_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_A9_prime_smg, d_over_dtau);

    //A10'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_A10_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_A10_prime_smg, d_over_dtau);

    //A12'
	  d_over_dtau = factor*pvecback_derivs[pba->index_bg_A12_smg];

    //A13'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_A13_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_A13_prime_smg, d_over_dtau);

		if (pba->field_evolution_smg == _TRUE_) {
	    //C9'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_C9_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_C9_prime_smg, d_over_dtau);

	    //C10'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_C10_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_C10_prime_smg, d_over_dtau);

	    //C12'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_C12_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_C12_prime_smg, d_over_dtau);

	    //C13'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_C13_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_C13_prime_smg, d_over_dtau);
		}

    //lambda_2'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_lambda_2_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_2_prime_smg, d_over_dtau);

    //lambda_8'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_lambda_8_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_8_prime_smg, d_over_dtau);

		//lambda_9'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_lambda_9_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_9_prime_smg, d_over_dtau);

		//lambda_11'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_lambda_11_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_11_prime_smg, d_over_dtau);


	  // check if any of the values becomes nan
	  int j = 0;
	  while (j < pba->bg_size){
	  	class_test_except(isnan(pvecback[j]) && (pba->parameters_tuned_smg == _TRUE_),
	             pba->error_message,
	             free(pvecback_derivs);free(pvecback);free(pvecback_integration);background_free(pba),
	             "pvecback[%i] = %e at a = %e in background!",j,pvecback[j],pvecback[pba->index_bg_a]);
	 		j++;
		}

	}

	free(pvecback_derivs);  //free the structure

	return _SUCCESS_;

}


/**
* Send _smg information to the standard output.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
* @return the error status
*/
int background_print_stdout_smg(
			  											  struct background *pba,
																double * pvecback,
																double * pvecback_integration
															  ) {

	printf(" -> Omega_smg = %f, wanted %f ",pvecback[pba->index_bg_rho_smg]/pvecback[pba->index_bg_rho_crit], pba->Omega0_smg);
	if(pba->has_lambda == _TRUE_)
		printf(", Omega_Lambda = %f", pba->Omega0_lambda);
	printf("\n");
	if (pba->background_verbose > 3) {
		printf("Minimal stability values: cs2 = %g, ct2 = %g, D = %g, M2 = %g \n",pba->min_cs2_smg,pba->min_ct2_smg,pba->min_D_smg,pba->min_M2_smg);
	}

	if (pba->field_evolution_smg == _TRUE_){
	    pba->xi_0_smg = pvecback_integration[pba->index_bi_phi_prime_smg]*pvecback[pba->index_bg_H]/pow(pba->H0,2);
	    pba->phi_0_smg = pvecback_integration[pba->index_bi_phi_smg];
	}

	gravity_models_print_stdout_smg(pba);

	return _SUCCESS_;
}


/**
* For each integrated variable fix the initial condition.
*
* @param pba                  Input: pointer to background structure
* @param a                    Input: scale factor
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
* @param ptr_rho_rad          Input: pointer to the density of radiation
* @return the error status
*/
int background_initial_conditions_smg(
        															struct background *pba,
																			double a,
																			double * pvecback,
        															double * pvecback_integration,
																			double * ptr_rho_rad
														  			  ) {

	double rho_rad = *ptr_rho_rad;
	double phi_scale, V_scale,p1,p2,p3; //smg related variables
	int i = 0;

	/** - fix initial value of modified gravity
	* run over all possible model cases
	*/

	pba->initial_conditions_set_smg = _FALSE_;

	//default value, can override later
	if (pba->M2_evolution_smg ==_TRUE_){
	  pvecback_integration[pba->index_bi_delta_M2_smg] = 0.;
	}

	class_call(gravity_models_initial_conditions_smg(pba, a, pvecback, pvecback_integration, &rho_rad),
						 pba->error_message,
						 pba->error_message);

	if (pba->field_evolution_smg == _TRUE_){
	  if (pba->background_verbose>3)
	    printf(" -> Initial conditions: phi = %e, phi' = %e \n",pvecback_integration[pba->index_bi_phi_smg],pvecback_integration[pba->index_bi_phi_prime_smg]);

	  class_test_except(!isfinite(pvecback_integration[pba->index_bi_phi_smg]) || !isfinite(pvecback_integration[pba->index_bi_phi_smg]),
		 								 pba->error_message,
	   						 	   free(pvecback);free(pvecback_integration);background_free(pba),
		 							   "initial phi = %e phi_prime = %e -> check initial conditions",
		 							 	 pvecback_integration[pba->index_bi_phi_smg],pvecback_integration[pba->index_bi_phi_smg]);
	 }
	if (pba->M2_evolution_smg == _TRUE_)
	  if (pba->background_verbose>3)
	    printf(" -> Initial conditions: delta_M2 = %e \n",pvecback_integration[pba->index_bi_delta_M2_smg]);

	return _SUCCESS_;
}


/**
 * Subroutine for formatting background _smg output
 *
 * @param pba                  Input: pointer to background structure
 * @param titles               Ouput: name of columns when printing the background table
 * @return the error status
 */
int background_store_columntitles_smg(
																		  struct background *pba,
																			char titles[_MAXTITLESTRINGLENGTH_]
																		  ) {

	class_store_columntitle(titles,"(.)rho_smg",_TRUE_);
  class_store_columntitle(titles,"(.)p_smg",_TRUE_);

	if (pba->output_background_smg >= 1){
    class_store_columntitle(titles,"M*^2_smg",_TRUE_);
    class_store_columntitle(titles,"D_M*^2_smg",_TRUE_);
    class_store_columntitle(titles,"kineticity_smg",_TRUE_);
    class_store_columntitle(titles,"braiding_smg",_TRUE_);
    class_store_columntitle(titles,"tensor_excess_smg",_TRUE_);
    class_store_columntitle(titles,"M2_running_smg",_TRUE_);
    class_store_columntitle(titles,"beyond_horndeski_smg",_TRUE_);
    class_store_columntitle(titles,"c_s^2",_TRUE_);
    class_store_columntitle(titles,"kin (D)",_TRUE_);
    class_store_columntitle(titles,"Current",pba->field_evolution_smg);
    class_store_columntitle(titles,"Shift",pba->field_evolution_smg);
  }

  if (pba->output_background_smg >= 2){
    class_store_columntitle(titles,"phi_smg",pba->field_evolution_smg);
    class_store_columntitle(titles,"phi'",pba->field_evolution_smg);
    class_store_columntitle(titles,"phi''",pba->field_evolution_smg);
    class_store_columntitle(titles,"E0",pba->field_evolution_smg);
    class_store_columntitle(titles,"E1",pba->field_evolution_smg);
    class_store_columntitle(titles,"E2",pba->field_evolution_smg);
    class_store_columntitle(titles,"E3",pba->field_evolution_smg);
    class_store_columntitle(titles,"P0",pba->field_evolution_smg);
    class_store_columntitle(titles,"P1",pba->field_evolution_smg);
    class_store_columntitle(titles,"P2",pba->field_evolution_smg);
    class_store_columntitle(titles,"R0",pba->field_evolution_smg);
    class_store_columntitle(titles,"R1",pba->field_evolution_smg);
    class_store_columntitle(titles,"R2",pba->field_evolution_smg);
    class_store_columntitle(titles,"G_eff_smg",_TRUE_);
    class_store_columntitle(titles,"slip_eff_smg",_TRUE_);
  }

  if (pba->output_background_smg >= 3){
    class_store_columntitle(titles,"kineticity_prime_smg",_TRUE_);
    class_store_columntitle(titles,"braiding_prime_smg",_TRUE_);
    class_store_columntitle(titles,"kineticity_over_phiphi_smg",pba->field_evolution_smg);
    class_store_columntitle(titles,"braiding_over_phi_smg",pba->field_evolution_smg);
    class_store_columntitle(titles,"beyond_horndeski_over_phi_smg",pba->field_evolution_smg);
    class_store_columntitle(titles,"kin_over_phiphi (D)",pba->field_evolution_smg);
    class_store_columntitle(titles,"A0",_TRUE_);
    class_store_columntitle(titles,"A1",_TRUE_);
    class_store_columntitle(titles,"A2",_TRUE_);
    class_store_columntitle(titles,"A3",_TRUE_);
    class_store_columntitle(titles,"A4",_TRUE_);
    class_store_columntitle(titles,"A5",_TRUE_);
    class_store_columntitle(titles,"A6",_TRUE_);
    class_store_columntitle(titles,"A7",_TRUE_);
    class_store_columntitle(titles,"A8",_TRUE_);
    class_store_columntitle(titles,"A9",_TRUE_);
    class_store_columntitle(titles,"A10",_TRUE_);
    class_store_columntitle(titles,"A11",_TRUE_);
    class_store_columntitle(titles,"A12",_TRUE_);
    class_store_columntitle(titles,"A13",_TRUE_);
    class_store_columntitle(titles,"A14",_TRUE_);
    class_store_columntitle(titles,"A15",_TRUE_);
    class_store_columntitle(titles,"A16",_TRUE_);
    class_store_columntitle(titles,"B0",pba->field_evolution_smg);
    class_store_columntitle(titles,"B1",pba->field_evolution_smg);
    class_store_columntitle(titles,"B2",pba->field_evolution_smg);
    class_store_columntitle(titles,"B3",pba->field_evolution_smg);
    class_store_columntitle(titles,"B4",pba->field_evolution_smg);
    class_store_columntitle(titles,"B5",pba->field_evolution_smg);
    class_store_columntitle(titles,"B6",pba->field_evolution_smg);
    class_store_columntitle(titles,"B7",pba->field_evolution_smg);
    class_store_columntitle(titles,"B8",pba->field_evolution_smg);
    class_store_columntitle(titles,"B9",pba->field_evolution_smg);
    class_store_columntitle(titles,"B10",pba->field_evolution_smg);
    class_store_columntitle(titles,"B11",pba->field_evolution_smg);
    class_store_columntitle(titles,"B12",pba->field_evolution_smg);
    class_store_columntitle(titles,"C0",pba->field_evolution_smg);
    class_store_columntitle(titles,"C1",pba->field_evolution_smg);
    class_store_columntitle(titles,"C2",pba->field_evolution_smg);
    class_store_columntitle(titles,"C3",pba->field_evolution_smg);
    class_store_columntitle(titles,"C4",pba->field_evolution_smg);
    class_store_columntitle(titles,"C5",pba->field_evolution_smg);
    class_store_columntitle(titles,"C6",pba->field_evolution_smg);
    class_store_columntitle(titles,"C7",pba->field_evolution_smg);
    class_store_columntitle(titles,"C8",pba->field_evolution_smg);
    class_store_columntitle(titles,"C9",pba->field_evolution_smg);
    class_store_columntitle(titles,"C10",pba->field_evolution_smg);
    class_store_columntitle(titles,"C11",pba->field_evolution_smg);
    class_store_columntitle(titles,"C12",pba->field_evolution_smg);
    class_store_columntitle(titles,"C13",pba->field_evolution_smg);
    class_store_columntitle(titles,"C14",pba->field_evolution_smg);
    class_store_columntitle(titles,"C15",pba->field_evolution_smg);
    class_store_columntitle(titles,"C16",pba->field_evolution_smg);
    class_store_columntitle(titles,"lambda_1",_TRUE_);
    class_store_columntitle(titles,"lambda_2",_TRUE_);
    class_store_columntitle(titles,"lambda_3",_TRUE_);
    class_store_columntitle(titles,"lambda_4",_TRUE_);
    class_store_columntitle(titles,"lambda_5",_TRUE_);
    class_store_columntitle(titles,"lambda_6",_TRUE_);
    class_store_columntitle(titles,"lambda_7",_TRUE_);
    class_store_columntitle(titles,"lambda_8",_TRUE_);
    class_store_columntitle(titles,"lambda_9",_TRUE_);
    class_store_columntitle(titles,"lambda_10",_TRUE_);
    class_store_columntitle(titles,"lambda_11",_TRUE_);
    class_store_columntitle(titles,"lambda_2_p",_TRUE_);
    class_store_columntitle(titles,"lambda_8_p",_TRUE_);
    class_store_columntitle(titles,"lambda_9_p",_TRUE_);
    class_store_columntitle(titles,"lambda_11_p",_TRUE_);
    class_store_columntitle(titles,"cs2num",_TRUE_);
    class_store_columntitle(titles,"cs2num_p",_TRUE_);
  }

	return _SUCCESS_;
}


/**
 * Subroutine for writing the background _smg output
 *
 * @param pba                  Input: pointer to background structure
 * @param pvecback             Input: vector of background quantities
 * @param dataptr              Ouput: pointer to 1d array storing all the background table
 * @param ptr_storeidx         Ouput: pointer to index with number of columns
 * @return the error status
 */
int background_output_data_smg(
															 struct background *pba,
        							 				 double * pvecback,
															 double * dataptr,
															 int * ptr_storeidx
														   ) {

	int storeidx = *ptr_storeidx;

	class_store_double(dataptr,pvecback[pba->index_bg_rho_smg],_TRUE_,storeidx);
  class_store_double(dataptr,pvecback[pba->index_bg_p_smg],_TRUE_,storeidx);

	if (pba->output_background_smg >= 1){
		class_store_double(dataptr,pvecback[pba->index_bg_M2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_delta_M2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_kineticity_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_braiding_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_tensor_excess_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_M2_running_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_beyond_horndeski_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_cs2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_kinetic_D_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_current_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_shift_smg],pba->field_evolution_smg,storeidx);
	}

	if (pba->output_background_smg >= 2){
		class_store_double(dataptr,pvecback[pba->index_bg_phi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_phi_prime_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_phi_prime_prime_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_E0_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_E1_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_E2_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_E3_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_P0_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_P1_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_P2_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_R0_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_R1_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_R2_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_G_eff_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_slip_eff_smg],_TRUE_,storeidx);
	}

	if (pba->output_background_smg >= 3){
		class_store_double(dataptr,pvecback[pba->index_bg_kineticity_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_braiding_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_kineticity_over_phiphi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_braiding_over_phi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_beyond_horndeski_over_phi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_kinetic_D_over_phiphi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A0_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A1_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A3_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A4_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A5_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A6_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A7_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A8_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A9_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A10_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A11_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A12_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A13_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A14_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A15_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A16_smg],_TRUE_,storeidx);

		if (pba->field_evolution_smg == _TRUE_) {
			class_store_double(dataptr,pvecback[pba->index_bg_B0_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B1_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B2_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B3_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B4_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B5_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B6_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B7_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B8_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B9_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B10_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B11_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B12_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C0_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C1_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C2_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C2_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C3_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C4_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C5_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C6_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C7_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C8_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C9_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C10_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C11_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C12_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C13_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C15_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C16_smg],_TRUE_,storeidx);
		}

		class_store_double(dataptr,pvecback[pba->index_bg_lambda_1_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_3_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_4_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_5_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_6_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_7_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_8_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_9_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_10_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_11_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_2_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_8_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_9_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_11_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_cs2num_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_cs2num_prime_smg],_TRUE_,storeidx);
	}

	*ptr_storeidx = storeidx;

	return _SUCCESS_;
}


/**
 * Derivatives of {B} _smg quantities
 *
 * @param pba                  Input: pointer to background structure
 * @param a                    Input: scale factor
 * @param pvecback             Input: vector of background quantities
 * @param y                    Input: current vector of integrated quantities (with index_bi)
 * @param dy                   Output: current derivative of y w.r.t log(a)
 * @return the error status
 */
int background_derivs_smg(
			  									struct background *pba,
													double a,
													double * pvecback,
													double * y,
													double * dy
												  ) {

	double H = pvecback[pba->index_bg_H];

	/** - calculate /f$ \rho'(\tau)= -3aH (1+w) \rho /f$ written as \f$ d\rho/dloga = \rho' / (aH) \f$ */
	if (pba->rho_evolution_smg == _TRUE_){
	  dy[pba->index_bi_rho_smg] = pvecback[pba->index_bg_rho_prime_smg]/a/H;
	}

	/** - Scalar field equation: \f$ \phi''(t) + 2 a H \phi'(t) + a^2 dV = 0 \f$
				(note H is wrt cosmic time) **/
	if(pba->field_evolution_smg){
	  dy[pba->index_bi_phi_smg] = y[pba->index_bi_phi_prime_smg]/a/H;
	  dy[pba->index_bi_phi_prime_smg] = pvecback[pba->index_bg_phi_prime_prime_smg]/a/H;
	}
	/** - Planck mass equation (if parameterization in terms of alpha_m **/
	if (pba->M2_evolution_smg == _TRUE_)
	  dy[pba->index_bi_delta_M2_smg] = pvecback[pba->index_bg_M2_running_smg]*(y[pba->index_bi_delta_M2_smg] + 1);   //in this case the running has to be integrated (eq 3.3 of 1404.3713 yields M2' = aH\alpha_M)

	return _SUCCESS_;
}


/**
* Stability tests for smg.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Input: vector of background quantities
* @param pvecback_integration Input: vector of background quantities to be integrated
* @return the error status
*/
int stability_tests_smg(
        								struct background *pba,
        								double * pvecback,
												double * pvecback_integration
											  ) {

	/* Horndeski stability tests
	* only if not overriden
	* and model is tuned!
	*/
	if ((pba->parameters_tuned_smg == _TRUE_) &&
	    (pba->skip_stability_tests_smg == _FALSE_)){

	  class_test(pba->min_D_smg <= -fabs(pba->D_safe_smg),
	      pba->error_message,
	      "Ghost instability for scalar field perturbations with minimum D=%g at a=%e\n", pba->min_D_smg, pba->a_min_D_smg);
	  class_test(pba->min_cs2_smg < -fabs(pba->cs2_safe_smg),
	      pba->error_message,
	      "Gradient instability for scalar field perturbations with minimum c_s^2=%g at a=%e\n", pba->min_cs2_smg, pba->a_min_cs2_smg);
	  class_test(pba->min_M2_smg < -fabs(pba->M2_safe_smg),
	      pba->error_message,
	      "Ghost instability for metric tensor perturbations with minimum M*^2=%g at a=%e\n", pba->min_M2_smg, pba->a_min_M2_smg);
	  class_test(pba->min_ct2_smg < -fabs(pba->ct2_safe_smg),
	      pba->error_message,
	      "Gradient instability for metric tensor perturbations with minimum c_t^2=%g at a=%e\n",pba->min_ct2_smg,pba->a_min_ct2_smg);

	 }

	return _SUCCESS_;
}


/**
* Numerical derivatives of the alphas.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Input: vector of background quantities
* @param pvecback_derivs      Output: vector of derivatives
* @param i                    Input: counter for the time step
* @return the error status
*/
int derivatives_alphas_smg(
													 struct background *pba,
													 double * pvecback,
													 double * pvecback_derivs,
													 int i
													 ) {

	/* -> write in the table (overwrite the alpha time derivatives, which were set to nan in background_functions)
	 * direction of copy: add the corresponding indices to the coordinates
	 * thing to be copied: the address (&) to the element of pvecback corresponding to the index we want
	 * size: just a single double number
	 * -> repeat for all necessary quantities
	 */

	/* necessary for calling array_interpolate(), but never used */
	int last_index;

	// NOTE: note that the derivative is now calculated w.r.t. log(a), while our _prime are w.r.t. tau
	double factor = pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H];
	double d_over_dtau;

	/*NOTE: here we compute the derivatives of quantities coputed in background_gravity_functions_smg during the integration.
	 * for quantities that depend on these derivatives (e.g. the gamma_i functions determining the effective mass)
	 * there is an additional loop at the end of background_solve
	 */

	// Kineticity'
	d_over_dtau = factor*pvecback_derivs[pba->index_bg_kineticity_smg];
	copy_to_background_table_smg(pba, i, pba->index_bg_kineticity_prime_smg, d_over_dtau);

   //Braiding'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_braiding_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_braiding_prime_smg, d_over_dtau);

   //Planck mass run rate'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_M2_running_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_M2_running_prime_smg, d_over_dtau);

   //Tensor excess'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_tensor_excess_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_tensor_excess_prime_smg, d_over_dtau);

  //Beyond horndeski'
	d_over_dtau = factor*pvecback_derivs[pba->index_bg_beyond_horndeski_smg];
	copy_to_background_table_smg(pba, i, pba->index_bg_beyond_horndeski_prime_smg, d_over_dtau);

	if (pba->field_evolution_smg == _TRUE_) {
		//Braiding_over_phi'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_braiding_over_phi_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_braiding_over_phi_prime_smg, d_over_dtau);

		//Beyond_horndeski_over_phi'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_beyond_horndeski_over_phi_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_beyond_horndeski_over_phi_prime_smg, d_over_dtau);
	}

   //H''
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_H_prime];
	 copy_to_background_table_smg(pba, i, pba->index_bg_H_prime_prime, d_over_dtau);

   // p_tot_wo_smg'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_p_tot_wo_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_p_tot_wo_prime_smg, d_over_dtau);

   // p_smg'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_p_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_p_prime_smg, d_over_dtau);

	// Planck's mass running
	// Only need to compute it if neither self consistent field evolution nor evolving M2 in terms of alpha_M
	// check equation 3.3 of Bellini & Sawicki 2014

	if (pba->field_evolution_smg == _FALSE_ && pba->M2_evolution_smg == _FALSE_){
		double alpha_M = pvecback_derivs[pba->index_bg_delta_M2_smg]/pvecback[pba->index_bg_M2_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_M2_running_smg, alpha_M);
	}

	if(pba->background_verbose > 15 && fabs(1. - pvecback[pba->index_bg_H_prime]/pvecback_derivs[pba->index_bg_H]/factor)>1e-8)
printf("a = %g, (delta H')/H' = %g \n", pvecback[pba->index_bg_a], 1. - pvecback[pba->index_bg_H_prime]/pvecback_derivs[pba->index_bg_H]/factor);

	return _SUCCESS_;
}


/**
* Copy to the background table _smg quantities.
*
* @param pba                  Input: pointer to background structure
* @param row                  Input: table row
* @param column               Input: table column
* @param value                Input: value to copy
* @return the error status
*/
int copy_to_background_table_smg(
																 struct background *pba,
                                 int row,
                                 int column,
                                 double value
															   ) {

	/* needed for growing table */
	void * memcopy_result;

	memcopy_result = memcpy(pba->background_table + row*pba->bg_size + column,
	                        &value, 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + row*pba->bg_size + column,
	           pba->error_message, "cannot copy data back to pba->background_table");

  return _SUCCESS_;
}
