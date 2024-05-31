/** @file input_smg.c Documented input_smg module
 *
 * Emilio Bellini, Ignacy Sawicki, Miguel Zumalacarregui, 2024
 *
 * Additional functions for the input module.
 * It contains all the hi_class related functions (_smg)
 * that are used by input.c. In this way the main hi_class
 * modifications are stored here and the standard Class modules
 * remain cleaner.
 *
 * The following nomenclature has been adopted:
 *
 * -# all the functions end with "_smg" to make them easily
 *    recognizable
 * -# all the functions starting with "input_" are
 *    directly called by input.c or the classy wrapper
 * -# all the functions that do not start with "input_"
 *    are only used internally in input_smg.c
 */

#include "input_smg.h"
#include "background.h"
#include "perturbations.h"


/**
 * Place to put the warnings related to the hi_class input parameters.
 *
 * @param ppt              Input: pointer to perturbation structure
 * @param input_verbose    Input: input verbosity
 * @return the error status
 */
int input_warnings_smg(
                       struct perturbations * ppt,
                       int input_verbose
                       ) {

  /* Here we put a warning as we want to encourage hi_class users to get
  h_prime from the trace of the Einstein ij equation than from the Einstein 00
  equation. This is because the Einstein 00 equation has a gauge dependent
  singularity that can be removed using the trace of the Einstein ij equation.
  */
  if (input_verbose > 10) {
    if (ppt->get_h_from_trace == _FALSE_) {
      printf("\n");
      printf("WARNING: you set get_h_from_trace to False.\n");
      printf("While this is still accepted in hi_class, it can cause gauge dependent\n");
      printf("singularities if your model crosses alphaB=2. For this reason in\n");
      printf("future versions of the code this option will be removed and the\n");
      printf("Einstein 00 equation will be used only to set the initial conditions\n");
      printf("for h_prime and as a test to check that it is satisfied during the evolution.\n");
      printf("On the other hand this is a safe choice if you want very large k modes\n");
      printf("(typically k>10 Mpc^-1), where the constraint and the dynamical equations\n");
      printf("disagree by a non negligible amount in some of the models studied.\n");
      printf("\n");
    }
    else if (ppt->get_h_from_trace == _TRUE_) {
      printf("\n");
      printf("WARNING: you set get_h_from_trace to True.\n");
      printf("While this will be the default option in future versions of hi_class it might\n");
      printf("be safer to set get_h_from_trace to False if you want very large k modes\n");
      printf("(typically k>10 Mpc^-1). In this regime the constraint and the dynamical \n");
      printf("equations disagree by a non negligible amount some of the models studied.\n");
      printf("\n");
    }
  }

  return _SUCCESS_;
}


/**
 * Parse the hi_class parameters.
 *
 * @param pfc              Input: pointer to local structure
 * @param ppr              Input/Output: pointer to precision structure
 * @param pba              Input/Output: pointer to background structure
 * @param ppt              Input/Output: pointer to perturbation structure
 * @param errmsg           Input: Error message
 * @return the error status
 */
int input_read_parameters_smg(
                              struct file_content * pfc,
                              struct precision * ppr,
                              struct background * pba,
                              struct perturbations * ppt,
                              ErrorMsg errmsg
                              ) {

  int flag1, flag2, flag3;
  double param1, param2;
  int param3;
  int entries_read;
  int int1, n;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_ARGUMENT_LENGTH_MAX_];

  pba->has_smg = _TRUE_;  //default is _FALSE_

  /** Main flag for the quasi-static approximation scheme */

  class_call(parser_read_string(pfc,"method_qs_smg",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"automatic") != NULL) || (strstr(string1,"a") != NULL) || (strstr(string1,"A") != NULL)) {
      ppt->method_qs_smg = automatic;
    }

    if ((strstr(string1,"fully_dynamic") != NULL) || (strstr(string1,"fd") != NULL) || (strstr(string1,"FD") != NULL)) {
      ppt->method_qs_smg = fully_dynamic;
    }

    if ((strstr(string1,"quasi_static") != NULL) || (strstr(string1,"qs") != NULL) || (strstr(string1,"QS") != NULL)) {
      ppt->method_qs_smg = quasi_static;
    }

    if ((strstr(string1,"fully_dynamic_debug") != NULL) || (strstr(string1,"fdd") != NULL) || (strstr(string1,"FDD") != NULL)) {
      ppt->method_qs_smg = fully_dynamic_debug;
    }

    if ((strstr(string1,"quasi_static_debug") != NULL) || (strstr(string1,"qsd") != NULL) || (strstr(string1,"QSD") != NULL)) {
      ppt->method_qs_smg = quasi_static_debug;
    }
  }

  class_call(parser_read_string(pfc, "use_pert_var_deltaphi_smg", &string1, &flag1, errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_){
    if((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)){
      ppt->use_pert_var_deltaphi_smg = _TRUE_;
    }
    else{
      ppt->use_pert_var_deltaphi_smg = _FALSE_;
    }
  }

  /** Read tuning parameter and guess for the parameter variation range
   * These can be adjusted later on a model basis
   */
  int has_tuning_index_smg, has_dxdy_guess_smg;

  class_call(parser_read_int(pfc,"tuning_index_smg",&param3,&flag3,errmsg),
             errmsg,
             errmsg);
  if (flag3 == _TRUE_){
    pba->tuning_index_smg = param3;
  }
  has_tuning_index_smg = flag3;

  class_call(parser_read_double(pfc,"tuning_dxdy_guess_smg",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  if (flag2 == _TRUE_){
    pba->tuning_dxdy_guess_smg = param2;
  }
  has_dxdy_guess_smg = flag2;
  if (has_dxdy_guess_smg == _FALSE_)
    pba->tuning_dxdy_guess_smg = 1;

  class_call(parser_read_string(pfc,"gravity_model",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_test(flag1 == _FALSE_,
              errmsg,
              "gravity_model not read, you should specify one!\n");

  /** Fix flags for each gravity model */
    class_call(gravity_models_gravity_properties_smg(pfc, pba, string1, has_tuning_index_smg, has_dxdy_guess_smg, errmsg),
              errmsg,
              errmsg);
  // end of loop over models

  if(pba->field_evolution_smg == _TRUE_){

    //TODO: include generic stuff for covariant theories

  }
  else {
    //if no self-consistent evolution, need a parameterization for Omega_smg

    class_test(ppt->use_pert_var_deltaphi_smg==_TRUE_,
               errmsg,
               "It is not consistent to evolve delta_phi_smg and choose parametrized models.");

    class_call(parser_read_string(pfc,"expansion_model",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_test(flag1 == _FALSE_,
                errmsg,
                "expansion_model not read, you should specify one!\n");

    /** Fix flags for each expansion model */
    class_call(gravity_models_expansion_properties_smg(pfc, pba, string1, errmsg),
               errmsg,
               errmsg);

  }

  /** Other generic specifications:
   * - whether stability tests are skipped (skip_stability_tests_smg) or softened (cs2_safe_smg)
   * - thresholds for approximations in the cubic Friedmann equation
   * - add a value to have better behaved perturbations
   * - approximations in the perturbations
   */

  class_read_double("cs2_safe_smg",pba->cs2_safe_smg);
  class_read_double("D_safe_smg",pba->D_safe_smg);
  class_read_double("ct2_safe_smg",pba->ct2_safe_smg);
  class_read_double("M2_safe_smg",pba->M2_safe_smg);

  class_read_double("pert_ic_tolerance_smg",ppr->pert_ic_tolerance_smg);
  class_read_double("pert_ic_ini_z_ref_smg",ppr->pert_ic_ini_z_ref_smg);
  class_read_double("pert_ic_regulator_smg",ppr->pert_ic_regulator_smg);
  class_read_double("pert_qs_ic_tolerance_test_smg",ppr->pert_qs_ic_tolerance_test_smg);

  class_read_double("a_min_stability_test_smg",pba->a_min_stability_test_smg);

  class_read_double("kineticity_safe_smg",pba->kineticity_safe_smg); // minimum value of the kineticity (to avoid trouble)
  class_read_double("min_a_pert_smg",ppr->min_a_pert_smg);

  class_call(parser_read_string(pfc,
			  "skip_stability_tests_smg",
			  &string1,
			  &flag1,
			  errmsg),
	errmsg,
	errmsg);

  if (flag1 == _TRUE_){
    if((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)){
      pba->skip_stability_tests_smg = _TRUE_;
    }
    else{
      pba->skip_stability_tests_smg = _FALSE_;
    }
  }

  //IC for perturbations
  class_call(parser_read_string(pfc,
			       "pert_initial_conditions_smg",
			       &string1,
			       &flag1,
			       errmsg),
	         errmsg,
	         errmsg);

  if (strcmp(string1,"single_clock") == 0) {
    ppt->pert_initial_conditions_smg = single_clock;
  }
  else if (strcmp(string1,"gravitating_attr") == 0) {
    ppt->pert_initial_conditions_smg = gravitating_attr;
  }
  else if (strcmp(string1,"zero") == 0) {
    ppt->pert_initial_conditions_smg = zero;
  }
  else if (strcmp(string1,"kin_only") == 0) {
    ppt->pert_initial_conditions_smg = kin_only;
  }
  else if (strcmp(string1,"ext_field_attr") == 0 ){//this is the default
    ppt->pert_initial_conditions_smg = ext_field_attr;
  }

  // else {
  //   if (ppt->perturbations_verbose > 1)
  //     printf(" Initial conditions for Modified gravity perturbations not specified, using default \n");
  // }

  /** re-assign shooting parameter (for no-tuning debug mode) */
  if (pba->Omega_smg_debug == 0)
    class_read_double("shooting_parameter_smg",pba->parameters_smg[pba->tuning_index_smg]);

  // test that the tuning is correct
  class_test(pba->tuning_index_smg >= pba->parameters_size_smg,
       errmsg,
       "Tuning index tuning_index_smg = %d is larger than the number of entries %d in parameters_smg. Check your .ini file.",
       pba->tuning_index_smg,pba->parameters_size_smg);

  //how much info on background.dat?
  class_read_double("output_background_smg",pba->output_background_smg);

  return _SUCCESS_;
}


/**
 * Sometimes hi_class has higher precision parameters. This is where
 * it is possible to readjust them.
 *
 * @param ppr              Input/Output: pointer to precision structure
 * @return the error status
 */
int input_readjust_precision_smg(
                                 struct precision * ppr
                                 ) {

  /** readjust some precision parameters for modified gravity */

  //otherwise problems with ISW effect
  if (ppr->perturbations_sampling_stepsize > 0.05)
    ppr->perturbations_sampling_stepsize=0.05;

  return _SUCCESS_;
}


/**
 * List of default hi_class parameters.
 *
 * @param pba              Input/Output: pointer to background structure
 * @param ppt              Input: pointer to perturbation structure
 * @return the error status
 */
int input_default_params_smg(
                             struct background * pba,
                             struct perturbations * ppt
                             ) {

  /** - background structure */

  pba->gravity_model_smg = propto_omega; /* gravitational model */
  pba->expansion_model_smg = lcdm; /*expansion model (only for parameterizations*/
  pba->Omega0_smg = 0.; /* Scalar field defaults */
  pba->M2_today_smg = 1.; //*Planck mass today*/
  pba->M2_tuning_smg = _FALSE_; //* Tune Planck mass?*/
  pba->Omega_smg_debug = 0;
  pba->field_evolution_smg = _FALSE_; /* does the model require solving the background equations? */
  pba->M2_evolution_smg = _FALSE_; /* does the model require integrating M2 from alpha_M? */
  pba->skip_stability_tests_smg = _FALSE_; /*if you want to skip the stability tests for the perturbations */
  pba->a_min_stability_test_smg = 0; /** < skip stability tests for a < a_min */

  pba->hubble_evolution = _TRUE_; /** dynamical evolution of Friedmann eq. */
  pba->hubble_friction = 3.; /** friction coefficient in H' equation: H' = ... + H_friction*(H^2 - rho_crit) [NOT ONLY IN SMG!] */
  pba->rho_evolution_smg = _FALSE_; /*does the model require to evolve the background energy density? (only for parameterizations)*/
  pba->wede_Omega_e_regularizer_smg = 1.e-10; /** regularize adding a tiny Omega_e */

  pba->kineticity_safe_smg = 0; /* value added to the kineticity, useful to cure perturbations at early time in some models */
  pba->cs2_safe_smg = 0; /* threshold to consider the sound speed of scalars negative in the stability check */
  pba->D_safe_smg = 0; /* threshold to consider the kinetic term of scalars negative in the stability check */
  pba->ct2_safe_smg = 0; /* threshold to consider the sound speed of tensors negative in the stability check */
  pba->M2_safe_smg = 0; /* threshold to consider the kinetic term of tensors (M2) negative in the stability check */


  /*set stability quantities to nonzero values*/
  pba->min_M2_smg = 1.e10;
  pba->min_ct2_smg = 1.e10;
  pba->min_D_smg = 1.e10;
  pba->min_cs2_smg = 1.e10;
  pba->a_min_M2_smg = 0.;
  pba->a_min_ct2_smg = 0.;
  pba->a_min_D_smg = 0.;
  pba->a_min_cs2_smg = 0.;

  pba->min_bra_smg = 4.;
  pba->max_bra_smg = 0.;
  pba->quintessence_w_safe_smg = 0;

  pba->parameters_smg = NULL;
  pba->parameters_size_smg = 0;
  pba->tuning_index_smg = 0;
  pba->tuning_dxdy_guess_smg = 1;

  pba->output_background_smg = 1; /**< amount of information printed onbackground.dat output */

  pba->has_smg= _FALSE_;
  pba->parameters_tuned_smg = _FALSE_;
  pba->shooting_failed = _FALSE_;
  pba->is_quintessence_smg = _FALSE_;
  pba->attractor_ic_smg = _TRUE_;  /* only read for those models in which it is implemented */
  pba->initial_conditions_set_smg = _FALSE_;


  /** - perturbations structure */

  ppt->get_h_from_trace=_FALSE_; /* Get h' from Einstein trace rather than 00 (not only _smg!!) */

  ppt->method_qs_smg=fully_dynamic;
  ppt->pert_initial_conditions_smg = ext_field_attr; /* default IC for perturbations in the scalar */

  ppt->use_pert_var_deltaphi_smg=_FALSE_;

  return _SUCCESS_;

}
