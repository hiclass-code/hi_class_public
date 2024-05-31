#ifndef __GRAVITY_MODELS_SMG__
#define __GRAVITY_MODELS_SMG__

#include "common.h"
#include "input.h"
#include "background.h"
#include "rootfinder.h"

/**
 * structure containing the G's and F's functions, for Horndeski
 * and beyond, as well as their derivatives w.r.t. the scalar field
 * and its canonical kinetic term.
 * All of them are inizialized to zero except G4 = 1./2.
 * This configuration reproduces General Relativity.
 */
struct G_functions_and_derivs {
  /* G2s */
  double G2;
  double G2_X;
  double G2_XX;
  double G2_XXX;
  double G2_phi;
  double G2_Xphi;
  double G2_XXphi;
  double G2_phiphi;
  double G2_Xphiphi;

  /* G3s */
  double G3_X;
  double G3_XX;
  double G3_XXX;
  double G3_phi;
  double G3_Xphi;
  double G3_XXphi;
  double G3_phiphi;
  double G3_Xphiphi;
  double G3_phiphiphi;

  /* G4s */
  double G4;
  double DG4;
  double G4_X;
  double G4_XX;
  double G4_XXX;
  double G4_XXXX;
  double G4_phi;
  double G4_Xphi;
  double G4_XXphi;
  double G4_XXXphi;
  double G4_phiphi;
  double G4_Xphiphi;
  double G4_XXphiphi;
  double G4_phiphiphi;
  double G4_Xphiphiphi;

  /* G5s */
  double G5_X;
  double G5_XX;
  double G5_XXX;
  double G5_XXXX;
  double G5_phi;
  double G5_Xphi;
  double G5_XXphi;
  double G5_XXXphi;
  double G5_phiphi;
  double G5_Xphiphi;
  double G5_XXphiphi;
  double G5_phiphiphi;
  double G5_Xphiphiphi;

  /* F4s */
  double F4;
  double F4_X;
  double F4_XX;
  double F4_XXX;
  double F4_phi;
  double F4_Xphi;
  double F4_XXphi;
  double F4_phiphi;
  double F4_Xphiphi;

  /* F4s */
  double F5;
  double F5_X;
  double F5_XX;
  double F5_XXX;
  double F5_phi;
  double F5_Xphi;
  double F5_XXphi;
  double F5_phiphi;
  double F5_Xphiphi;
};

/**
 * Default values for the G_functions_and_derivs structure.
 */
#define DEFAULT_G_FUNCTIONS_AND_DERIVS {                         \
  /* G2s */                                                      \
  0., 0., 0., 0., 0., 0., 0., 0., 0.,                            \
  /* G3s */                                                      \
  0., 0., 0., 0., 0., 0., 0., 0., 0.,                            \
  /* G4s */                                                      \
  1./2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
  /* G5s */                                                      \
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,            \
  /* F4s */                                                      \
  0., 0., 0., 0., 0., 0., 0., 0., 0.,                            \
  /* F5s */                                                      \
  0., 0., 0., 0., 0., 0., 0., 0., 0.                             \
}


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


int gravity_models_gravity_properties_smg(struct file_content * pfc,
                                          struct background *pba,
                                          char string1[_ARGUMENT_LENGTH_MAX_],
                                          int has_tuning_index_smg,
                                          int has_dxdy_guess_smg,
                                          ErrorMsg errmsg);

int gravity_models_expansion_properties_smg(struct file_content * pfc,
                                            struct background *pba,
                                            char string1[_ARGUMENT_LENGTH_MAX_],
                                            ErrorMsg errmsg);

int gravity_models_get_Gs_smg(struct background *pba,
                              double a,
                              double * pvecback_B,
                              struct G_functions_and_derivs *g_fun);

int gravity_models_get_back_par_smg(struct background *pba,
                                    double a,
                                    double * pvecback,
                                    double * pvecback_B);

int gravity_models_get_alphas_par_smg(struct background *pba,
                                      double a,
                                      double * pvecback,
                                      double * pvecback_B);

int gravity_models_initial_conditions_smg(struct background *pba,
											    								double a,
																			    double * pvecback,
        															    double * pvecback_integration,
																			    double * ptr_rho_rad);

int gravity_models_print_stdout_smg(struct background *pba);


#ifdef __cplusplus
}
#endif

#endif
