#ifndef __GRAVITY_FUNCTIONS_SMG__
#define __GRAVITY_FUNCTIONS_SMG__

#include "common.h"
#include "background.h"
#include "gravity_models_smg.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int gravity_functions_Es_from_Gs_smg(struct background *pba,
                                     double a,
                                     double * pvecback_B,
                                     double * pvecback,
                                     struct G_functions_and_derivs *pgf);

int gravity_functions_Ps_and_Rs_from_Gs_smg(struct background *pba,
                                            double a,
                                            double * pvecback_B,
                                            double * pvecback,
                                            struct G_functions_and_derivs *pgf);

int gravity_functions_building_blocks_from_Gs_smg(struct background *pba,
                                                  double a,
                                                  double * pvecback_B,
                                                  double * pvecback,
                                                  struct G_functions_and_derivs *pgf);

int gravity_functions_As_from_alphas_smg(struct background *pba,
                                         double * pvecback,
                                         double * pvecback_derivs);

int gravity_functions_Bs_from_Gs_smg(struct background *pba,
                                     double a,
                                     double * pvecback_B,
                                     double * pvecback,
                                     struct G_functions_and_derivs *pgf);

int gravity_functions_Cs_from_Bs_smg(struct background *pba,
                                     double * pvecback,
                                     double * pvecback_derivs);


#ifdef __cplusplus
}
#endif

#endif
