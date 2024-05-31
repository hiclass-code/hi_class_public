#ifndef __INPUT_SMG__
#define __INPUT_SMG__

#include "common.h"
#include "input.h"
#include "gravity_models_smg.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int input_warnings_smg(struct perturbations * ppt,
                       int input_verbose);

int input_read_parameters_smg(struct file_content * pfc,
                              struct precision * ppr,
                              struct background * pba,
                              struct perturbations * ppt,
                              ErrorMsg errmsg);

int input_readjust_precision_smg(struct precision * ppr);

int input_default_params_smg(struct background * pba,
                             struct perturbations * ppt);


#ifdef __cplusplus
}
#endif

#endif
