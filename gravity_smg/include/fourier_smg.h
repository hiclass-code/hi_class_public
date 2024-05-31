#ifndef __FOURIER_SMG__
#define __FOURIER_SMG__

#include "common.h"
#include "fourier.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int fourier_hmcode_Delta_v_0_smg(
  struct background *pba,
  double z_at_tau,
  double * Delta_v_0
);

#ifdef __cplusplus
}
#endif

#endif
