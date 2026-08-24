#ifndef __PERTURBATIONS_SMG__
#define __PERTURBATIONS_SMG__

#include "common.h"
#include "perturbations.h"
#include "rootfinder.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int get_gravity_coefficients_smg(struct background * pba,
                                 struct perturbations * ppt,
                                 double * pvecback,
                                 double * delM2, double * M2, double * kin, double * bra,
                                 double * ten, double * run, double * beh, double * res,
                                 double * cD, double * cK, double * cB, double * cH, double * c0,
                                 double * c1, double * c2, double * c3, double * c4,
                                 double * c5, double * c6, double * c7, double * c8,
                                 double * c9, double * c10, double * c11, double * c12,
                                 double * c13, double * c14, double * c15, double * c16,
                                 double * res_p, double *  cD_p, double *  cB_p, double *  cH_p,
                                 double * c9_p, double * c10_p, double * c12_p, double * c13_p);

int perturbations_tests_smg(struct precision * ppr,
                            struct background * pba,
                            struct perturbations * ppt);

int perturbations_define_indices_tp_smg(struct perturbations * ppt,
				                                int * index_type);

int perturbations_define_indices_mt_smg(struct perturbations_workspace * ppw,
				                                int * index_mt);

int perturbations_define_indices_ap_smg(struct perturbations_workspace * ppw,
  			                                int * index_ap);

int perturbations_define_indices_pt_smg(struct perturbations_workspace * ppw,
                                        struct perturbations_vector * ppv,
			                                  int * index_pt);

int perturbations_prepare_k_output_smg(struct perturbations * ppt);

int perturbations_print_variables_smg(struct precision * ppr,
                                      struct background * pba,
                                      struct perturbations * ppt,
                                      struct perturbations_workspace * ppw,
                                      double k,
                                      double tau,
                                      double * dataptr,
                                      int * ptr_storeidx);

int perturbations_get_h_prime_ic_from_00_smg(struct background * pba,
                                             struct perturbations_workspace * ppw,
                                             double k,
                                             double eta,
                                             double delta_rho_tot);

int perturbations_einstein_scalar_smg(struct precision * ppr,
                                      struct background * pba,
                                      struct thermodynamics * pth,
                                      struct perturbations * ppt,
                                      struct perturbations_workspace * ppw,
                                      double k,
                                      double tau,
                                      double * y);

int perturbations_einstein_tensor_smg(struct background * pba,
                                      struct perturbations_workspace * ppw,
                                      double k,
                                      double tau,
                                      double * y);

int perturbations_derivs_smg(struct perturbations * ppt,
                             struct perturbations_workspace * ppw,
                             struct perturbations_vector * pv,
                             double * dy,
                             double * pvecmetric);

int perturbations_get_approximation_qs_smg(struct precision * ppr,
                                           struct background * pba,
                                           struct perturbations * ppt,
                                           struct perturbations_workspace * ppw,
                                           double k,
                                           double * tau_ini,
                                           double tau_end);

int perturbations_switch_approximation_qs_smg(struct perturbations * ppt,
                                              struct perturbations_workspace * ppw,
                                              double tau);

int perturbations_qs_functions_at_tau_and_k_qs_smg(struct precision * ppr,
                                                   struct background * pba,
                                                   struct perturbations * ppt,
                                                   double k,
                                                   double tau,
                                                   double *mass2,
                                                   double *mass2_p,
                                                   double *rad2,
                                                   double *friction,
                                                   double *slope,
                                                   short *approx);

int perturbations_verbose_qs_smg(struct perturbations_workspace * ppw,
                                 double k,
                                 double tau_switch,
                                 int * ap_ini,
                                 int * ap_end);

int perturbations_vector_init_qs_smg(struct perturbations_workspace * ppw,
                                     struct perturbations_vector * ppv,
                                     int * pa_old);

int get_x_x_prime_qs_smg(struct precision * ppr,
                         struct background * pba,
                         struct perturbations * ppt,
                         struct perturbations_workspace * ppw,
                         double k,
                         double * x_qs_smg,
                         double * x_prime_qs_smg);

int sample_approximation_qs_smg(struct precision * ppr,
                                struct background * pba,
                                struct perturbations * ppt,
                                double k,
                                double tau_ini,
                                double tau_end,
                                double * tau_sample,
                                double * slope_sample,
                                int * approx_sample,
                                int *size_sample);

int shorten_first_qs_smg(double * tau_sample,
                         double * slope_sample,
                         int * approx_sample,
                         int size_sample,
                         double * tau_array,
                         double * slope_array,
                         int * approx_array,
                         int *size_array,
                         double tau_end);

int correct_with_slope_qs_smg(struct precision * ppr,
                              struct background * pba,
                              struct perturbations * ppt,
                              double tau_ini,
                              double tau_end,
                              double * tau_array,
                              double * slope_array,
                              int * approx_array,
                              int size_array);

int shorten_second_qs_smg(double * tau_array,
                          int * approx_array,
                          int size_array,
                          double * tau_scheme,
                          int * approx_scheme,
                          int *size_scheme);

int fit_real_scheme_qs_smg(double tau_end,
                           int * approx_scheme,
                           double * tau_scheme,
                           int size_scheme,
                           double * tau_export);

int perturbations_adiabatic_ic_smg(struct precision * ppr,
                                   struct background * pba,
                                   struct perturbations * ppt,
                                   struct perturbations_workspace * ppw,
                                   double * ptr_eta,
                                   double * ptr_delta_ur,
                                   double * ptr_theta_ur,
                                   double * ptr_shear_ur,
                                   double * ptr_l3_ur,
                                   double * ptr_delta_dr,
                                   double tau,
                                   double k,
                                   double fracnu,
                                   double om,
                                   double rho_r);

int perturbations_isocurvature_cdm_ic_smg(struct precision * ppr,
                                          struct background * pba,
                                          struct perturbations * ppt,
                                          struct perturbations_workspace * ppw,
                                          double tau,
                                          double k,
                                          double fraccdm,
                                          double om);

int perturbations_isocurvature_b_ic_smg(struct precision * ppr,
                                        struct background * pba,
                                        struct perturbations * ppt,
                                        struct perturbations_workspace * ppw,
                                        double tau,
                                        double k,
                                        double fracb,
                                        double om);

int perturbations_isocurvature_urd_ic_smg(struct precision * ppr,
                                          struct background * pba,
                                          struct perturbations * ppt,
                                          struct perturbations_workspace * ppw,
                                          double tau,
                                          double k,
                                          double fracnu,
                                          double fracg,
                                          double fracb,
                                          double om);

int perturbations_isocurvature_urv_ic_smg(struct precision * ppr,
                                          struct background * pba,
                                          struct perturbations * ppt,
                                          struct perturbations_workspace * ppw,
                                          double tau,
                                          double k,
                                          double fracnu,
                                          double fracg,
                                          double fracb,
                                          double om);

int test_ini_grav_ic_smg(struct precision * ppr,
                         struct background * pba,
                         struct perturbations * ppt);

int test_ini_extfld_ic_smg(struct precision * ppr,
                           struct background * pba,
                           struct perturbations * ppt);

int calc_extfld_ampl_smg(int n,  double kin, double bra, double dbra,
                         double run, double ten, double DelM2, double Omx,
                         double wx, double l1, double l2, double l3,
                         double l4, double l5, double l6,double l7, double l8,
                         double cs2num, double Dd, double ic_regulator_smg,
                         double * amplitude);


#ifdef __cplusplus
}
#endif

#endif
