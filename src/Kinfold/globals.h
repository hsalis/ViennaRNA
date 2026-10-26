/*
  Last changed Time-stamp: <2006-10-03 10:53:27 xtof>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: globals.h,v 1.3 2006/10/04 12:45:13 xtof Exp $
*/

#ifndef GLOBDEFS_H
#define GLOBDEFS_H

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include "config.h"

#if HAVE_LIBRNA_API3
#include <ViennaRNA/model.h>
#include <ViennaRNA/data_structures.h>
#else
#include <fold_vars.h>
#include <params.h>
#endif

typedef struct _GlobVars {
  int len;
  int full_len;
  int tx_len;
  int trajectory_index;
  int RNAP_exit_channel_length;
  int channel_left;
  int channel_length;
  int num;
  int maxS;
  int steps;
  float cut;
  float Temp;
  float startE;
  float stopE;
  float currE;
  double grow;
  int    glen;
  double time;
  double time_step;
  double phi;
  double simTime;
  double next_sample_time;
  double last_sample_time;
  double transcription_elongation_rate_nt_per_s;
  double seconds_per_internal_time_unit;
  double elongation_propensity_internal;
  double elongation_rate_current_nt_per_s;
  double elongation_propensity_internal_current;
  int max_bubble_width;
  int bubble_left;
  int bubble_width;
  int hybrid_left;
  int hybrid_len;
} GlobVars;

typedef struct _GlobArrays {
  char *ParamFile;
  char *ProgramName;
  char *BaseName;      /* output file basename */
  char *CSVOutputFile; /* optional CSV trajectory export */
  char *farbe;         /* sequence */
  char *farbe_full;    /* full sequence (for chain growth simulation) */
  char *template_dna_full; /* implicit DNA template */
  char *startform;     /* start structure */
  char **stopform;     /* stop structure(s) */
  char *currform;      /* current structure */
  char *prevform;      /* current structure of previous time step */
  float *sE;           /* energy(s) of stop structure(s) */
  double phi_bounds[3];   /* phi_min, phi_inc, phi_max */
  unsigned short subi[3]; /* seeds for random-number-generator */

#if HAVE_LIBRNA_API3
  vrna_md_t md;
  vrna_fold_compound_t *vc;
#else
  model_detailsT  md;
  paramT *params;     /* pointer to ViennaRNA energy parameters */
#endif

} GlobArrays;

typedef struct _GlobToogles {
  int Par;
  int seed;
  int dangle;
  int logML;
  int noLP;
  int noShift;
  int start;
  int stop;
  int silent;
  int phi;
  int lmin;
  int fpt;
  int rect;
  int mc;
  int verbose;
  int transcription;
  int max_bubble;
  int structure_time_only;
  int dump_transcription_neighbors;
} GlobToggles;

typedef struct {
  double terminal_hairpin_energy;
  double terminal_hairpin_dg_dt;
  double terminal_hairpin_plus1_energy;
  double rna_dna_minus1_ddg;
  double plus1_ddg;
  int    exit_channel_bp_count;
  double terminal_energy;
  double folding_inside_channel_energy;
  double elongation_rate_current;
  double rna_dna_1nt_shift_energy;
  double invasion_1nt_energy;
  double invasion_2nt_energy;
  double rloop_1nt_energy;
  double rloop_2nt_energy;
} KinfoldInformationalMetrics;

typedef struct {
  int initialized;
  double min_rna_energy, max_rna_energy;
  double min_hybrid_energy, max_hybrid_energy;
  double min_dna_duplex_energy, max_dna_duplex_energy;
  double min_total_energy, max_total_energy;
  double min_terminal_hairpin_energy, max_terminal_hairpin_energy;
  double min_terminal_hairpin_dg_dt, max_terminal_hairpin_dg_dt;
  double min_terminal_hairpin_plus1_energy, max_terminal_hairpin_plus1_energy;
  double min_rna_dna_minus1_ddg, max_rna_dna_minus1_ddg;
  double min_plus1_ddg, max_plus1_ddg;
  double min_exit_channel_bp_count, max_exit_channel_bp_count;
  double min_terminal_energy, max_terminal_energy;
  double min_folding_inside_channel_energy, max_folding_inside_channel_energy;
  double min_elongation_rate_current, max_elongation_rate_current;
  double min_rna_dna_1nt_shift_energy, max_rna_dna_1nt_shift_energy;
  double min_invasion_1nt_energy, max_invasion_1nt_energy;
  double min_invasion_2nt_energy, max_invasion_2nt_energy;
  double min_rloop_1nt_energy, max_rloop_1nt_energy;
  double min_rloop_2nt_energy, max_rloop_2nt_energy;
} KinfoldIntervalMetrics;

typedef struct {
  double last_time;
  double smoothed_terminal_hairpin;
  int    initialized;
} KinfoldMetricSmootherState;

void decode_switches(int argc, char *argv[]);
void clean_up_globals(void);
void log_prog_params(FILE *FP);
void log_start_stop(FILE *FP);
void kinfold_finalize_transcription_input(void);
int  kinfold_transcription_enabled(void);
int  kinfold_bubble_enabled(void);
const char *kinfold_full_sequence_state(void);
const char *kinfold_full_structure_state(void);
const char *kinfold_cli_structure_state(void);
const char *kinfold_state_cache_key(void);
double kinfold_current_hybrid_energy(void);
int kinfold_current_hybrid_energy_dcal(void);
double kinfold_current_dna_duplex_energy(void);
int kinfold_current_dna_duplex_energy_dcal(void);
double kinfold_current_total_energy(void);
double kinfold_current_elongation_rate(void);
double kinfold_current_elongation_propensity(void);
int kinfold_channel_left(void);
int kinfold_channel_length(void);
void kinfold_reset_sample_metrics(void);
void kinfold_prepare_sample_metrics(double sample_time);
void kinfold_get_sample_metrics(KinfoldInformationalMetrics *metrics);
void kinfold_metric_smoother_reset(KinfoldMetricSmootherState *state);
void kinfold_compute_metrics_for_time(double time_value,
                                      KinfoldMetricSmootherState *state,
                                      KinfoldInformationalMetrics *metrics);
int kinfold_stop_matches(const char *stopform);

extern GlobVars GSV;
extern GlobArrays GAV;
extern GlobToggles GTV;

#endif


/* End of file */
