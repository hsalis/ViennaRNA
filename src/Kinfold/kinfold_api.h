/*
  Reusable Kinfold execution and capture API
*/

#ifndef KINFOLD_API_H
#define KINFOLD_API_H

#include <stdio.h>

typedef struct {
  double  time;
  int     transcribed_length;
  double  energy;
  double  total_energy;
  double  dna_duplex_energy;
  double  terminal_hairpin_energy;
  double  terminal_hairpin_dg_dt;
  double  terminal_hairpin_plus1_energy;
  double  rna_dna_minus1_ddg;
  double  plus1_ddg;
  int     exit_channel_left;
  int     exit_channel_length;
  int     exit_channel_bp_count;
  double  terminal_energy;
  double  folding_inside_channel_energy;
  double  elongation_rate_current;
  double  rna_dna_1nt_shift_energy;
  double  invasion_1nt_energy;
  double  invasion_2nt_energy;
  double  rloop_1nt_energy;
  double  rloop_2nt_energy;
  double  min_rna_energy, max_rna_energy;
  double  min_hybrid_energy, max_hybrid_energy;
  double  min_dna_duplex_energy, max_dna_duplex_energy;
  double  min_total_energy, max_total_energy;
  double  min_terminal_hairpin_energy, max_terminal_hairpin_energy;
  double  min_terminal_hairpin_dg_dt, max_terminal_hairpin_dg_dt;
  double  min_terminal_hairpin_plus1_energy, max_terminal_hairpin_plus1_energy;
  double  min_rna_dna_minus1_ddg, max_rna_dna_minus1_ddg;
  double  min_plus1_ddg, max_plus1_ddg;
  double  min_exit_channel_bp_count, max_exit_channel_bp_count;
  double  min_terminal_energy, max_terminal_energy;
  double  min_folding_inside_channel_energy, max_folding_inside_channel_energy;
  double  min_elongation_rate_current, max_elongation_rate_current;
  double  min_rna_dna_1nt_shift_energy, max_rna_dna_1nt_shift_energy;
  double  min_invasion_1nt_energy, max_invasion_1nt_energy;
  double  min_invasion_2nt_energy, max_invasion_2nt_energy;
  double  min_rloop_1nt_energy, max_rloop_1nt_energy;
  double  min_rloop_2nt_energy, max_rloop_2nt_energy;
  int     bubble_length;
  int     bubble_left;
  int     hybrid_length;
  int     hybrid_left;
  char    move;
  char    *structure;
  char    *sequence_state;
} kinfold_step_t;

typedef struct {
  int             trajectory_index;
  unsigned short  seed_initial[3];
  unsigned short  seed_final[3];
  char            *termination;
  int             stop_index;
  int             step_count;
  int             step_capacity;
  kinfold_step_t  *steps;
} kinfold_trajectory_t;

typedef struct {
  char                  *sequence;
  char                  *start_structure;
  char                  **stop_structures;
  int                   stop_count;
  int                   trajectory_count;
  int                   trajectory_capacity;
  kinfold_trajectory_t  *trajectories;
} kinfold_result_t;

void
kinfold_result_init(kinfold_result_t *result);

void
kinfold_result_free(kinfold_result_t *result);

int
kinfold_run_with_args(int               argc,
                      char              **argv,
                      FILE              *input,
                      kinfold_result_t  *result,
                      int               suppress_output);

void
kinfold_capture_prepare_result(void);

void
kinfold_capture_begin_trajectory(int trajectory_index);

void
kinfold_capture_initial_state(void);

void
kinfold_capture_transition(char move);

void
kinfold_capture_sample(double time_value);

void
kinfold_capture_terminal(int found_stop);

int
kinfold_capture_suppress_output(void);

#endif
