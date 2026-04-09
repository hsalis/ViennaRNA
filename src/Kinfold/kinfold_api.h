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
