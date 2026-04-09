/*
  Reusable Kinfold execution and capture API
*/

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>

#if HAVE_LIBRNA_API3
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/string_utils.h>
#else
#include <fold_vars.h>
#include <fold.h>
#include <utils.h>
#endif

#include "baum.h"
#include "nachbar.h"
#include "cache_util.h"
#include "globals.h"
#include "kinfold_api.h"

extern void  read_parameter_file(const char fname[]);
extern void  get_from_cache(cache_entry *c);

static void   ini_energy_model(void);
static void   read_data(FILE *input);
static void   clean_up(void);
static void   capture_step_at_time(kinfold_trajectory_t *trajectory, char move, double time_value);
static void   ensure_trajectory_capacity(kinfold_result_t *result, int need);
static void   ensure_step_capacity(kinfold_trajectory_t *trajectory, int need);
static char **duplicate_stop_structures(int count);

static kinfold_result_t      *active_result        = NULL;
static kinfold_trajectory_t  *active_trajectory    = NULL;
static int                   active_suppress_output = 0;


void
kinfold_result_init(kinfold_result_t *result)
{
  if (!result)
    return;

  memset(result, 0, sizeof(*result));
}


void
kinfold_result_free(kinfold_result_t *result)
{
  int i, j;

  if (!result)
    return;

  free(result->sequence);
  free(result->start_structure);

  if (result->stop_structures) {
    for (i = 0; i < result->stop_count; i++)
      free(result->stop_structures[i]);
    free(result->stop_structures);
  }

  if (result->trajectories) {
    for (i = 0; i < result->trajectory_count; i++) {
      kinfold_trajectory_t *trajectory = &(result->trajectories[i]);
      free(trajectory->termination);
      if (trajectory->steps) {
        for (j = 0; j < trajectory->step_count; j++) {
          free(trajectory->steps[j].structure);
          free(trajectory->steps[j].sequence_state);
        }
        free(trajectory->steps);
      }
    }
    free(result->trajectories);
  }

  memset(result, 0, sizeof(*result));
}


int
kinfold_capture_suppress_output(void)
{
  return active_suppress_output;
}


void
kinfold_capture_prepare_result(void)
{
  if (!active_result)
    return;

  free(active_result->sequence);
  active_result->sequence = strdup(GAV.farbe_full ? GAV.farbe_full : GAV.farbe);

  free(active_result->start_structure);
  active_result->start_structure = strdup(GAV.startform ? GAV.startform : "");

  if (active_result->stop_structures) {
    int i;
    for (i = 0; i < active_result->stop_count; i++)
      free(active_result->stop_structures[i]);
    free(active_result->stop_structures);
  }

  active_result->stop_count = GSV.maxS;
  active_result->stop_structures = duplicate_stop_structures(GSV.maxS);
}


void
kinfold_capture_begin_trajectory(int trajectory_index)
{
  kinfold_trajectory_t *trajectory;

  if (!active_result)
    return;

  ensure_trajectory_capacity(active_result, active_result->trajectory_count + 1);

  trajectory = &(active_result->trajectories[active_result->trajectory_count++]);
  memset(trajectory, 0, sizeof(*trajectory));

  trajectory->trajectory_index = trajectory_index;
  trajectory->seed_initial[0] = GAV.subi[0];
  trajectory->seed_initial[1] = GAV.subi[1];
  trajectory->seed_initial[2] = GAV.subi[2];

  active_trajectory = trajectory;
}


void
kinfold_capture_initial_state(void)
{
  if (!active_trajectory)
    return;

  capture_step_at_time(active_trajectory, '\0', 0.0);
}


void
kinfold_capture_transition(char move)
{
  if (!active_trajectory)
    return;

  capture_step_at_time(active_trajectory, move, GSV.simTime);
}


void
kinfold_capture_sample(double time_value)
{
  if (!active_trajectory)
    return;

  capture_step_at_time(active_trajectory, '\0', time_value);
}


void
kinfold_capture_terminal(int found_stop)
{
  char marker[16];

  if (!active_trajectory)
    return;

  if ((active_trajectory->step_count == 0) ||
      fabs(active_trajectory->steps[active_trajectory->step_count - 1].time - GSV.simTime) > 1e-9)
    capture_step_at_time(active_trajectory, '\0', GSV.simTime);

  active_trajectory->seed_final[0] = GAV.subi[0];
  active_trajectory->seed_final[1] = GAV.subi[1];
  active_trajectory->seed_final[2] = GAV.subi[2];
  active_trajectory->stop_index = found_stop;

  free(active_trajectory->termination);
  if (found_stop > 0) {
    snprintf(marker, sizeof(marker), "X%d", found_stop);
    active_trajectory->termination = strdup(marker);
  } else {
    active_trajectory->termination = strdup("O");
  }
}


int
kinfold_run_with_args(int               argc,
                      char              **argv,
                      FILE              *input,
                      kinfold_result_t  *result,
                      int               suppress_output)
{
  int           i;
  int           rect;
  char          *start, *tmp;

  active_result = result;
  active_trajectory = NULL;
  active_suppress_output = suppress_output ? 1 : 0;

  decode_switches(argc, argv);
  ini_energy_model();
  read_data(input);

#if HAVE_LIBRNA_API3
  tmp     = vrna_cut_point_insert(GAV.farbe_full ? GAV.farbe_full : GAV.farbe, cut_point);
  GAV.vc  = vrna_fold_compound(tmp, &(GAV.md), VRNA_OPTION_DEFAULT);
  free(tmp);
  GAV.vc->length = GSV.len;
#endif

  initialize_cache();

  if (active_result)
    kinfold_capture_prepare_result();

  start = strdup(GAV.startform);
  rect = GTV.rect;

  for (i = 0; i < GSV.num; i++) {
    GTV.rect = rect;

    ini_or_reset_rl();
    if (!kinfold_transcription_enabled() && GSV.grow > 0) {
      if ((int)strlen(GAV.farbe) > GSV.glen) {
        start[GSV.glen] = '\0';
        GAV.farbe[GSV.glen] = '\0';
        strcpy(GAV.startform, start);
        strcpy(GAV.currform, start);
        GSV.len = GSV.glen;
#if HAVE_LIBRNA_API3
        GAV.vc->length = GSV.len;
#endif
      }
      clean_up_rl();
      ini_or_reset_rl();
    }

    GSV.simTime = 0.0;
    GSV.last_sample_time = 0.0;
    GSV.next_sample_time = GSV.time_step;
    kinfold_reset_sample_metrics();

    if (active_result) {
      kinfold_capture_begin_trajectory(i);
      kinfold_capture_initial_state();
    }

    for (GSV.steps = 1;; GSV.steps++) {
      cache_entry *c;

      if (GTV.dump_transcription_neighbors && kinfold_transcription_enabled()) {
        move_it();
      } else if ((c = lookup_cache((char *)kinfold_state_cache_key()))) {
        get_from_cache(c);
      } else {
        move_it();
      }

      if (sel_nb() > 0)
        break;
    }
  }

  free(start);
  clean_up();
  active_result = NULL;
  active_trajectory = NULL;
  active_suppress_output = 0;

  return 0;
}


static void
ensure_trajectory_capacity(kinfold_result_t *result,
                           int              need)
{
  int new_capacity;

  if (need <= result->trajectory_capacity)
    return;

  new_capacity = (result->trajectory_capacity > 0) ? result->trajectory_capacity : 4;
  while (new_capacity < need)
    new_capacity *= 2;

  result->trajectories = (kinfold_trajectory_t *)realloc(result->trajectories,
                                                          sizeof(kinfold_trajectory_t) *
                                                          (size_t)new_capacity);
  assert(result->trajectories != NULL);
  memset(result->trajectories + result->trajectory_capacity,
         0,
         sizeof(kinfold_trajectory_t) * (size_t)(new_capacity - result->trajectory_capacity));
  result->trajectory_capacity = new_capacity;
}


static void
ensure_step_capacity(kinfold_trajectory_t *trajectory,
                     int                  need)
{
  int new_capacity;

  if (need <= trajectory->step_capacity)
    return;

  new_capacity = (trajectory->step_capacity > 0) ? trajectory->step_capacity : 64;
  while (new_capacity < need)
    new_capacity *= 2;

  trajectory->steps = (kinfold_step_t *)realloc(trajectory->steps,
                                                sizeof(kinfold_step_t) *
                                                (size_t)new_capacity);
  assert(trajectory->steps != NULL);
  memset(trajectory->steps + trajectory->step_capacity,
         0,
         sizeof(kinfold_step_t) * (size_t)(new_capacity - trajectory->step_capacity));
  trajectory->step_capacity = new_capacity;
}


static void
capture_step_at_time(kinfold_trajectory_t *trajectory,
                     char                 move,
                     double               time_value)
{
  kinfold_step_t *step;
  KinfoldInformationalMetrics metrics;

  kinfold_prepare_sample_metrics(time_value);
  ensure_step_capacity(trajectory, trajectory->step_count + 1);
  step = &(trajectory->steps[trajectory->step_count++]);
  kinfold_get_sample_metrics(&metrics);

  step->time              = time_value;
  step->transcribed_length = GSV.tx_len;
  step->energy            = GSV.currE;
  step->total_energy      = kinfold_current_total_energy();
  step->dna_duplex_energy = kinfold_current_dna_duplex_energy();
  step->terminal_hairpin_energy = metrics.terminal_hairpin_energy;
  step->terminal_hairpin_dg_dt = metrics.terminal_hairpin_dg_dt;
  step->terminal_hairpin_plus1_energy = metrics.terminal_hairpin_plus1_energy;
  step->rna_dna_minus1_ddg = metrics.rna_dna_minus1_ddg;
  step->plus1_ddg = metrics.plus1_ddg;
  step->bubble_length     = GSV.bubble_width;
  step->bubble_left       = GSV.bubble_left;
  step->hybrid_length     = GSV.hybrid_len;
  step->hybrid_left       = GSV.hybrid_left;
  step->move              = move;
  step->structure         = strdup(kinfold_full_structure_state());
  step->sequence_state    = strdup(kinfold_full_sequence_state());
}


static char **
duplicate_stop_structures(int count)
{
  char  **copy;
  int   i;

  if (count <= 0)
    return NULL;

  copy = (char **)calloc((size_t)count, sizeof(char *));
  assert(copy != NULL);
  for (i = 0; i < count; i++)
    copy[i] = strdup(GAV.stopform[i] ? GAV.stopform[i] : "");

  return copy;
}


static void
ini_energy_model(void)
{
#if HAVE_LIBRNA_API3
  vrna_md_set_default(&(GAV.md));
#else
  set_model_details(&(GAV.md));
#endif

  if (!GTV.seed) {
    init_rand();
    GAV.subi[0] = xsubi[0];
    GAV.subi[1] = xsubi[1];
    GAV.subi[2] = xsubi[2];
  } else {
    xsubi[0] = GAV.subi[0];
    xsubi[1] = GAV.subi[1];
    xsubi[2] = GAV.subi[2];
  }

  GAV.md.logML        = logML = GTV.logML;
  GAV.md.dangles      = dangles = GTV.dangle;
  GAV.md.temperature  = temperature = GSV.Temp;
  GAV.md.noLP         = GTV.noLP;

  if (GTV.Par)
    read_parameter_file(GAV.ParamFile);

#if HAVE_LIBRNA_API2
  GAV.params = get_scaled_parameters(temperature, GAV.md);
#endif
}


static void
read_data(FILE *input)
{
  char  *ctmp, *c, **s;
  int   i, len;

  ctmp = get_line(input);
  len = strlen(ctmp);
  if ((c = strchr(ctmp, '&'))) {
    cut_point = (int)(c - ctmp) + 1;
    for (; *c; c++)
      *c = *(c + 1);
  }

  if (kinfold_transcription_enabled() && cut_point > 0) {
    fprintf(stderr,
            "Transcription-aware mode does not support cofold / cut-point inputs in this tranche\n");
    exit(EXIT_FAILURE);
  }

  GSV.full_len = len;
  GSV.tx_len   = kinfold_transcription_enabled() ? GSV.glen : len;
  if (GSV.tx_len < 0)
    GSV.tx_len = 0;
  if (GSV.tx_len > GSV.full_len)
    GSV.tx_len = GSV.full_len;

  GAV.farbe = (char *)calloc((size_t)len + 1, sizeof(char));
  assert(GAV.farbe != NULL);
  sscanf(ctmp, "%s", GAV.farbe);
  GSV.len = GSV.tx_len;
  for (i = 0; i < len; i++)
    GAV.farbe[i] = (char)toupper((unsigned char)GAV.farbe[i]);
  GAV.farbe_full = strdup(GAV.farbe);
  GAV.farbe[GSV.len] = '\0';
  free(ctmp);

  GAV.currform = (char *)calloc((size_t)GSV.full_len + 1, sizeof(char));
  assert(GAV.currform != NULL);
  GAV.prevform = (char *)calloc((size_t)GSV.full_len + 1, sizeof(char));
  assert(GAV.prevform != NULL);
  GAV.startform = (char *)calloc((size_t)GSV.full_len + 1, sizeof(char));
  assert(GAV.startform != NULL);

  if (GTV.start) {
    ctmp = get_line(input);
    sscanf(ctmp, "%s", GAV.startform);
    if ((c = strchr(GAV.startform, '&'))) {
      for (; *c; c++)
        *c = *(c + 1);
    }

    if (!kinfold_transcription_enabled() && ((int)strlen(GAV.startform) != GSV.len)) {
      fprintf(stderr,
              "read_data():\n%s\n%s\n unequal length!\n",
              GAV.farbe,
              GAV.startform);
      exit(EXIT_FAILURE);
    }
    free(ctmp);
  } else {
    for (i = 0; i < GSV.len; i++)
      GAV.startform[i] = '.';
  }

  if (GTV.stop) {
    s = GAV.stopform;
    while ((ctmp = get_line(input))) {
      *s = (char *)calloc((size_t)GSV.full_len + 1, sizeof(char));
      sscanf(ctmp, "%s", *s);
      if ((c = strchr(ctmp, '&'))) {
        for (; *c; c++)
          *c = *(c + 1);
      }

      if ((s - GAV.stopform) >= GSV.maxS) {
        fprintf(stderr,
                "WARNING: Can handle a maximum of %d stop structures\n",
                GSV.maxS);
        break;
      }

      if ((!kinfold_transcription_enabled() && ((int)strlen(*s) != GSV.len)) ||
          (kinfold_transcription_enabled() && ((int)strlen(*s) != GSV.full_len))) {
        fprintf(stderr,
                "read_data():\n%s\n%s\n unequal length!\n",
                GAV.farbe,
                *s);
        exit(EXIT_FAILURE);
      }

      s++;
      free(ctmp);
    }
    GSV.maxS = (int)(s - GAV.stopform);
  } else {
    GSV.maxS = 1;
    GAV.stopform[0] = (char *)calloc((size_t)GSV.full_len + 1, sizeof(char));
  }

  if (kinfold_transcription_enabled()) {
    free(GAV.farbe);
    GAV.farbe = (char *)calloc((size_t)GSV.full_len + 1, sizeof(char));
    assert(GAV.farbe != NULL);
    memcpy(GAV.farbe, GAV.farbe_full, (size_t)GSV.tx_len);
    GAV.farbe[GSV.tx_len] = '\0';
    kinfold_finalize_transcription_input();
    if (kinfold_bubble_enabled()) {
      GSV.bubble_width = (GSV.tx_len < GSV.max_bubble_width) ? GSV.tx_len : GSV.max_bubble_width;
      GSV.bubble_left  = GSV.tx_len - GSV.bubble_width;
      GSV.hybrid_len   = GSV.bubble_width;
      GSV.hybrid_left  = GSV.bubble_left;
    }
  }

  GAV.sE = (float *)calloc((size_t)GSV.maxS, sizeof(float));
}


static void
clean_up(void)
{
  clean_up_globals();
  clean_up_rl();
  clean_up_nbList();
  kill_cache();
}
