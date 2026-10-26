/*
  Last changed Time-stamp: <2011-03-10 18:02:10 ivo>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: nachbar.c,v 1.8 2008/06/03 21:55:11 ivo Exp $
*/

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "assert.h"

#if HAVE_LIBRNA_API3
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/utils.h>
#else
#include <fold_vars.h>
#include <energy_const.h>
#include <utils.h>
#endif

#include "cache_util.h"
#include "baum.h"
#include "nachbar.h"
#include "kinfold_api.h"

static char UNUSED rcsid[]="$Id: nachbar.c,v 1.8 2008/06/03 21:55:11 ivo Exp $";

/* arrays */
static short *neighbor_list=NULL;
static float *bmf=NULL; /* boltzmann weight of structure */
static unsigned char *move_kinds=NULL;
static const char *costring(const char *str);

/* globals for laplace stuff */
static double L = 0.0;
static double D = 0.0;
static double sumT = 0.0;
static double sumK = 0.0;
static double sumKK = 0.0;
static double sumD = 0.0;
static double *energies=NULL; /* energies of neighbors */

/* variables */
/*  static double highestE = -1000.0; */
/*  static double OhighestE = -1000.0; */
/*  static char *highestS, *OhighestS; */
static int lmin = 1;
static int top = 0;
static int is_from_cache = 0;
/*  static double meanE = 0.0; */
static double totalflux = 0.0;
static double Zeit = 0.0;
static double zeitInc = 0.0;
static double _RT = 0.6;

enum {
  MOVE_RNA = 0,
  MOVE_ELONGATE,
  MOVE_HYBRID_FORM,
  MOVE_HYBRID_BREAK,
  MOVE_BUBBLE_GROW,
  MOVE_BUBBLE_INVASION
};

/* public functiones */
void ini_nbList(int chords);
void update_nbList(int i, int j, int iE);
static void update_nbList_kind(unsigned char kind, int i, int j, int rna_iE, double total_E);
int sel_nb(void);
void clean_up_nbList(void);
extern void update_tree(int i, int j);
char kinfold_move_char(int next);

/* privat functiones */
static void reset_nbList(void);
static void grow_chain(void);
static void append_transcription_neighbors(void);
static int transcription_active(void);
static double displayed_energy(void);
static double displayed_structure_energy(void);
static double displayed_hybrid_energy(void);
static double displayed_dna_duplex_energy(void);
static int should_emit_to_stdout(void);
static void emit_state_line(double out_time, int found_stop, int terminal);
static void emit_sampled_states(double limit_time, int include_limit_time);
static void print_console_legend(void);
static void csv_write_escaped(FILE *fp, const char *text);
static void csv_write_header(FILE *fp);
static void emit_csv_state_line(double out_time, int found_stop, int terminal, const KinfoldInformationalMetrics *metrics);
static void interval_metrics_reset(KinfoldIntervalMetrics *interval);
static void interval_metrics_update(KinfoldIntervalMetrics *interval, const KinfoldInformationalMetrics *metrics);
static void interval_seed_from_emitted_state(const KinfoldInformationalMetrics *metrics);
static void track_hidden_state(double time_value);
static void begin_trajectory_if_needed(void);
static double compute_transition_rate(unsigned char kind, double curr_total, double total_E);
static int position_in_exit_channel_current(int idx0);
static void dump_transcription_state(FILE *fp);
static void dump_transcription_neighbor(FILE *fp,
                                       unsigned char kind,
                                       int i,
                                       int j,
                                       int next_tx_len,
                                       int next_bubble_left,
                                       int next_bubble_width,
                                       int next_hybrid_left,
                                       int next_hybrid_len,
                                       int rna_iE,
                                       double total_E,
                                       double rate);
static int eval_extended_structure_dcal(int new_len, const char *structure);
static void set_active_prefix_length(int new_len);
static void apply_transcription_move(int next);
static FILE *logFP=NULL;
static FILE *csvFP=NULL;
static KinfoldIntervalMetrics current_interval_metrics;
static KinfoldIntervalMetrics last_exported_interval_metrics;
static KinfoldMetricSmootherState csv_metric_state;
static int csv_tracker_active = 0;
static int csv_tracker_trajectory = -1;
static double csv_last_emitted_time = -1.0;
static int structure_time_header_printed = 0;
static int console_legend_printed = 0;

/**/
void ini_nbList(int chords) {
  char logFN[256];

  _RT = (((temperature + K0) * GASCONST) / 1000.0);
  structure_time_header_printed = 0;
  console_legend_printed = 0;
  if (neighbor_list!=NULL) return;
  /*
    list for move coding
    make room for 2*chords neighbors (safe bet)
  */
  if (chords == 0) chords = 1;
  neighbor_list = (short *)calloc(4*chords, sizeof(short));
  assert(neighbor_list != NULL);
  /*
    list for Boltzmann-factors
  */
  bmf = (float *)calloc(2*chords, sizeof(double));
  assert(bmf != NULL);

  /* list of neighbor energies */
  energies = (double*)calloc(2*chords, sizeof(double));
  assert(energies != NULL);
  move_kinds = (unsigned char *)calloc(2*chords, sizeof(unsigned char));
  assert(move_kinds != NULL);
  
  if (!kinfold_capture_suppress_output()) {
    /* open log-file */
    logFP = fopen(strcat(strcpy(logFN, GAV.BaseName), ".log"), "a+");
    assert(logFP != NULL);

    /* log initial condition */
    log_prog_params(logFP);
    log_start_stop(logFP);
  } else {
    logFP = NULL;
  }

  if (GAV.CSVOutputFile && GAV.CSVOutputFile[0] != '\0') {
    csvFP = fopen(GAV.CSVOutputFile, "w");
    if (csvFP == NULL) {
      fprintf(stderr, "Failed to open CSV output file '%s'\n", GAV.CSVOutputFile);
      exit(EXIT_FAILURE);
    }
    csv_write_header(csvFP);
    fflush(csvFP);
  } else {
    csvFP = NULL;
  }
}

/**/
void update_nbList(int i, int j, int iE) {
  double total_E;

  total_E = (double)iE / 100.;
  if (kinfold_bubble_enabled())
    total_E += kinfold_current_hybrid_energy();

  update_nbList_kind(MOVE_RNA, i, j, iE, total_E);
}

static void update_nbList_kind(unsigned char kind, int i, int j, int rna_iE, double total_E) {
  double rna_E, p, curr_total;

  rna_E = (double)rna_iE / 100.;
  neighbor_list[2*top] = (short)i;
  neighbor_list[2*top+1] = (short)j;
  move_kinds[top] = kind;

  curr_total = transcription_active() ? kinfold_current_total_energy() : GSV.currE;

  energies[top] = rna_E;
  L += curr_total - total_E;
  D++;

  p = compute_transition_rate(kind, curr_total, total_E);

  totalflux += p;
  bmf[top++] = (float)p;
  if ((total_E - curr_total) < 0) lmin = 0;
  if (((total_E - curr_total) == 0) && (lmin==1)) lmin = 2;
}

/**/
void get_from_cache(cache_entry *c) {
  top = c->top;
  totalflux = c->flux;
  GSV.currE = c->energy;
  lmin = c->lmin;
  memcpy(neighbor_list, c->neighbors, 2*top*sizeof(short));
  memcpy(bmf, c->rates, top*sizeof(float));
  memcpy(energies, c->energies, top*sizeof(double));
  memcpy(move_kinds, c->move_kinds, top*sizeof(unsigned char));
  is_from_cache = 1;
}

/**/
void put_in_cache(void) {
  cache_entry *c;

  if ((c = (cache_entry *) malloc(sizeof(cache_entry)))==NULL) {
    fprintf(stderr, "out of memory\n"); exit(255);
  }
  c->structure = strdup(kinfold_state_cache_key());
  c->neighbors = (short *) malloc(top*2*sizeof(short));
  memcpy(c->neighbors,neighbor_list,top*2*sizeof(short));
  c->rates = (float *) malloc(top*sizeof(float));
  memcpy(c->rates, bmf, top*sizeof(float));
  c->energies = (double*)malloc(top*sizeof(double));
  memcpy(c->energies, energies, top*sizeof(double));
  c->move_kinds = (unsigned char *)malloc(top*sizeof(unsigned char));
  memcpy(c->move_kinds, move_kinds, top*sizeof(unsigned char));
  c->top = top;
  c->lmin = lmin;
  c->flux = totalflux;
  c->energy = GSV.currE;
  write_cache(c);
}

/*============*/

int sel_nb(void) {

  char trans, **s;
  int next, i;
  double pegel = 0.0, schwelle = 0.0, zufall = 0.0;
  int found_stop=0;

  begin_trajectory_if_needed();

  if (!is_from_cache && transcription_active())
    append_transcription_neighbors();

  /* before we select a move, store current conformation in cache */
  /* ... unless it just came from there */
  if (GTV.dump_transcription_neighbors && transcription_active()) {
    /* keep neighbor dumps complete and state-local during regression runs */
  } else if ( !is_from_cache ) put_in_cache();
  else
    /* laplace stuff */
    for (i=0; i<top; i++) {
      L += (GSV.currE - energies[i]);
      D++;
    }
  is_from_cache = 0;

  /* draw 2 different a random number */
  schwelle = urn();
  while ( zufall==0 ) zufall = urn();

  /* advance internal clock */
  if (totalflux>0)
    zeitInc = (log(1. / zufall) / totalflux);
  else {
    if (transcription_active() && (kinfold_current_elongation_propensity() > 0.0))
      zeitInc = 1.0 / kinfold_current_elongation_propensity();
    else if (GSV.grow>0) zeitInc=GSV.grow;
    else zeitInc = GSV.time;
  }

  Zeit += zeitInc;
  GSV.simTime = Zeit;

  /* laplace stuff */
  sumK  += L*zeitInc;
  sumKK += L*L*zeitInc;
  sumD  += D*zeitInc;
  
  if (!transcription_active() && GSV.grow>0 && GSV.len < strlen(GAV.farbe_full)) grow_chain();

  /* meanE /= (double)top; */

  /* normalize boltzmann weights */
  schwelle *=totalflux;

  /* and choose a neighbour structure next */
  for (next = 0; next < top; next++) {
    pegel += bmf[next];
    if (pegel > schwelle) break;
  }

  /* in case of rounding errors */
  if (next==top) next=top-1;

  /*
    process termination contitiones
  */
  /* is current structure identical to a stop structure ?*/
  for (found_stop = 0, s = GAV.stopform; *s; s++) {
    if (kinfold_stop_matches(*s)) {
      found_stop = (s - GAV.stopform) + 1;
      break;
    }
  }

  /* Recurrence time: Ignore when you observe the start structure for the first time. */
  if ((found_stop > 0) && (GTV.rect == 1) && (strcmp(GAV.startform, GAV.currform) == 0)) {
    GTV.rect = 0; found_stop = 0;
  }

  if ( ((found_stop > 0) && (GTV.fpt == 1)) || (Zeit > GSV.time) ) {
    /* met condition to stop simulation */
    double terminal_time = (found_stop > 0) ? Zeit : GSV.time;

    /* laplace stuff */
    double K, KK, N, sigma;
    K = sumK/Zeit;
    KK = sumKK/Zeit;
    N = sumD/Zeit;
    /* graph Laplacian is - Laplace-Beltrami operator */
    sigma = -1.0*sqrt((KK-K*K)/N)/(K/N);
    
    emit_sampled_states(terminal_time, 0);
    GSV.simTime = terminal_time;
    kinfold_prepare_sample_metrics(terminal_time);
    emit_state_line(terminal_time, found_stop, 1);

    /* laplace stuff */
    if ((!GTV.silent) && (!kinfold_capture_suppress_output()) && GTV.phi)
      printf("Curvature fluctuation sigma = %7.5f\n", sigma);

    /* this goes to log */
    if (logFP != NULL) {
      fprintf(logFP, "(%5hu %5hu %5hu)", GAV.subi[0], GAV.subi[1], GAV.subi[2]);
      /* comment log steps of simulation as well !!! %6.2f  round */
      if ( found_stop ) {
        fprintf(logFP," X%02d %12.3f", found_stop, terminal_time);

        /* laplace stuff */
        if (GTV.phi) fprintf(logFP, " %3g %7.5f", GSV.phi, sigma);

        fprintf(logFP,"\n");
      }
      else {
        fprintf(logFP," O   %12.3f", terminal_time);

        /* laplace stuff */
        if (GTV.phi) fprintf(logFP, " %3g %7.5f", GSV.phi, sigma);

        fprintf(logFP," %d %s\n", lmin,
                transcription_active() ? kinfold_full_structure_state() : costring(GAV.currform));
      }
      fflush(logFP);
    }

    /* set random number for next round */
    GAV.subi[0] = xsubi[0];
    GAV.subi[1] = xsubi[1];
    GAV.subi[2] = xsubi[2];
    kinfold_capture_terminal(found_stop);
    
    Zeit = 0.0;
    GSV.simTime = 0.0;

    /* reset laplace stuff for next trajectory */
    sumT = 0.0;
    sumK = 0.0;
    sumKK = 0.0;
    sumD = 0.0;
    L = 0.0;
    D = 0.0;
    
    /*  highestE = OhighestE = -1000.0; */
    reset_nbList();
    costring(NULL);
    return(1);
  }
  else {
    emit_sampled_states(Zeit, 1);
  }
  /* store last lmin seen, so we can avoid printing the same lmin twice */
  if (lmin==1)
    strcpy(GAV.prevform, GAV.currform);

#if 0
  if (lmin==1) {
    /* went back to previous lmin */
    if (strcmp(GAV.prevform, GAV.currform) == 0) {
      if (OhighestE < highestE) {
	highestE = OhighestE;  /* delete loop */
	strcpy(highestS, OhighestS);
      }
    } else {
      strcpy(GAV.prevform, GAV.currform);
      OhighestE = 10000.;
    }
  }

  if ( strcmp(GAV.currform, GAV.startform)==0 ) {
    OhighestE = highestE = -1000.;
    highestS[0] = 0;
  }

  /* log highes energy */
  if (GSV.currE > highestE) {
    OhighestE = highestE;
    highestE = GSV.currE;
    strcpy(OhighestS, highestS);
    strcpy(highestS, GAV.currform);
  }
#endif

  if (next>=0) {
    if (move_kinds[next] == MOVE_RNA) {
      GSV.currE = energies[next];
      update_tree(neighbor_list[2*next], neighbor_list[2*next+1]);
    } else {
      apply_transcription_move(next);
    }
    if (csv_tracker_active) {
      if (fabs(csv_last_emitted_time - Zeit) <= 1e-9)
        interval_metrics_reset(&current_interval_metrics);
      track_hidden_state(Zeit);
    }
  } else {
    clean_up_rl(); ini_or_reset_rl();
  }

  reset_nbList();
  return(0);
}

/*==========================*/
static void reset_nbList(void) {

  top = 0;
  totalflux = 0.0;
  /*    meanE = 0.0; */
  lmin = 1;
}

static int transcription_active(void) {
  return kinfold_transcription_enabled();
}

static double displayed_energy(void) {
  return transcription_active() ? kinfold_current_total_energy() : GSV.currE;
}

static double displayed_structure_energy(void) {
  return GSV.currE;
}

static double displayed_hybrid_energy(void) {
  return kinfold_bubble_enabled() ? kinfold_current_hybrid_energy() : 0.0;
}

static double displayed_dna_duplex_energy(void) {
  return kinfold_bubble_enabled() ? kinfold_current_dna_duplex_energy() : 0.0;
}

static int position_in_exit_channel_current(int idx0) {
  int left = kinfold_channel_left();
  int length = kinfold_channel_length();

  return (length > 0) && (idx0 >= left) && (idx0 < left + length);
}

static int should_emit_to_stdout(void) {
  return (!GTV.silent) &&
         (!kinfold_capture_suppress_output()) &&
         (displayed_structure_energy() <= GSV.stopE + GSV.cut);
}

static void csv_write_escaped(FILE *fp, const char *text) {
  const char *p;

  fputc('"', fp);
  if (text) {
    for (p = text; *p; p++) {
      if (*p == '"')
        fputc('"', fp);
      fputc(*p, fp);
    }
  }
  fputc('"', fp);
}

static void csv_write_header(FILE *fp) {
  fprintf(fp,
          "trajectory_index,step,time,terminal,stop_index,transcribed_length,"
          "bubble_left,bubble_width,hybrid_left,hybrid_len,exit_channel_left,exit_channel_length,exit_channel_bp_count,"
          "structure,sequence_state,"
          "dG_RNA_RNA,dG_RNA_DNA,dG_DNA_DNA,dG_total,"
          "dG_terminal_hairpin,dG_hairpin_dt,dG_terminal_hairpin_plus1,ddG_RNA_DNA_minus1,ddG_plus1,"
          "dG_terminal,dG_folding_inside_channel,elongation_rate_current,dG_RNA_DNA_1nt_shift,"
          "dG_invasion_1nt,dG_invasion_2nt,dG_Rloop_1nt,dG_Rloop_2nt,"
          "min_dG_RNA_RNA,max_dG_RNA_RNA,min_dG_RNA_DNA,max_dG_RNA_DNA,min_dG_DNA_DNA,max_dG_DNA_DNA,min_dG_total,max_dG_total,"
          "min_dG_terminal_hairpin,max_dG_terminal_hairpin,min_dG_hairpin_dt,max_dG_hairpin_dt,"
          "min_dG_terminal_hairpin_plus1,max_dG_terminal_hairpin_plus1,min_ddG_RNA_DNA_minus1,max_ddG_RNA_DNA_minus1,min_ddG_plus1,max_ddG_plus1,"
          "min_exit_channel_bp_count,max_exit_channel_bp_count,min_dG_terminal,max_dG_terminal,min_dG_folding_inside_channel,max_dG_folding_inside_channel,"
          "min_elongation_rate_current,max_elongation_rate_current,min_dG_RNA_DNA_1nt_shift,max_dG_RNA_DNA_1nt_shift,"
          "min_dG_invasion_1nt,max_dG_invasion_1nt,min_dG_invasion_2nt,max_dG_invasion_2nt,min_dG_Rloop_1nt,max_dG_Rloop_1nt,min_dG_Rloop_2nt,max_dG_Rloop_2nt\n");
}

static void print_console_legend(void) {
  if (console_legend_printed)
    return;

  console_legend_printed = 1;

  printf("# Console output legend\n");
  printf("# structure: RNA secondary structure at the sampled timepoint; x marks untranscribed nucleotides, b marks bubble boundaries, c marks the RNAP exit-channel start when present.\n");

  if (GTV.structure_time_only) {
    printf("# 1. time: simulation time for the sampled state.\n");
  } else if (kinfold_bubble_enabled()) {
    printf("# 1. dG_RNA_RNA: RNA folding free energy of the current RNA secondary structure.\n");
    printf("# 2. dG_RNA_DNA: RNA:DNA hybrid free energy for the currently paired hybrid region.\n");
    printf("# 3. dG_DNA_DNA: informational DNA:DNA duplex free energy for base-paired DNA nts in the bubble.\n");
    printf("# 4. dG_terminal_hairpin: folding free energy of the 3'-most exterior RNA paired component upstream of the bubble.\n");
    printf("# 5. dG_hairpin_dt: causal smoothed time derivative of dG_terminal_hairpin.\n");
    printf("# 6. dG_terminal_hairpin_plus1: hypothetical terminal-hairpin energy after adding one upstream hybrid RNA nt to the hairpin.\n");
    printf("# 7. ddG_RNA_DNA_minus1: change in RNA:DNA hybrid energy after removing the upstream-most hybrid nucleotide.\n");
    printf("# 8. ddG_plus1: informational net free-energy change for the plus1 hairpin/hybrid displacement state.\n");
    printf("# 9. channel_bp_count: number of RNA-RNA base pairs with at least one endpoint inside the RNAP exit channel.\n");
    printf("# 10. dG_terminal: folding free energy of the most downstream RNA paired component intersecting the exit channel.\n");
    printf("# 11. dG_folding_inside_channel: RNA folding free energy evaluated only inside the exit channel, ignoring external pairs.\n");
    printf("# 12. elongation_rate_current: sequence-dependent RNAP elongation rate at the current position.\n");
    printf("# 13. dG_RNA_DNA_1nt_shift: RNA:DNA hybrid free energy after a hypothetical 1-nt upstream RNA register shift.\n");
    printf("# 14. dG_invasion_1nt: minimum RNA folding free energy when the 3'-most hybrid RNA nt may invade RNA structure.\n");
    printf("# 15. dG_invasion_2nt: minimum RNA folding free energy when the two 3'-most hybrid RNA nts may invade RNA structure.\n");
    printf("# 16. dG_Rloop_1nt: net free-energy change for one 5' bubble RNA nt invading the exit-channel RNA structure.\n");
    printf("# 17. dG_Rloop_2nt: net free-energy change for two 5' bubble RNA nts invading the exit-channel RNA structure.\n");
    printf("# 18. time: simulation time for the sampled state.\n");
  } else {
    printf("# 1. dG_RNA_RNA: RNA folding free energy of the current RNA secondary structure.\n");
    printf("# 2. time: simulation time for the sampled state.\n");
  }

  if (GTV.phi) {
    printf("# phi extras: zeitInc is the internal time increment, L is the accumulated Laplacian term, and D is the neighbor-count accumulator.\n");
  }
  if (GTV.verbose) {
    printf("# terminal verbose extras: top is the number of enumerated neighbors and lmin is the local-minimum status flag.\n");
  }
  printf("# terminal marker: O means the trajectory ended by time limit; Xn means stop structure n was reached.\n");
  fflush(stdout);
}

static void interval_metrics_reset(KinfoldIntervalMetrics *interval) {
  if (!interval)
    return;

  memset(interval, 0, sizeof(*interval));
  interval->initialized = 0;
}

static void interval_metrics_update(KinfoldIntervalMetrics *interval, const KinfoldInformationalMetrics *metrics) {
  double rna_energy, hybrid_energy, dna_duplex_energy, total_energy;

  if (!interval || !metrics)
    return;

  rna_energy = displayed_structure_energy();
  hybrid_energy = displayed_hybrid_energy();
  dna_duplex_energy = displayed_dna_duplex_energy();
  total_energy = displayed_energy();

#define UPDATE_RANGE(field_min, field_max, value) \
  do { \
    if (!interval->initialized) { \
      interval->field_min = (value); \
      interval->field_max = (value); \
    } else { \
      if ((value) < interval->field_min) interval->field_min = (value); \
      if ((value) > interval->field_max) interval->field_max = (value); \
    } \
  } while (0)

  UPDATE_RANGE(min_rna_energy, max_rna_energy, rna_energy);
  UPDATE_RANGE(min_hybrid_energy, max_hybrid_energy, hybrid_energy);
  UPDATE_RANGE(min_dna_duplex_energy, max_dna_duplex_energy, dna_duplex_energy);
  UPDATE_RANGE(min_total_energy, max_total_energy, total_energy);
  UPDATE_RANGE(min_terminal_hairpin_energy, max_terminal_hairpin_energy, metrics->terminal_hairpin_energy);
  UPDATE_RANGE(min_terminal_hairpin_dg_dt, max_terminal_hairpin_dg_dt, metrics->terminal_hairpin_dg_dt);
  UPDATE_RANGE(min_terminal_hairpin_plus1_energy, max_terminal_hairpin_plus1_energy, metrics->terminal_hairpin_plus1_energy);
  UPDATE_RANGE(min_rna_dna_minus1_ddg, max_rna_dna_minus1_ddg, metrics->rna_dna_minus1_ddg);
  UPDATE_RANGE(min_plus1_ddg, max_plus1_ddg, metrics->plus1_ddg);
  UPDATE_RANGE(min_exit_channel_bp_count, max_exit_channel_bp_count, (double)metrics->exit_channel_bp_count);
  UPDATE_RANGE(min_terminal_energy, max_terminal_energy, metrics->terminal_energy);
  UPDATE_RANGE(min_folding_inside_channel_energy, max_folding_inside_channel_energy, metrics->folding_inside_channel_energy);
  UPDATE_RANGE(min_elongation_rate_current, max_elongation_rate_current, metrics->elongation_rate_current);
  UPDATE_RANGE(min_rna_dna_1nt_shift_energy, max_rna_dna_1nt_shift_energy, metrics->rna_dna_1nt_shift_energy);
  UPDATE_RANGE(min_invasion_1nt_energy, max_invasion_1nt_energy, metrics->invasion_1nt_energy);
  UPDATE_RANGE(min_invasion_2nt_energy, max_invasion_2nt_energy, metrics->invasion_2nt_energy);
  UPDATE_RANGE(min_rloop_1nt_energy, max_rloop_1nt_energy, metrics->rloop_1nt_energy);
  UPDATE_RANGE(min_rloop_2nt_energy, max_rloop_2nt_energy, metrics->rloop_2nt_energy);

#undef UPDATE_RANGE

  interval->initialized = 1;
}

static void interval_seed_from_emitted_state(const KinfoldInformationalMetrics *metrics) {
  interval_metrics_reset(&current_interval_metrics);
  interval_metrics_update(&current_interval_metrics, metrics);
}

static void begin_trajectory_if_needed(void) {
  KinfoldInformationalMetrics metrics;

  if (csv_tracker_active && (csv_tracker_trajectory == GSV.trajectory_index))
    return;

  csv_tracker_active = 1;
  csv_tracker_trajectory = GSV.trajectory_index;
  csv_last_emitted_time = -1.0;
  interval_metrics_reset(&current_interval_metrics);
  interval_metrics_reset(&last_exported_interval_metrics);
  kinfold_metric_smoother_reset(&csv_metric_state);
  memset(&metrics, 0, sizeof(metrics));
  kinfold_compute_metrics_for_time(0.0, &csv_metric_state, &metrics);
  interval_metrics_update(&current_interval_metrics, &metrics);
}

static void track_hidden_state(double time_value) {
  KinfoldInformationalMetrics metrics;

  begin_trajectory_if_needed();
  memset(&metrics, 0, sizeof(metrics));
  kinfold_compute_metrics_for_time(time_value, &csv_metric_state, &metrics);
  interval_metrics_update(&current_interval_metrics, &metrics);
}

static void emit_csv_state_line(double out_time, int found_stop, int terminal, const KinfoldInformationalMetrics *metrics) {
  if (!csvFP || !metrics)
    return;

  fprintf(csvFP,
          "%d,%d,%.10g,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,",
          GSV.trajectory_index,
          GSV.steps,
          out_time,
          terminal ? 1 : 0,
          found_stop,
          GSV.tx_len,
          GSV.bubble_left,
          GSV.bubble_width,
          GSV.hybrid_left,
          GSV.hybrid_len,
          kinfold_channel_left(),
          kinfold_channel_length(),
          metrics->exit_channel_bp_count);
  csv_write_escaped(csvFP, transcription_active() ? kinfold_cli_structure_state() : costring(GAV.currform));
  fputc(',', csvFP);
  csv_write_escaped(csvFP, transcription_active() ? kinfold_full_sequence_state() : GAV.farbe);
  fprintf(csvFP,
          ",%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g",
          displayed_structure_energy(),
          displayed_hybrid_energy(),
          displayed_dna_duplex_energy(),
          displayed_energy(),
          metrics->terminal_hairpin_energy,
          metrics->terminal_hairpin_dg_dt,
          metrics->terminal_hairpin_plus1_energy,
          metrics->rna_dna_minus1_ddg,
          metrics->plus1_ddg,
          metrics->terminal_energy,
          metrics->folding_inside_channel_energy,
          metrics->elongation_rate_current,
          metrics->rna_dna_1nt_shift_energy,
          metrics->invasion_1nt_energy,
          metrics->invasion_2nt_energy,
          metrics->rloop_1nt_energy,
          metrics->rloop_2nt_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_rna_energy, current_interval_metrics.max_rna_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_hybrid_energy, current_interval_metrics.max_hybrid_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_dna_duplex_energy, current_interval_metrics.max_dna_duplex_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_total_energy, current_interval_metrics.max_total_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_terminal_hairpin_energy, current_interval_metrics.max_terminal_hairpin_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_terminal_hairpin_dg_dt, current_interval_metrics.max_terminal_hairpin_dg_dt);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_terminal_hairpin_plus1_energy, current_interval_metrics.max_terminal_hairpin_plus1_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_rna_dna_minus1_ddg, current_interval_metrics.max_rna_dna_minus1_ddg);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_plus1_ddg, current_interval_metrics.max_plus1_ddg);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_exit_channel_bp_count, current_interval_metrics.max_exit_channel_bp_count);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_terminal_energy, current_interval_metrics.max_terminal_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_folding_inside_channel_energy, current_interval_metrics.max_folding_inside_channel_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_elongation_rate_current, current_interval_metrics.max_elongation_rate_current);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_rna_dna_1nt_shift_energy, current_interval_metrics.max_rna_dna_1nt_shift_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_invasion_1nt_energy, current_interval_metrics.max_invasion_1nt_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_invasion_2nt_energy, current_interval_metrics.max_invasion_2nt_energy);
  fprintf(csvFP, ",%.10g,%.10g", current_interval_metrics.min_rloop_1nt_energy, current_interval_metrics.max_rloop_1nt_energy);
  fprintf(csvFP, ",%.10g,%.10g\n", current_interval_metrics.min_rloop_2nt_energy, current_interval_metrics.max_rloop_2nt_energy);
  fflush(csvFP);
  last_exported_interval_metrics = current_interval_metrics;
  interval_seed_from_emitted_state(metrics);
  csv_last_emitted_time = out_time;
}

static void emit_state_line(double out_time, int found_stop, int terminal) {
  int flag = 0;
  KinfoldInformationalMetrics metrics;
  const char *structure;

  if (!should_emit_to_stdout())
    return;

  if (!terminal && GTV.lmin && !(lmin == 1 && strcmp(GAV.prevform, GAV.currform) != 0))
    return;

  kinfold_get_sample_metrics(&metrics);
  emit_csv_state_line(out_time, found_stop, terminal, &metrics);
  structure = transcription_active() ? kinfold_cli_structure_state() : costring(GAV.currform);

  print_console_legend();

  if (GTV.structure_time_only && !structure_time_header_printed) {
    printf("%s\n", GAV.farbe_full ? GAV.farbe_full : GAV.farbe);
    structure_time_header_printed = 1;
    if (!transcription_active())
      structure = costring(GAV.currform);
  }

  if (GTV.structure_time_only) {
    char format[64];
    sprintf(format, "%%-%lus %%10.3f",
            (unsigned long)strlen(GAV.farbe_full ? GAV.farbe_full : GAV.farbe) + 1UL);
    printf(format, structure, out_time);
  } else if (kinfold_bubble_enabled()) {

    char format[256];
    sprintf(format, "%%-%lus %%6.2f %%6.2f %%6.2f %%6.2f %%8.4f %%6.2f %%6.2f %%6.2f %%4d %%6.2f %%6.2f %%9.3f %%6.2f %%6.2f %%6.2f %%6.2f %%6.2f %%10.3f",
            (unsigned long)strlen(GAV.farbe_full) + 1UL);
    printf(format,
           structure,
           displayed_structure_energy(),
           displayed_hybrid_energy(),
           displayed_dna_duplex_energy(),
           metrics.terminal_hairpin_energy,
           metrics.terminal_hairpin_dg_dt,
           metrics.terminal_hairpin_plus1_energy,
           metrics.rna_dna_minus1_ddg,
           metrics.plus1_ddg,
           metrics.exit_channel_bp_count,
           metrics.terminal_energy,
           metrics.folding_inside_channel_energy,
           metrics.elongation_rate_current,
           metrics.rna_dna_1nt_shift_energy,
           metrics.invasion_1nt_energy,
           metrics.invasion_2nt_energy,
           metrics.rloop_1nt_energy,
           metrics.rloop_2nt_energy,
           out_time);
  } else {
    char format[64];
    sprintf(format, "%%-%lus %%6.2f %%10.3f",
            (unsigned long)strlen(GAV.farbe_full) + 1UL);
    printf(format,
           structure,
           displayed_structure_energy(),
           out_time);
  }

  if (GTV.phi)
    printf(" %8.3f %8.3f %3g", zeitInc, L, D);

  if (GTV.verbose && terminal)
    printf(" %4d _ %d", top, lmin);

  if (terminal) {
    if (found_stop)
      printf(" X%d\n", found_stop);
    else
      printf(" O\n");
  } else {
    printf("\n");
  }

  fflush(stdout);
}

static void emit_sampled_states(double limit_time, int include_limit_time) {
  double limit = include_limit_time ? (limit_time + 1e-9) : (limit_time - 1e-9);

  while (GSV.next_sample_time <= limit) {
    kinfold_prepare_sample_metrics(GSV.next_sample_time);
    kinfold_capture_sample(GSV.next_sample_time);
    emit_state_line(GSV.next_sample_time, 0, 0);
    GSV.last_sample_time = GSV.next_sample_time;
    GSV.next_sample_time += GSV.time_step;
  }
}

static double compute_transition_rate(unsigned char kind, double curr_total, double total_E) {
  double dE;
  double elongation_propensity;

  dE = total_E - curr_total;
  if (kind == MOVE_ELONGATE) {
    elongation_propensity = kinfold_current_elongation_propensity();
    if (elongation_propensity > 0.0)
      return elongation_propensity;
  }

  if (GTV.mc) {
    if (dE < 0.0)
      return 1.0;
    return exp(-(dE / _RT * GSV.phi));
  }

  return exp(-0.5 * (dE / _RT * GSV.phi));
}

static void dump_transcription_state(FILE *fp) {
  if (!GTV.dump_transcription_neighbors)
    return;

  fprintf(fp,
          "#TXSTATE step=%d tx_len=%d full_len=%d channel_left=%d channel_length=%d elongation_rate=%.6f bubble_left=%d bubble_width=%d hybrid_left=%d hybrid_len=%d curr_rna_E=%.2f curr_total_E=%.2f seq=%s struct=%s\n",
          GSV.steps,
          GSV.tx_len,
          GSV.full_len,
          kinfold_channel_left(),
          kinfold_channel_length(),
          kinfold_current_elongation_rate(),
          GSV.bubble_left,
          GSV.bubble_width,
          GSV.hybrid_left,
          GSV.hybrid_len,
          GSV.currE,
          displayed_energy(),
          kinfold_full_sequence_state(),
          kinfold_full_structure_state());
}

static void dump_transcription_neighbor(FILE *fp,
                                       unsigned char kind,
                                       int i,
                                       int j,
                                       int next_tx_len,
                                       int next_bubble_left,
                                       int next_bubble_width,
                                       int next_hybrid_left,
                                       int next_hybrid_len,
                                       int rna_iE,
                                       double total_E,
                                       double rate) {
  const char *kind_name = "?";

  if (!GTV.dump_transcription_neighbors)
    return;

  switch (kind) {
    case MOVE_ELONGATE:      kind_name = "t"; break;
    case MOVE_HYBRID_FORM:   kind_name = "h"; break;
    case MOVE_HYBRID_BREAK:  kind_name = "H"; break;
    case MOVE_BUBBLE_GROW:   kind_name = "b"; break;
    case MOVE_BUBBLE_INVASION: kind_name = "B"; break;
    default: break;
  }

  fprintf(fp,
          "#TXNEIGH kind=%s i=%d j=%d next_tx_len=%d next_bubble_left=%d next_bubble_width=%d next_hybrid_left=%d next_hybrid_len=%d next_rna_E=%.2f next_total_E=%.2f rate=%.10g\n",
          kind_name,
          i,
          j,
          next_tx_len,
          next_bubble_left,
          next_bubble_width,
          next_hybrid_left,
          next_hybrid_len,
          (double)rna_iE / 100.0,
          total_E,
          rate);
}

static int eval_extended_structure_dcal(int new_len, const char *structure) {
  int old_len;
  char old_term;
  int eval_dcal;

  old_len  = GSV.len;
  old_term = GAV.farbe[new_len];

  memcpy(GAV.farbe, GAV.farbe_full, (size_t)new_len);
  GAV.farbe[new_len] = '\0';
  GSV.len = new_len;
#if HAVE_LIBRNA_API3
  GAV.vc->length = new_len;
#endif
  eval_dcal = kinfold_eval_structure_dcal(structure);

  GSV.len = old_len;
  memcpy(GAV.farbe, GAV.farbe_full, (size_t)old_len);
  GAV.farbe[old_len] = '\0';
  GAV.farbe[new_len] = old_term;
#if HAVE_LIBRNA_API3
  GAV.vc->length = old_len;
#endif

  return eval_dcal;
}

static void set_active_prefix_length(int new_len) {
  GSV.len = GSV.tx_len = new_len;
  memcpy(GAV.farbe, GAV.farbe_full, (size_t)new_len);
  GAV.farbe[new_len] = '\0';
#if HAVE_LIBRNA_API3
  GAV.vc->length = new_len;
#endif
}

static int current_rna_energy_dcal(void) {
  return (int)(GSV.currE * 100.0 + ((GSV.currE < 0.0) ? -0.4 : 0.4));
}

static void append_transcription_neighbors(void) {
  double total_E;
  double curr_total_before;
  int    rna_iE;
  int    i, j;
  int    saved_tx_len, saved_bubble_left, saved_bubble_width, saved_hybrid_left, saved_hybrid_len;

  if (!transcription_active())
    return;

  saved_tx_len       = GSV.tx_len;
  saved_bubble_left  = GSV.bubble_left;
  saved_bubble_width = GSV.bubble_width;
  saved_hybrid_left  = GSV.hybrid_left;
  saved_hybrid_len   = GSV.hybrid_len;
  curr_total_before  = displayed_energy();

  dump_transcription_state(stderr);

  if ((GSV.tx_len < GSV.full_len) && (kinfold_current_elongation_propensity() > 0.0)) {
    int new_len;
    int next_bubble_left = 0, next_bubble_width = 0, next_hybrid_left = 0, next_hybrid_len = 0;
    char *temp_structure;

    new_len = GSV.tx_len + 1;
    temp_structure = (char *)calloc((size_t)new_len + 1, sizeof(char));
    memcpy(temp_structure, GAV.currform, (size_t)GSV.tx_len);
    temp_structure[GSV.tx_len] = '.';
    temp_structure[new_len] = '\0';

    rna_iE = eval_extended_structure_dcal(new_len, temp_structure);
    total_E = (double)rna_iE / 100.0;

    if (kinfold_bubble_enabled()) {
      int old_bl, old_bw, old_hl, old_hlen;

      old_bl   = GSV.bubble_left;
      old_bw   = GSV.bubble_width;
      old_hl   = GSV.hybrid_left;
      old_hlen = GSV.hybrid_len;

      GSV.tx_len = new_len;
      if (old_bw <= 0) {
        GSV.bubble_width = 1;
        GSV.bubble_left  = new_len - 1;
        GSV.hybrid_len   = 1;
        GSV.hybrid_left  = new_len - 1;
      } else if (old_bw < GSV.max_bubble_width) {
        GSV.bubble_width = old_bw + 1;
        GSV.bubble_left  = new_len - GSV.bubble_width;
        GSV.hybrid_len   = old_hlen + 1;
        GSV.hybrid_left  = new_len - GSV.hybrid_len;
      } else {
        GSV.bubble_width = GSV.max_bubble_width;
        GSV.bubble_left  = new_len - GSV.bubble_width;
        if (old_hlen <= 0) {
          GSV.hybrid_len  = 1;
          GSV.hybrid_left = new_len - 1;
        } else if (old_hl == old_bl) {
          GSV.hybrid_len  = old_hlen;
          GSV.hybrid_left = old_hl + 1;
        } else {
          GSV.hybrid_len  = old_hlen + 1;
          if (GSV.hybrid_len > GSV.bubble_width)
            GSV.hybrid_len = GSV.bubble_width;
          GSV.hybrid_left = new_len - GSV.hybrid_len;
        }
      }
      next_bubble_left = GSV.bubble_left;
      next_bubble_width = GSV.bubble_width;
      next_hybrid_left = GSV.hybrid_left;
      next_hybrid_len = GSV.hybrid_len;
      total_E += kinfold_current_hybrid_energy();
      GSV.tx_len        = saved_tx_len;
      GSV.bubble_left   = old_bl;
      GSV.bubble_width  = old_bw;
      GSV.hybrid_left   = old_hl;
      GSV.hybrid_len    = old_hlen;
    }

    dump_transcription_neighbor(stderr,
                                MOVE_ELONGATE,
                                0,
                                0,
                                new_len,
                                next_bubble_left,
                                next_bubble_width,
                                next_hybrid_left,
                                next_hybrid_len,
                                rna_iE,
                                total_E,
                                compute_transition_rate(MOVE_ELONGATE, curr_total_before, total_E));
    update_nbList_kind(MOVE_ELONGATE, 0, 0, rna_iE, total_E);
    free(temp_structure);
  }

  if (!kinfold_bubble_enabled())
    return;

  if (GSV.bubble_left < GSV.hybrid_left) {
    GSV.hybrid_left = saved_hybrid_left - 1;
    GSV.hybrid_len  = saved_hybrid_len + 1;
    total_E = GSV.currE + kinfold_current_hybrid_energy();
    dump_transcription_neighbor(stderr,
                                MOVE_HYBRID_FORM,
                                0,
                                0,
                                GSV.tx_len,
                                GSV.bubble_left,
                                GSV.bubble_width,
                                GSV.hybrid_left,
                                GSV.hybrid_len,
                                current_rna_energy_dcal(),
                                total_E,
                                compute_transition_rate(MOVE_HYBRID_FORM, curr_total_before, total_E));
    update_nbList_kind(MOVE_HYBRID_FORM, 0, 0, current_rna_energy_dcal(), total_E);
    GSV.hybrid_left = saved_hybrid_left;
    GSV.hybrid_len  = saved_hybrid_len;
  }

  if (GSV.hybrid_len > 0) {
    GSV.hybrid_left = saved_hybrid_left + 1;
    GSV.hybrid_len  = saved_hybrid_len - 1;
    total_E = GSV.currE + kinfold_current_hybrid_energy();
    dump_transcription_neighbor(stderr,
                                MOVE_HYBRID_BREAK,
                                0,
                                0,
                                GSV.tx_len,
                                GSV.bubble_left,
                                GSV.bubble_width,
                                GSV.hybrid_left,
                                GSV.hybrid_len,
                                current_rna_energy_dcal(),
                                total_E,
                                compute_transition_rate(MOVE_HYBRID_BREAK, curr_total_before, total_E));
    update_nbList_kind(MOVE_HYBRID_BREAK, 0, 0, current_rna_energy_dcal(), total_E);
    GSV.hybrid_left = saved_hybrid_left;
    GSV.hybrid_len  = saved_hybrid_len;
  }

  if ((GSV.bubble_width < GSV.max_bubble_width) &&
      (GSV.bubble_left > 0) &&
      kinfold_is_unpaired_position(GSV.bubble_left - 1)) {
    GSV.bubble_left  = saved_bubble_left - 1;
    GSV.bubble_width = saved_bubble_width + 1;
    GSV.hybrid_left  = saved_hybrid_left - 1;
    GSV.hybrid_len   = saved_hybrid_len + 1;
    total_E = GSV.currE + kinfold_current_hybrid_energy();
    dump_transcription_neighbor(stderr,
                                MOVE_BUBBLE_GROW,
                                0,
                                0,
                                GSV.tx_len,
                                GSV.bubble_left,
                                GSV.bubble_width,
                                GSV.hybrid_left,
                                GSV.hybrid_len,
                                current_rna_energy_dcal(),
                                total_E,
                                compute_transition_rate(MOVE_BUBBLE_GROW, curr_total_before, total_E));
    update_nbList_kind(MOVE_BUBBLE_GROW, 0, 0, current_rna_energy_dcal(), total_E);
    GSV.bubble_left  = saved_bubble_left;
    GSV.bubble_width = saved_bubble_width;
    GSV.hybrid_left  = saved_hybrid_left;
    GSV.hybrid_len   = saved_hybrid_len;
  }

  if ((GSV.hybrid_len > 0) &&
      kinfold_is_unpaired_position(GSV.hybrid_left)) {
    char *temp_structure;

    temp_structure = strdup(GAV.currform);
    assert(temp_structure != NULL);
    for (i = 0; i < GSV.hybrid_left; i++) {
      if (!position_in_exit_channel_current(i))
        continue;
      if (!kinfold_is_unpaired_position(i))
        continue;
      if (!kinfold_can_pair_positions_raw(i, saved_hybrid_left))
        continue;

      if (i < saved_hybrid_left) {
        temp_structure[i] = '(';
        temp_structure[saved_hybrid_left] = ')';
      } else {
        temp_structure[saved_hybrid_left] = '(';
        temp_structure[i] = ')';
      }

      rna_iE = kinfold_eval_structure_dcal(temp_structure);
      GSV.hybrid_left = saved_hybrid_left + 1;
      GSV.hybrid_len = saved_hybrid_len - 1;
      if (GSV.hybrid_len <= 0) {
        GSV.hybrid_len = 0;
        GSV.hybrid_left = GSV.len;
      }
      total_E = (double)rna_iE / 100.0 + kinfold_current_hybrid_energy();
      dump_transcription_neighbor(stderr,
                                  MOVE_BUBBLE_INVASION,
                                  i + 1,
                                  saved_hybrid_left + 1,
                                  GSV.tx_len,
                                  GSV.bubble_left,
                                  GSV.bubble_width,
                                  GSV.hybrid_left,
                                  GSV.hybrid_len,
                                  rna_iE,
                                  total_E,
                                  compute_transition_rate(MOVE_BUBBLE_INVASION, curr_total_before, total_E));
      update_nbList_kind(MOVE_BUBBLE_INVASION,
                         i + 1,
                         saved_hybrid_left + 1,
                         rna_iE,
                         total_E);

      GSV.hybrid_left = saved_hybrid_left;
      GSV.hybrid_len = saved_hybrid_len;
      temp_structure[i] = '.';
      temp_structure[saved_hybrid_left] = '.';
    }
    free(temp_structure);
  }

  if ((GSV.bubble_left < GSV.hybrid_left) &&
      kinfold_is_unpaired_position(GSV.bubble_left)) {
    char *temp_structure;

    temp_structure = strdup(GAV.currform);
    for (j = 0; j < GSV.len; j++) {
      i = GSV.bubble_left;
      if (j == i)
        continue;
      if ((j >= GSV.bubble_left) && (j < GSV.len))
        continue;
      if (!kinfold_is_unpaired_position(j))
        continue;
      GSV.bubble_left  = saved_bubble_left + 1;
      GSV.bubble_width = saved_bubble_width - 1;
      if (!kinfold_can_pair_positions(i, j)) {
        GSV.bubble_left  = saved_bubble_left;
        GSV.bubble_width = saved_bubble_width;
        continue;
      }

      if (i < j) {
        temp_structure[i] = '(';
        temp_structure[j] = ')';
      } else {
        temp_structure[j] = '(';
        temp_structure[i] = ')';
      }

      rna_iE = kinfold_eval_structure_dcal(temp_structure);
      GSV.bubble_left  = saved_bubble_left + 1;
      GSV.bubble_width = saved_bubble_width - 1;
      if (GSV.bubble_width <= 0) {
        GSV.bubble_width = 0;
        GSV.bubble_left  = GSV.len;
        GSV.hybrid_left  = GSV.len;
        GSV.hybrid_len   = 0;
      }
      total_E = (double)rna_iE / 100.0 + kinfold_current_hybrid_energy();
      dump_transcription_neighbor(stderr,
                                  MOVE_BUBBLE_INVASION,
                                  i + 1,
                                  j + 1,
                                  GSV.tx_len,
                                  GSV.bubble_left,
                                  GSV.bubble_width,
                                  GSV.hybrid_left,
                                  GSV.hybrid_len,
                                  rna_iE,
                                  total_E,
                                  compute_transition_rate(MOVE_BUBBLE_INVASION, curr_total_before, total_E));
      update_nbList_kind(MOVE_BUBBLE_INVASION, i + 1, j + 1, rna_iE, total_E);

      GSV.bubble_left  = saved_bubble_left;
      GSV.bubble_width = saved_bubble_width;
      GSV.hybrid_left  = saved_hybrid_left;
      GSV.hybrid_len   = saved_hybrid_len;
      temp_structure[i] = '.';
      temp_structure[j] = '.';
    }
    free(temp_structure);
  }
}

static void apply_transcription_move(int next) {
  unsigned char kind;

  kind = move_kinds[next];

  switch (kind) {
    case MOVE_ELONGATE:
      set_active_prefix_length(GSV.tx_len + 1);
      GAV.currform[GSV.len - 1] = '.';
      GAV.currform[GSV.len] = '\0';
      if (kinfold_bubble_enabled()) {
        if (GSV.bubble_width <= 0) {
          GSV.bubble_width = 1;
          GSV.bubble_left  = GSV.len - 1;
          GSV.hybrid_len   = 1;
          GSV.hybrid_left  = GSV.len - 1;
        } else if (GSV.bubble_width < GSV.max_bubble_width) {
          GSV.bubble_width++;
          GSV.bubble_left = GSV.len - GSV.bubble_width;
          GSV.hybrid_len++;
          GSV.hybrid_left = GSV.len - GSV.hybrid_len;
        } else {
          GSV.bubble_left = GSV.len - GSV.bubble_width;
          if (GSV.hybrid_len <= 0) {
            GSV.hybrid_len = 1;
          } else if (GSV.hybrid_left == GSV.bubble_left - 1) {
            GSV.hybrid_left++;
          } else {
            GSV.hybrid_len++;
            if (GSV.hybrid_len > GSV.bubble_width)
              GSV.hybrid_len = GSV.bubble_width;
          }
          GSV.hybrid_left = GSV.len - GSV.hybrid_len;
        }
      }
      GSV.currE = energies[next];
      kinfold_rebuild_current_state();
      break;

    case MOVE_HYBRID_FORM:
      GSV.hybrid_left--;
      GSV.hybrid_len++;
      break;

    case MOVE_HYBRID_BREAK:
      GSV.hybrid_left++;
      GSV.hybrid_len--;
      if (GSV.hybrid_len <= 0) {
        GSV.hybrid_len = 0;
        GSV.hybrid_left = GSV.len;
      }
      break;

    case MOVE_BUBBLE_GROW:
      GSV.bubble_left--;
      GSV.bubble_width++;
      GSV.hybrid_left--;
      GSV.hybrid_len++;
      break;

    case MOVE_BUBBLE_INVASION: {
      int i, j;
      int targets_hybrid_left;

      i = neighbor_list[2*next] - 1;
      j = neighbor_list[2*next+1] - 1;
      targets_hybrid_left = ((i == GSV.hybrid_left) || (j == GSV.hybrid_left)) ? 1 : 0;
      if (i < j) {
        GAV.currform[i] = '(';
        GAV.currform[j] = ')';
      } else {
        GAV.currform[j] = '(';
        GAV.currform[i] = ')';
      }
      if (targets_hybrid_left) {
        GSV.hybrid_left++;
        GSV.hybrid_len--;
        if (GSV.hybrid_len <= 0) {
          GSV.hybrid_len = 0;
          GSV.hybrid_left = GSV.len;
        }
      } else {
        GSV.bubble_left++;
        GSV.bubble_width--;
        if (GSV.bubble_width <= 0) {
          GSV.bubble_width = 0;
          GSV.bubble_left = GSV.len;
          GSV.hybrid_left = GSV.len;
          GSV.hybrid_len = 0;
        }
      }
      GSV.currE = energies[next];
      kinfold_rebuild_current_state();
      break;
    }
  }
}

/*======================*/
void clean_up_nbList(void){

  free(neighbor_list);
  free(bmf);
  free(energies);
  free(move_kinds);
  neighbor_list = NULL;
  bmf = NULL;
  energies = NULL;
  move_kinds = NULL;
  if (logFP != NULL) {
    fprintf(logFP,"\n");
    fclose(logFP);
    logFP = NULL;
  }
  if (csvFP != NULL) {
    fclose(csvFP);
    csvFP = NULL;
  }
}

/*======================*/
static void grow_chain(void){
  int newl;
  /* note Zeit=0 corresponds to chain length GSV.glen */
  if (Zeit<(GSV.len+1-GSV.glen) * GSV.grow) return;
  newl = GSV.len+1;
  Zeit = (newl-GSV.glen) * GSV.grow;
  top=0; /* prevent structure move in sel_nb */

  if (GSV.len<newl) {
    strncpy(GAV.farbe, GAV.farbe_full, newl);
    GAV.farbe[newl] = '\0';
    strcpy(GAV.startform, GAV.currform);
    strcat(GAV.startform, ".");

    GSV.len = newl;
#if HAVE_LIBRNA_API3
    /* fake actual length of sequence in GAV.vc */
    GAV.vc->length = newl;
#endif
  }
}

char
kinfold_move_char(int next)
{
  int ii, jj;

  if (next < 0)
    return 'g';

  if (move_kinds[next] == MOVE_ELONGATE)
    return 't';
  if (move_kinds[next] == MOVE_HYBRID_FORM)
    return 'h';
  if (move_kinds[next] == MOVE_HYBRID_BREAK)
    return 'H';
  if (move_kinds[next] == MOVE_BUBBLE_GROW)
    return 'b';
  if (move_kinds[next] == MOVE_BUBBLE_INVASION)
    return 'B';

  ii = neighbor_list[2 * next];
  jj = neighbor_list[2 * next + 1];

  if (abs(ii) < GSV.len) {
    if ((ii > 0) && (jj > 0))
      return 'i';
    if ((ii < 0) && (jj < 0))
      return 'd';
    if ((ii > 0) && (jj < 0))
      return 's';
    return 'S';
  }

  if ((ii > 0) && (jj > 0))
    return 'I';

  return 'D';
}

void
kinfold_get_interval_metrics_for_time(double time_value,
                                      KinfoldIntervalMetrics *metrics)
{
  if (!metrics)
    return;

  if (csv_tracker_active && (fabs(csv_last_emitted_time - time_value) <= 1e-9))
    *metrics = last_exported_interval_metrics;
  else
    *metrics = current_interval_metrics;
}

static const char *costring(const char *str) {
  static char* buffer=NULL;
  static int size=0;
  int n;
  if (str==NULL) {
    if (buffer) {
      /* make it possible to free buffer */
      free(buffer);
      size = 0; buffer = NULL;
    }
    return NULL;
  }
  n=strlen(str);
  if (n>=size) {
    size = n+2;
    buffer = realloc(buffer, size);
  }
  if ((cut_point>0)&&(cut_point<=n)) {
    strncpy(buffer, str, cut_point-1);
    buffer[cut_point-1] = '&';
    strncpy(buffer+cut_point, str+cut_point-1, n-cut_point+1);
    buffer[n+1] = '\0';
  } else {
    strncpy(buffer, str, n+1);
  }
  return buffer;
}
