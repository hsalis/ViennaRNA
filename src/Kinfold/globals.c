/*
  Last changed Time-stamp: <2015-04-13 16:08:40 ivo>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: globals.c,v 1.8 2008/10/07 09:03:14 ivo Exp $
*/

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <float.h>
#include <math.h>
#include <ctype.h>

#if HAVE_LIBRNA_API3
#include <ViennaRNA/utils.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold_compound.h>
#else
#include <utils.h>
#include <fold_vars.h>
#include <fold.h>
#endif

#include "globals.h"
#include "cmdline.h"
#include "baum.h"


GlobVars GSV;
GlobArrays GAV;
GlobToggles GTV;

/* forward declarations privat functions */
static void ini_globs(void);
static void ini_gtoggles (void);
static void ini_gvars(void);
static void ini_garrays (void);
static void usage(int status);
static void display_settings(void);
static void display_fileformat(void);
static char *verbose(int optval, const char *which);
static int process_options (int argc, char *argv[]);
static void process_options_gg (int argc, char *argv[]);

static const char *costring(const char *str);
static void derive_template_dna(void);
static void update_transcription_mode_defaults(void);
static void update_exit_channel_state(void);
static void normalize_transcription_start_structure(void);
static int  is_transcribed_position(int idx0);
static char dna_complement(char rna_nt);
static char coding_dna_nt(int idx0);
static int  hybrid_segment_energy_dcal(int left, int len);
static int  hybrid_shifted_upstream_energy_dcal(int left, int len, int *supported);
static int  dna_duplex_segment_energy_dcal(int left, int len);
static int  eval_substructure_dcal(const char *sequence, const char *structure);
static int  build_pair_table(const char *structure, int len, int *pair_table, int *parent_open);
static int  compute_terminal_hairpin_metrics_raw(KinfoldInformationalMetrics *metrics);
static int  build_substructure_with_outside_pairs_dotted(const char *sequence,
                                                         const char *structure,
                                                         int        left,
                                                         int        right,
                                                         char       **subseq,
                                                         char       **substruc);
static int  pair_crosses_existing(const int *pair_table, int len, int i, int j);
static int  evaluate_invasion_energy(int invade_count);
static double terminal_component_energy_for_structure(const char *structure);
static int rewound_dna_duplex_energy_dcal(int released_hybrid_nts);
static double evaluate_rloop_terminal_plus_energy(int pair_count, int *found);
static double pause_rate_scale_for_state(void);
static int position_in_exit_channel(int idx0);
static int exit_channel_left_state(void);
static int exit_channel_right_state(void);
static const char *full_sequence_state(void);
static const char *full_structure_state(void);
static const char *full_structure_state_from(const char *structure);
static const char *state_cache_key(void);

static char UNUSED rcsid[] ="$Id: globals.c,v 1.8 2008/10/07 09:03:14 ivo Exp $";
static KinfoldInformationalMetrics current_sample_metrics;
static double sample_metrics_time = 0.0;
static double sample_metrics_smoothed_hairpin = 0.0;
static int sample_metrics_ready = 0;
static int sample_metrics_initialized = 0;
#define MAXMSG 8
static char msg[MAXMSG][60] =
{{"off"},
 {"on"},
 {"current time"},
 {"take from input file"},
 {"open chain"},
 {"mfe structure of sequence"},
 {"ViennaRNA-Package-1.4 defaults"},
 {""}
};

static struct option const long_options[] =
{ {"dangle",  required_argument, 0, 0},
  {"Temp",    required_argument, 0, 0},
  {"Par",     required_argument, 0, 0},
  {"phi",     required_argument, 0, 0},
  {"pbounds", required_argument, 0, 0},
  {"logML",   no_argument,       &GTV.logML, 0},
  {"noShift", no_argument,       &GTV.noShift, 1},
  {"noLP",    no_argument,       &GTV.noLP, 1},
  {"seed",    required_argument, 0, 0},
  {"time",    required_argument, 0, 0},
  {"num",     required_argument, 0, 0},
  {"start",   no_argument,       &GTV.start, 1},
  {"stop",    no_argument,       &GTV.stop, 1},
  {"fpt",     no_argument,       &GTV.fpt, 0},
  {"rect",     no_argument,       &GTV.rect, 1},
  {"met",      no_argument,       &GTV.mc, 1},
  {"grow",    required_argument, 0,  0},
  {"glen",    required_argument, 0,  0},
  {"log",     required_argument, 0,  0},
  {"silent",  no_argument,       0,  0},
  {"lmin",    no_argument,       &GTV.lmin, 1},
  {"cut",     required_argument, 0, 0},
  {"help",    no_argument,       0, 'h'},
  {"verbose", no_argument,       0, 0},
  {NULL, 0, NULL, 0}
};

/**/
void decode_switches(int argc, char *argv[]) {
  ini_globs();
  strcpy(GAV.ProgramName, argv[0]);
  process_options_gg(argc, argv);
}

/**/
void clean_up_globals(void) {
  int i;
  free(GAV.ProgramName);
  free(GAV.BaseName);
  free(GAV.CSVOutputFile);
  free(GAV.farbe);
  free(GAV.farbe_full);
  free(GAV.template_dna_full);
  free(GAV.startform);
  free(GAV.currform);
  free(GAV.prevform);
  for (i = 0; i < GSV.maxS; i++) free(GAV.stopform[i]);
  free(GAV.stopform);
  free(GAV.sE);
#if HAVE_RNALIB_API3
  vrna_fold_compound_free(GAV.vc);
#else
# if HAVE_RNALIB2
  free(GAV.params);
# endif
#endif
}

/**/
static void usage(int status) {
  fprintf(stderr, "\n%s - Kinetic Folding Program for Nucleic Acids -\n",
	  GAV.ProgramName);
  fprintf(stderr, "Usage: %s [OPTION] < FILE\n", GAV.ProgramName);
  fprintf(stderr,
	  "Options:\n"
	  " EnergyModel\n"
	  "  --dangle <0|1|2>  set dangling end model to (non|normal|double)\n"
	  "  --Temp <float>    set simulation temperature to <float>\n"
	  "  --Par <string>    use energy-parameter-file <string>\n"
	  "  --logML           use linear multiloop-function not logarithmic\n"
	  " MoveSet\n"
	  "  --noShift         turn off shift-moves\n"
	  "  --noLP            forbit structures with isolated base-pairs\n"
	  " Simulation\n"
	  "  --seed <int=int=int>  set random seed to <int=int=int>\n"
	  "  --time <float>        set maxtime of simulation to <float>\n"
	  "  --time-step <float>   save/output trajectory state every <float> seconds\n"
	  "  --num <int>           set number of simulations to <int>\n"
	  "  --start               set start structure\n"
	  "  --stop                set stop structure(s)\n"
	  "  --met                 use Metropolis rule not Kawasaki rule\n"
	  "  --fpt                 stop stop structure(s) is reached\n"
	  "  --rect                stop start structure is reached\n"
	  "  --grow <float>        deprecated alias for transcription elongation\n"
	  "  --transcription_elongation_rate <float>\n"
	  "                        elongate transcription at <float> nt/s\n"
	  "  --kinfold_seconds_per_time_unit <float>\n"
	  "                        physical seconds per Kinfold time unit\n"
	  "  --max_bubble_width <int>\n"
	  "                        maximum transcription bubble width\n"
	  "  --RNAP_exit_channel_length <int>\n"
	  "                        RNAP RNA exit-channel length in nucleotides\n"
	  "  --phi <double>        set phi value to <double>\n"
	  "  --pbounds <d1=d2=d3>  set phi_min to d1\n"
	  "                            phi_inc to d2\n"
	  "                            phi_max to d2\n"
	  "                            (d? is a double value)\n"
	  " Output\n"
	  "  --log <string>  set basename of log-file to <string>\n"
	  "  --csv-output <string>\n"
	  "                  write sampled trajectory states and metrics to CSV file\n"
	  "  --structure-time-only\n"
	  "                  print sequence once, then only structure and time to stdout\n"
	  "  --err <string>  set basename of error-log-file to <string>\n"
	  "  --silent        no output to stdout\n"
	  "  --verbose       more information to stdout\n"
	  "  --lmin          output only local minima to stdout\n"
	  "  --cut <float>   output structures with E <= <float> to stdout (default: no cutoff)\n");
  display_settings();
  display_fileformat();
  exit (status);
}

/**/
static void display_fileformat(void) {
  fprintf(stderr,
	  "Input File Format:\n"
	  "1st line sequence\n"
	  "2nd line start structure (if option --start is used)\n"
	  "following lines stop structures\n\n");
}

/**/
void log_prog_params(FILE *FP) {
    const char *cut_value = (GSV.cut >= (float)(DBL_MAX / 2.0)) ? "inf" : NULL;
    fprintf( FP,
	     "#<\n#Date: %s"
	     "#EnergyModel: dangle=%d Temp=%.1f logML=%s Par=%s\n"
	     "#MoveSet: noShift=%s noLP=%s\n"
	     "#Simulation: num=%d time=%.2f time_step=%.2f seed=%s fpt=%s rect=%s mc=%s\n"
	     "#Simulation: transcription=%s elongation_rate_nt_s=%.6g seconds_per_time_unit=%.6g max_bubble_width=%d exit_channel_length=%d tx_len=%d bubble=[%d,%d) hybrid=[%d,%d) channel=[%d,%d)\n"
	     "#Simulation: phi=%g pbounds=%s\n"
	     "#Output: log=%s csv=%s structure_time_only=%s silent=%s lmin=%s cut=%s",
	     time_stamp(),
	     GTV.dangle,
	     GSV.Temp,
	     verbose(GTV.logML, "logML"),
	     GAV.ParamFile,
	     verbose(GTV.noShift, NULL),
	     verbose(GTV.noLP, NULL),
	     GSV.num,
	     GSV.time,
             GSV.time_step,
	     verbose(GTV.seed, "seed"),
	     verbose(GTV.fpt, NULL),
	     verbose(GTV.rect, NULL),
	     verbose(GTV.mc, "met"),
	     verbose(GTV.transcription, NULL),
	     GSV.transcription_elongation_rate_nt_per_s,
	     GSV.seconds_per_internal_time_unit,
	     GSV.max_bubble_width,
	     GSV.RNAP_exit_channel_length,
	     GSV.tx_len,
	     GSV.bubble_left,
	     GSV.bubble_left + GSV.bubble_width,
	     GSV.hybrid_left,
	     GSV.hybrid_left + GSV.hybrid_len,
	     GSV.channel_left,
	     GSV.channel_left + GSV.channel_length,
	     GSV.phi,
	     verbose(GTV.phi, "pbounds"),
	     GAV.BaseName,
	     GAV.CSVOutputFile,
	     verbose(GTV.structure_time_only, NULL),
	     verbose(GTV.silent, NULL),
	     verbose(GTV.lmin, NULL),
	     cut_value ? cut_value : "");
    if (!cut_value)
      fprintf(FP, "%.2f", GSV.cut);
    fprintf(FP, "\n");
    fflush(FP);
}

/**/
void log_start_stop(FILE *FP) {
  int i;
  fprintf(FP, "#%s\n", GTV.transcription ? kinfold_full_sequence_state() : costring(GAV.farbe));
  /* SB: I don't know why, but having two costring calls in one print function yields
   * the same return value twice, so let's just use two lines of code. */
  fprintf(FP, "#%s (%6.2f)\n",
          GTV.transcription ? full_structure_state_from(GAV.startform) : costring(GAV.startform),
          GSV.startE);
  for (i = 0; i < GSV.maxS; i++) {
    fprintf(FP, "#%s (%6.2f) X%02d\n", costring(GAV.stopform[i]), GAV.sE[i], i+1);
  }
  costring(NULL);
  fflush(FP);
}

void kinfold_finalize_transcription_input(void) {
  if (!GTV.transcription)
    return;

  normalize_transcription_start_structure();
  derive_template_dna();
}

/**/
static void display_settings(void) {
  fprintf(stderr,
	  "Default Settings:\n"
	  " EnergyModel\n"
	  "  --dangle  = %d\n"
	  "  --Temp    = %.2f\n"
	  "  --Par     = %s\n"
	  "  --logML   = %s\n"
	  " MoveSet\n"
	  "  --noShift = %s\n"
	  "  --noLP    = %s\n"
	  " Simulation\n"
	  "  --seed    = %s\n"
	  "  --time    = %.1f\n"
	  "  --time-step = %.2f\n"
	  "  --phi     = %g\n"
	  "  --pbounds = %s\n"
	  "  --num     = %d\n"
	  "  --start   = %s\n"
	  "  --stop    = %s\n"
	  "  --met     = %s\n"
	  "  --transcription_elongation_rate = %.6g\n"
	  "  --kinfold_seconds_per_time_unit = %.6g\n"
	  "  --max_bubble_width = %d\n"
	  "  --RNAP_exit_channel_length = %d\n"
		  "  --fpt     = %s\n"
		  "  --rect     = %s\n"
		  " Output\n"
		  "  --log     = %s\n"
		  "  --csv-output = %s\n"
		  "  --structure-time-only = %s\n"
		  "  --silent  = %s\n"
	  "  --verbose = %s\n"
	  "  --lmin    = %s\n"
	  "  --cut     = %s",
	  GTV.dangle,
	  GSV.Temp,
	  GAV.ParamFile,
	  verbose(GTV.logML, "logML"),
	  verbose(GTV.noShift, NULL),
	  verbose(GTV.noLP, NULL),
	  verbose(GTV.seed, "seed"),
	  GSV.time,
	  GSV.time_step,
	  GSV.phi,
	  verbose(GTV.phi, "pbounds"),
	  GSV.num,
	  verbose(GTV.start, "start"),
	  verbose(GTV.stop, "stop"),
	  verbose(GTV.mc, "met"),
	  GSV.transcription_elongation_rate_nt_per_s,
	  GSV.seconds_per_internal_time_unit,
	  GSV.max_bubble_width,
	  GSV.RNAP_exit_channel_length,
		  verbose(GTV.fpt, NULL),
		  verbose(GTV.rect, NULL),
		  GAV.BaseName,
		  GAV.CSVOutputFile,
		  verbose(GTV.structure_time_only, NULL),
		  verbose(GTV.silent, NULL),
	  verbose(GTV.verbose, NULL),
	  verbose(GTV.lmin, NULL),
	  (GSV.cut >= (float)(DBL_MAX / 2.0)) ? "inf" : "");
  if (GSV.cut < (float)(DBL_MAX / 2.0))
    fprintf(stderr, "%.2f", GSV.cut);
  fprintf(stderr, "\n");
}

/**/
static char *verbose (int optval, const char *which) {
  if (which == NULL) return msg[optval];
  else {
    if ( strcmp(which, "seed") == 0 ) {
      if (optval == 0 ) return "clock";
      else {
	sprintf(msg[MAXMSG - 1],
		"%hu %hu %hu", GAV.subi[0], GAV.subi[1],GAV.subi[2]);
	return msg[MAXMSG - 1];
      }
    }

    if (strcmp(which, "pbounds") == 0) {
      sprintf(msg[MAXMSG - 1],
	      "%g %g %g",
	      GAV.phi_bounds[0], GAV.phi_bounds[1], GAV.phi_bounds[2]);
      return msg[MAXMSG - 1];
    }

    if ( strcmp(which, "met") == 0 ) {
      if ( optval == 0 ) return "Kawasaki";
      else return "Metropolis";
    }

    if ( strcmp(which, "logML") == 0 ) {
      if ( optval == 0 ) return "linear";
      else return "logarithmic";
    }

    if ( strcmp(which, "start") == 0 ) {
      if ( optval == 0 ) return "OpenChain";
      else return "input file";
    }

    if ( strcmp(which, "stop") == 0 ) {
      if ( optval == 0 ) return "Mfe";
      else return "input file";
    }
  }
  return "hae?";
}

/**/
static void ini_globs (void) {
  ini_gtoggles();
  ini_gvars();
  ini_garrays();
}

static void process_options_gg (int argc, char *argv[]) {
  struct gengetopt_args_info args_info;
  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(1) ;
  if (args_info.help_given) cmdline_parser_print_help();
  if (args_info.full_help_given) cmdline_parser_print_full_help();

  GTV.dangle = args_info.dangle_arg;
  GSV.Temp = args_info.Temp_arg;
  if (args_info.Par_given) 
    strncpy(GAV.ParamFile,args_info.Par_arg,255);
  if (args_info.csv_output_given)
    strncpy(GAV.CSVOutputFile, args_info.csv_output_arg, 255);
  GTV.Par = args_info.Par_given;
  GTV.logML = args_info.logML_flag;

  GTV.noShift = args_info.noShift_flag;
  GTV.noLP = args_info.noLP_flag;

  GTV.start= args_info.start_flag;
  GTV.stop = args_info.stop_flag;
  GTV.mc   = args_info.met_flag;
  GTV.lmin = args_info.lmin_flag;

  GTV.verbose = args_info.verbose_flag;
  GTV.structure_time_only = args_info.structure_time_only_flag;
  GTV.silent = (GTV.verbose==0)?
    args_info.silent_flag : 0;
  if (args_info.seed_given)
    if (sscanf(args_info.seed_arg, "%hu=%hu=%hu",
	       &GAV.subi[0], &GAV.subi[1], &GAV.subi[2]) != 3)
      usage(EXIT_FAILURE);
    else GTV.seed = 1;

  /* laplace stuff */
  if (args_info.phi_given) {
    if (args_info.phi_arg>0) {
      GTV.phi = 1;
      GSV.phi = args_info.phi_arg;
    }
    else {
      fprintf(stderr, "Value of --phi must be > 0 >%lf<\n", args_info.phi_arg);
      exit(EXIT_FAILURE);
    }
  }
  if (args_info.pbounds_given) {
    if (sscanf(args_info.pbounds_arg, "%g=%g=%g",
	       &GAV.phi_bounds[0],
	       &GAV.phi_bounds[1],
	       &GAV.phi_bounds[2]) == 0)
      usage(EXIT_FAILURE);
    else
      GTV.phi = 1;
    /* check if values are proper */
    if (GAV.phi_bounds[0] > GAV.phi_bounds[2]
	|| GAV.phi_bounds[1] > GAV.phi_bounds[2]) {
      fprintf(stderr,
	      "Unmet requirements for pbounds:\n"
	      "phi_min < phi_max && phi_inc < phi_max\n"
	      "phi_min: %g phi_inc: %g phi_max: %g\n",
	      GAV.phi_bounds[0], GAV.phi_bounds[1], GAV.phi_bounds[2]);
      exit(EXIT_FAILURE);
    }
  }
  GSV.time = args_info.time_arg;
  GSV.time_step = args_info.time_step_arg;
  GSV.num = args_info.num_arg;
  strncpy(GAV.BaseName, args_info.log_arg, 255);
  GSV.cut = args_info.cut_given ? args_info.cut_arg : (float)DBL_MAX;
  GSV.grow = args_info.grow_arg;
  GSV.glen = args_info.glen_arg;
  GSV.transcription_elongation_rate_nt_per_s = args_info.transcription_elongation_rate_arg;
  GSV.seconds_per_internal_time_unit = args_info.kinfold_seconds_per_time_unit_arg;
  GSV.max_bubble_width = args_info.max_bubble_width_arg;
  GSV.RNAP_exit_channel_length = args_info.RNAP_exit_channel_length_arg;
  GTV.lmin = args_info.lmin_flag;
  GTV.fpt  = args_info.fpt_flag;
  GTV.rect = args_info.rect_flag;
  GTV.dump_transcription_neighbors = args_info.dump_transcription_neighbors_flag;
  if (GSV.seconds_per_internal_time_unit <= 0.0) {
    fprintf(stderr,
            "Value of --kinfold_seconds_per_time_unit must be > 0 >%g<\n",
            GSV.seconds_per_internal_time_unit);
    exit(EXIT_FAILURE);
  }
  if (GSV.max_bubble_width < 0) {
    fprintf(stderr,
            "Value of --max_bubble_width must be >= 0 >%d<\n",
            GSV.max_bubble_width);
    exit(EXIT_FAILURE);
  }
  if (GSV.RNAP_exit_channel_length < 0) {
    fprintf(stderr,
            "Value of --RNAP_exit_channel_length must be >= 0 >%d<\n",
            GSV.RNAP_exit_channel_length);
    exit(EXIT_FAILURE);
  }
  if (GSV.transcription_elongation_rate_nt_per_s < 0.0) {
    fprintf(stderr,
            "Value of --transcription_elongation_rate must be >= 0 >%g<\n",
            GSV.transcription_elongation_rate_nt_per_s);
    exit(EXIT_FAILURE);
  }
  if (GSV.time_step <= 0.0) {
    fprintf(stderr, "Value of --time-step must be > 0 >%g<\n", GSV.time_step);
    exit(EXIT_FAILURE);
  }
  if ((GSV.transcription_elongation_rate_nt_per_s > 0.0) && (!args_info.glen_given))
    GSV.glen = 0;
  update_transcription_mode_defaults();
  cmdline_parser_free(&args_info);
}
/**/
static int process_options (int argc, char *argv[]) {
  int c, itmp;
  float ftmp;
  double dtmp;
  int option_index = 0;
  while ((c = getopt_long (argc, argv, "h",
			     long_options, &option_index)) != EOF) {
      switch (c) {
      case 0:
	if (strcmp(long_options[option_index].name,"dangle")==0) {
	  itmp = -1;
	  if (sscanf(optarg, "%d", &itmp) == 0)
	    usage(EXIT_FAILURE);
	  else if (itmp == 0 || itmp == 1 || itmp == 2 || itmp == 3)
	    GTV.dangle = itmp;
	  else {
	    fprintf(stderr, "Value of --dangle must be 0|1|2 >%d<\n", itmp);
	    usage (EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"Temp")==0) {
	  ftmp = -1.0;
	  if (sscanf(optarg, "%f", &ftmp) == 0)
	    usage(EXIT_FAILURE);
	  else if ( ftmp >= 0 )
	    GSV.Temp = ftmp;
	  else {
	    fprintf(stderr, "Value of --Temp must be >= 0 >%.2f<\n", ftmp);
	    usage (EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"Par")==0) {
	  if (sscanf(optarg, "%s", GAV.ParamFile) == 0)
	    usage(EXIT_FAILURE);
	  else {
	    GTV.Par = 1;
	  }
	}

	if (strcmp(long_options[option_index].name,"silent")==0) {
	  GTV.silent = 1;
	  GTV.verbose = 0;
	}

	if (strcmp(long_options[option_index].name,"verbose")==0) {
	  if (GTV.silent == 0)
	    GTV.verbose = 1;
	}

	if (strcmp(long_options[option_index].name,"seed")==0) {
	  if (sscanf(optarg, "%hu=%hu=%hu",
		     &GAV.subi[0], &GAV.subi[1], &GAV.subi[2]) == 0)
	    usage(EXIT_FAILURE);
	  else GTV.seed = 1;
	}

	/* laplace stuff */
	if (strcmp(long_options[option_index].name,"pbounds")==0) {
	  if (sscanf(optarg, "%g=%g=%g",
		     &GAV.phi_bounds[0],
		     &GAV.phi_bounds[1],
		     &GAV.phi_bounds[2]) == 0)
	    usage(EXIT_FAILURE);
	  else
	    GTV.phi = 1;
	  /* check if values are proper */
	  if (GAV.phi_bounds[0] > GAV.phi_bounds[2]
	      || GAV.phi_bounds[1] > GAV.phi_bounds[2]) {
	    fprintf(stderr,
		    "Unmet requirements for pbounds:\n"
		    "phi_min < phi_max && phi_inc < phi_max\n"
		    "phi_min: %g phi_inc: %g phi_max: %g\n",
		    GAV.phi_bounds[0], GAV.phi_bounds[1], GAV.phi_bounds[2]);
	    exit(EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"time")==0) {
	  dtmp = -1.0;
	  if (sscanf(optarg, "%lf", &dtmp) == 0)
	    usage(EXIT_FAILURE);
	  else if ( dtmp > 0 )
	    GSV.time = dtmp;
	  else {
	    fprintf(stderr, "Value of --time must be > 0 >%lf<\n", dtmp);
	    usage(EXIT_FAILURE);
	  }
	}

	/* laplace stuff */
	if (strcmp(long_options[option_index].name, "phi")==0) {
	  dtmp = -1.0;
	  if (sscanf(optarg, "%lf", &dtmp) == 0)
	    usage(EXIT_FAILURE);
	  else if ( dtmp > 0 ) {
	    GSV.phi = dtmp;
	    GTV.phi = 1;
	  }
	  else {
	    fprintf(stderr, "Value of --phi must be > 0 >%lf<\n", dtmp);
	    exit(EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"num")==0) {
	  itmp = -1;
	  if (sscanf(optarg, "%d", &itmp) == 0)
	    usage(EXIT_FAILURE);
	  else if ( itmp > 0 )
	    GSV.num = itmp;
	  else {
	    fprintf(stderr, "Value of --num must be > 0 >%d<\n", itmp);
	    usage(EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"log")==0)
	  if (sscanf(optarg, "%s", GAV.BaseName) == 0)
	    usage(EXIT_FAILURE);

	if (strcmp(long_options[option_index].name,"cut")==0)
	  if (sscanf(optarg, "%f", &GSV.cut) == 0)
	    usage(EXIT_FAILURE);

	if (strcmp(long_options[option_index].name,"grow")==0)
	  if (sscanf(optarg, "%lf", &GSV.grow) == 0)
	    usage(EXIT_FAILURE);

	if (strcmp(long_options[option_index].name,"glen")==0)
	  if (sscanf(optarg, "%d", &GSV.glen) == 0)
	    usage(EXIT_FAILURE);

	break;

      case 'h':
	usage (0);
      default:
	usage (EXIT_FAILURE);
      }
  }
  return optind;
}

/**/
static void ini_gtoggles(void) {
  GTV.Par = 0;
  GTV.seed = 0;
  GTV.dangle = 2;
  GTV.logML = 1;
  GTV.noLP = 0;
  GTV.noShift = 0;
  GTV.start = 0;
  GTV.stop = 0;
  GTV.silent = 0;
  GTV.phi = 0;
  GTV.verbose = 0;
  GTV.lmin = 0;
  GTV.fpt = 1;
  GTV.rect = 0;
  GTV.mc = 0;
  GTV.transcription = 0;
  GTV.max_bubble = 0;
  GTV.structure_time_only = 0;
  GTV.dump_transcription_neighbors = 0;
}

/**/
static void ini_gvars(void) {
  GSV.len = 0;
  GSV.num = 1;
  GSV.trajectory_index = 0;
  GSV.maxS = 99;
  GSV.cut = (float)DBL_MAX;
  GSV.Temp = 37.0;
  GSV.startE = 0.0;
  GSV.stopE = 0.0;
  GSV.currE = 0.0;
  GSV.time = 500.0;
  GSV.time_step = 1.0;
  GSV.phi = 1.0;
  GSV.simTime = 0.0;
  GSV.next_sample_time = 1.0;
  GSV.last_sample_time = 0.0;
  GSV.glen = 15;
  GSV.grow = 0.0;
  GSV.full_len = 0;
  GSV.tx_len = 0;
  GSV.RNAP_exit_channel_length = 20;
  GSV.channel_left = 0;
  GSV.channel_length = 0;
  GSV.transcription_elongation_rate_nt_per_s = 0.0;
  GSV.seconds_per_internal_time_unit = 1e-5;
  GSV.elongation_propensity_internal = 0.0;
  GSV.elongation_rate_current_nt_per_s = 0.0;
  GSV.elongation_propensity_internal_current = 0.0;
  GSV.max_bubble_width = 0;
  GSV.bubble_left = 0;
  GSV.bubble_width = 0;
  GSV.hybrid_left = 0;
  GSV.hybrid_len = 0;
}

/**/
static void ini_garrays(void) {
  GAV.ProgramName = (char *)calloc((size_t)256, sizeof(char));
  assert(GAV.ProgramName != NULL);
  GAV.ParamFile   = (char *)calloc((size_t)256, sizeof(char));
  assert(GAV.ParamFile != NULL);
  strcpy(GAV.ParamFile, "None");
  GAV.BaseName    = (char *)calloc((size_t)256, sizeof(char));
  assert(GAV.BaseName != NULL);
  strcpy(GAV.BaseName, "kinout");
  GAV.CSVOutputFile = (char *)calloc((size_t)256, sizeof(char));
  assert(GAV.CSVOutputFile != NULL);
  strcpy(GAV.CSVOutputFile, "");
  GAV.stopform = (char **)calloc(GSV.maxS + 1, sizeof(char *));
  assert(GAV.stopform != NULL);
  GAV.farbe = NULL;
  GAV.farbe_full = NULL;
  GAV.template_dna_full = NULL;
  GAV.startform = NULL;
  GAV.currform = NULL;
  GAV.prevform = NULL;
  GAV.phi_bounds[0] = 0.1;
  GAV.phi_bounds[1] = 0.1;
  GAV.phi_bounds[2] = 2.0;
}

int kinfold_transcription_enabled(void) {
  return GTV.transcription;
}

int kinfold_bubble_enabled(void) {
  return GTV.transcription && GTV.max_bubble;
}

const char *kinfold_full_sequence_state(void) {
  return full_sequence_state();
}

const char *kinfold_full_structure_state(void) {
  return full_structure_state();
}

const char *kinfold_cli_structure_state(void) {
  static char *buffer = NULL;
  static int  size = 0;
  const char  *state;
  int         channel_left;
  int         right;

  state = full_structure_state();
  update_exit_channel_state();

  if ((!GTV.transcription) || (!kinfold_bubble_enabled()) || (GSV.bubble_width <= 0))
    return state;

  if (GSV.full_len + 1 > size) {
    size = GSV.full_len + 1;
    buffer = (char *)realloc(buffer, (size_t)size);
  }

  memcpy(buffer, state, (size_t)GSV.full_len + 1);

  channel_left = exit_channel_left_state();
  if ((channel_left >= 0) &&
      (channel_left < GSV.full_len) &&
      (channel_left < GSV.bubble_left))
    buffer[channel_left] = 'c';

  if ((GSV.bubble_left >= 0) && (GSV.bubble_left < GSV.full_len))
    buffer[GSV.bubble_left] = 'b';

  right = GSV.bubble_left + GSV.bubble_width - 1;
  if ((right >= 0) && (right < GSV.full_len))
    buffer[right] = 'b';

  return buffer;
}

const char *kinfold_state_cache_key(void) {
  return state_cache_key();
}

double kinfold_current_hybrid_energy(void) {
  if (!kinfold_bubble_enabled() || (GSV.hybrid_len <= 0))
    return 0.0;

  return (double)hybrid_segment_energy_dcal(GSV.hybrid_left, GSV.hybrid_len) / 100.0;
}

int kinfold_current_hybrid_energy_dcal(void) {
  if (!kinfold_bubble_enabled() || (GSV.hybrid_len <= 0))
    return 0;

  return hybrid_segment_energy_dcal(GSV.hybrid_left, GSV.hybrid_len);
}

double kinfold_current_dna_duplex_energy(void) {
  if (!kinfold_bubble_enabled() || (GSV.bubble_width <= 0))
    return 0.0;

  return (double)dna_duplex_segment_energy_dcal(GSV.bubble_left, GSV.bubble_width) / 100.0;
}

int kinfold_current_dna_duplex_energy_dcal(void) {
  if (!kinfold_bubble_enabled() || (GSV.bubble_width <= 0))
    return 0;

  return dna_duplex_segment_energy_dcal(GSV.bubble_left, GSV.bubble_width);
}

double kinfold_current_total_energy(void) {
  return GSV.currE + kinfold_current_hybrid_energy();
}

double kinfold_current_elongation_rate(void) {
  if (!kinfold_transcription_enabled() || (GSV.transcription_elongation_rate_nt_per_s <= 0.0))
    return 0.0;

  update_exit_channel_state();
  GSV.elongation_rate_current_nt_per_s =
    GSV.transcription_elongation_rate_nt_per_s * pause_rate_scale_for_state();
  return GSV.elongation_rate_current_nt_per_s;
}

double kinfold_current_elongation_propensity(void) {
  GSV.elongation_propensity_internal_current =
    kinfold_current_elongation_rate() * GSV.seconds_per_internal_time_unit;
  return GSV.elongation_propensity_internal_current;
}

int kinfold_channel_left(void) {
  update_exit_channel_state();
  return GSV.channel_left;
}

int kinfold_channel_length(void) {
  update_exit_channel_state();
  return GSV.channel_length;
}

void kinfold_reset_sample_metrics(void) {
  memset(&current_sample_metrics, 0, sizeof(current_sample_metrics));
  sample_metrics_time = 0.0;
  sample_metrics_smoothed_hairpin = 0.0;
  sample_metrics_ready = 0;
  sample_metrics_initialized = 0;
}

void kinfold_metric_smoother_reset(KinfoldMetricSmootherState *state) {
  if (!state)
    return;

  state->last_time = 0.0;
  state->smoothed_terminal_hairpin = 0.0;
  state->initialized = 0;
}

void kinfold_compute_metrics_for_time(double                    time_value,
                                      KinfoldMetricSmootherState *state,
                                      KinfoldInformationalMetrics *metrics) {
  KinfoldInformationalMetrics raw_metrics;
  double dt, tau, alpha, smoothed, previous_smoothed;

  if (!metrics)
    return;

  memset(metrics, 0, sizeof(*metrics));

  memset(&raw_metrics, 0, sizeof(raw_metrics));
  compute_terminal_hairpin_metrics_raw(&raw_metrics);
  *metrics = raw_metrics;

  if ((!state) || (!state->initialized) || (time_value <= state->last_time)) {
    metrics->terminal_hairpin_dg_dt = 0.0;
    if (state) {
      state->smoothed_terminal_hairpin = raw_metrics.terminal_hairpin_energy;
      state->last_time = time_value;
      state->initialized = 1;
    }
  } else {
    dt = time_value - state->last_time;
    if (dt <= 0.0) {
      metrics->terminal_hairpin_dg_dt = 0.0;
    } else {
      tau = 5.0 * GSV.time_step;
      if (tau <= 0.0)
        tau = 1.0;
      alpha = 1.0 - exp(-dt / tau);
      previous_smoothed = state->smoothed_terminal_hairpin;
      smoothed = alpha * raw_metrics.terminal_hairpin_energy +
                 (1.0 - alpha) * previous_smoothed;
      metrics->terminal_hairpin_dg_dt =
        (smoothed - previous_smoothed) / dt;
      state->smoothed_terminal_hairpin = smoothed;
    }
    state->last_time = time_value;
  }
}

void kinfold_prepare_sample_metrics(double sample_time) {
  KinfoldMetricSmootherState state;

  if (sample_metrics_ready && (fabs(sample_time - sample_metrics_time) <= 1e-9))
    return;

  state.last_time = sample_metrics_time;
  state.smoothed_terminal_hairpin = sample_metrics_smoothed_hairpin;
  state.initialized = sample_metrics_initialized;
  kinfold_compute_metrics_for_time(sample_time, &state, &current_sample_metrics);
  sample_metrics_time = sample_time;
  sample_metrics_smoothed_hairpin = state.smoothed_terminal_hairpin;
  sample_metrics_ready = 1;
  sample_metrics_initialized = state.initialized;
}

void kinfold_get_sample_metrics(KinfoldInformationalMetrics *metrics) {
  if (!metrics)
    return;

  *metrics = current_sample_metrics;
}

int kinfold_stop_matches(const char *stopform) {
  int i;
  const char *current;

  if (!GTV.transcription)
    return (strcmp(stopform, GAV.currform) == 0);

  current = full_structure_state();
  for (i = 0; i < GSV.full_len; i++) {
    if (!is_transcribed_position(i) && (stopform[i] != '.'))
      return 0;
  }

  return (strcmp(stopform, current) == 0);
}

static void update_exit_channel_state(void) {
  int left, right;

  if (!kinfold_bubble_enabled() || (GSV.RNAP_exit_channel_length <= 0) || (GSV.bubble_left <= 0)) {
    GSV.channel_left = GSV.bubble_left;
    GSV.channel_length = 0;
    return;
  }

  left = GSV.bubble_left - GSV.RNAP_exit_channel_length;
  if (left < 0)
    left = 0;
  right = GSV.bubble_left - 1;
  if (right < left) {
    GSV.channel_left = GSV.bubble_left;
    GSV.channel_length = 0;
    return;
  }

  GSV.channel_left = left;
  GSV.channel_length = right - left + 1;
}

static const char *costring(const char *str) {
  static char* buffer=NULL;
  static int size=0;
  int n;
  if ((str==NULL) && (buffer)) {
    /* make it possible to free buffer */
    free(buffer);
    return NULL;
  }
  n=strlen(str);
  if (n>size) {
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

static void derive_template_dna(void) {
  int i;

  free(GAV.template_dna_full);
  GAV.template_dna_full = NULL;

  if (!GAV.farbe_full)
    return;

  GAV.template_dna_full = (char *)calloc((size_t)GSV.full_len + 1, sizeof(char));
  assert(GAV.template_dna_full != NULL);

  for (i = 0; i < GSV.full_len; i++)
    GAV.template_dna_full[i] = dna_complement(GAV.farbe_full[i]);
}

static void update_transcription_mode_defaults(void) {
  if ((GSV.grow > 0.0) && (GSV.transcription_elongation_rate_nt_per_s <= 0.0)) {
    GSV.transcription_elongation_rate_nt_per_s = 1.0 / GSV.grow;
    fprintf(stderr,
            "WARNING: --grow is deprecated; translating it to --transcription_elongation_rate=%g nt/s using --kinfold_seconds_per_time_unit=%g\n",
            GSV.transcription_elongation_rate_nt_per_s,
            GSV.seconds_per_internal_time_unit);
  } else if ((GSV.grow > 0.0) && (GSV.transcription_elongation_rate_nt_per_s > 0.0)) {
    fprintf(stderr,
            "ERROR: --grow may not be combined with --transcription_elongation_rate\n");
    exit(EXIT_FAILURE);
  }

  if (GSV.transcription_elongation_rate_nt_per_s > 0.0) {
    GTV.transcription = 1;
    GSV.elongation_propensity_internal =
      GSV.transcription_elongation_rate_nt_per_s * GSV.seconds_per_internal_time_unit;
    GSV.elongation_rate_current_nt_per_s = GSV.transcription_elongation_rate_nt_per_s;
    GSV.elongation_propensity_internal_current = GSV.elongation_propensity_internal;
  } else {
    GTV.transcription = 0;
    GSV.elongation_propensity_internal = 0.0;
    GSV.elongation_rate_current_nt_per_s = 0.0;
    GSV.elongation_propensity_internal_current = 0.0;
  }

  GTV.max_bubble = (GTV.transcription && (GSV.max_bubble_width > 0)) ? 1 : 0;
  if (!GTV.transcription) {
    GSV.bubble_left = 0;
    GSV.bubble_width = 0;
    GSV.hybrid_left = 0;
    GSV.hybrid_len = 0;
  }
  update_exit_channel_state();
}

static void normalize_transcription_start_structure(void) {
  int i;
  char *normalized;

  if (!GTV.transcription)
    return;

  if (!GTV.start) {
    memset(GAV.startform, '.', (size_t)GSV.tx_len);
    GAV.startform[GSV.tx_len] = '\0';
    return;
  }

  if ((int)strlen(GAV.startform) == GSV.tx_len)
    return;

  if ((int)strlen(GAV.startform) != GSV.full_len) {
    fprintf(stderr,
            "Transcription-aware mode requires start structures of length %d or %d\n",
            GSV.tx_len,
            GSV.full_len);
    exit(EXIT_FAILURE);
  }

  for (i = GSV.tx_len; i < GSV.full_len; i++) {
    if (GAV.startform[i] != '.') {
      fprintf(stderr,
              "Untranscribed suffix of start structure must contain only dots in transcription-aware mode\n");
      exit(EXIT_FAILURE);
    }
  }

  normalized = (char *)calloc((size_t)GSV.tx_len + 1, sizeof(char));
  assert(normalized != NULL);
  memcpy(normalized, GAV.startform, (size_t)GSV.tx_len);
  free(GAV.startform);
  GAV.startform = normalized;
}

static int is_transcribed_position(int idx0) {
  return (idx0 >= 0) && (idx0 < GSV.tx_len);
}

static char dna_complement(char rna_nt) {
  switch (rna_nt) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'U':
    case 'T': return 'A';
    default:  return 'N';
  }
}

static char coding_dna_nt(int idx0) {
  char nt;

  if (!GAV.farbe_full || (idx0 < 0) || (idx0 >= GSV.full_len))
    return 'N';

  nt = GAV.farbe_full[idx0];
  if (nt == 'U')
    return 'T';
  return nt;
}

static int hybrid_segment_energy_dcal(int left, int len) {
  int energy = 0;
  int i;

  if (!GAV.template_dna_full || (len <= 1))
    return 0;

  for (i = left; i < left + len - 1; i++) {
    const char *step = GAV.template_dna_full + i;

    switch (step[0]) {
      case 'T':
        switch (step[1]) {
          case 'T': energy += -70; break;
          case 'G': energy += -120; break;
          case 'C': energy += -150; break;
          case 'A': energy += -50; break;
        }
        break;
      case 'G':
        switch (step[1]) {
          case 'T': energy += -150; break;
          case 'G': energy += -170; break;
          case 'C': energy += -200; break;
          case 'A': energy += -140; break;
        }
        break;
      case 'C':
        switch (step[1]) {
          case 'T': energy += -130; break;
          case 'G': energy += -140; break;
          case 'C': energy += -230; break;
          case 'A': energy += -160; break;
        }
        break;
      case 'A':
        switch (step[1]) {
          case 'T': energy += -40; break;
          case 'G': energy += -40; break;
          case 'C': energy += -140; break;
          case 'A': energy += 20; break;
        }
        break;
    }
  }

  return energy;
}

static int hybrid_shifted_upstream_energy_dcal(int left, int len, int *supported) {
  int energy = 0;
  int i;

  if (supported)
    *supported = 0;

  if (!GAV.template_dna_full || !GAV.farbe_full || (len <= 1))
    return 0;

  if ((left < 0) || ((left + len - 1) >= GSV.full_len))
    return 0;

  for (i = left; i < left + len - 1; i++) {
    if (dna_complement(GAV.farbe_full[i + 1]) != GAV.template_dna_full[i])
      return 0;
  }

  if (supported)
    *supported = 1;

  for (i = left; i < left + len - 2; i++) {
    const char *step = GAV.template_dna_full + i;

    switch (step[0]) {
      case 'T':
        switch (step[1]) {
          case 'T': energy += -70; break;
          case 'G': energy += -120; break;
          case 'C': energy += -150; break;
          case 'A': energy += -50; break;
        }
        break;
      case 'G':
        switch (step[1]) {
          case 'T': energy += -150; break;
          case 'G': energy += -170; break;
          case 'C': energy += -200; break;
          case 'A': energy += -140; break;
        }
        break;
      case 'C':
        switch (step[1]) {
          case 'T': energy += -130; break;
          case 'G': energy += -140; break;
          case 'C': energy += -230; break;
          case 'A': energy += -160; break;
        }
        break;
      case 'A':
        switch (step[1]) {
          case 'T': energy += -40; break;
          case 'G': energy += -40; break;
          case 'C': energy += -140; break;
          case 'A': energy += 20; break;
        }
        break;
    }
  }

  return energy;
}

static int dna_duplex_segment_energy_dcal(int left, int len) {
  int energy = 0;
  int i;

  if (!GAV.farbe_full || (len <= 1))
    return 0;

  for (i = left; i < left + len - 1; i++) {
    const char *step = GAV.farbe_full + i;

    switch (step[0]) {
      case 'A':
        switch (step[1]) {
          case 'A': energy += -100; break; /* AA/TT */
          case 'G': energy += -130; break; /* AG/TC */
          case 'C': energy += -145; break; /* AC/TG */
          case 'U':
          case 'T': energy += -88; break;  /* AT/TA */
        }
        break;
      case 'G':
        switch (step[1]) {
          case 'A': energy += -130; break; /* GA/CT */
          case 'G': energy += -184; break; /* GG/CC */
          case 'C': energy += -224; break; /* GC/CG */
          case 'U':
          case 'T': energy += -144; break; /* GT/CA */
        }
        break;
      case 'C':
        switch (step[1]) {
          case 'A': energy += -145; break; /* CA/GT */
          case 'G': energy += -217; break; /* CG/GC */
          case 'C': energy += -184; break; /* CC/GG */
          case 'U':
          case 'T': energy += -128; break; /* CT/GA */
        }
        break;
      case 'U':
      case 'T':
        switch (step[1]) {
          case 'A': energy += -58; break;  /* TA/AT */
          case 'G': energy += -144; break; /* TG/AC */
          case 'C': energy += -128; break; /* TC/AG */
          case 'U':
          case 'T': energy += -100; break; /* TT/AA */
        }
        break;
    }
  }

  return energy;
}

static int eval_substructure_dcal(const char *sequence, const char *structure) {
  int i, n;
  char *seq_copy;

  if ((!sequence) || (!structure))
    return 0;

  n = (int)strlen(sequence);
  if ((n <= 0) || ((int)strlen(structure) != n))
    return 0;

  seq_copy = strdup(sequence);
  assert(seq_copy != NULL);
  for (i = 0; i < n; i++) {
    if (seq_copy[i] == 'T')
      seq_copy[i] = 'U';
  }

#if HAVE_LIBRNA_API3
  {
    float e;
    vrna_md_t md;
    vrna_fold_compound_t *fc;

    md = GAV.md;
    fc = vrna_fold_compound(seq_copy, &md, VRNA_OPTION_EVAL_ONLY);
    e = vrna_eval_structure(fc, structure);
    vrna_fold_compound_free(fc);
    free(seq_copy);
    return (int)(e * 100.0 + ((e < 0.0) ? -0.4 : 0.4));
  }
#else
  {
    float e;
    e = energy_of_structure(seq_copy, (char *)structure, 0);
    free(seq_copy);
    return (int)(e * 100.0 + ((e < 0.0) ? -0.4 : 0.4));
  }
#endif
}

static int build_pair_table(const char *structure, int len, int *pair_table, int *parent_open) {
  int *stack;
  int top, i;

  stack = (int *)calloc((size_t)len, sizeof(int));
  assert(stack != NULL);

  for (i = 0; i < len; i++) {
    pair_table[i] = -1;
    parent_open[i] = -1;
  }

  top = -1;
  for (i = 0; i < len; i++) {
    switch (structure[i]) {
      case '(':
        parent_open[i] = (top >= 0) ? stack[top] : -1;
        stack[++top] = i;
        break;
      case ')':
        if (top < 0) {
          free(stack);
          return 0;
        }
        pair_table[i] = stack[top];
        pair_table[stack[top]] = i;
        top--;
        break;
      default:
        break;
    }
  }

  free(stack);

  return (top == -1) ? 1 : 0;
}

static int build_substructure_with_outside_pairs_dotted(const char *sequence,
                                                        const char *structure,
                                                        int        left,
                                                        int        right,
                                                        char       **subseq,
                                                        char       **substruc) {
  int i;
  int sublen;
  int *pair_table, *parent_open;

  if ((!sequence) || (!structure) || (!subseq) || (!substruc) || (left < 0) || (right < left))
    return 0;

  sublen = right - left + 1;
  pair_table = (int *)calloc((size_t)GSV.len, sizeof(int));
  parent_open = (int *)calloc((size_t)GSV.len, sizeof(int));
  assert(pair_table != NULL);
  assert(parent_open != NULL);

  if (!build_pair_table(structure, GSV.len, pair_table, parent_open)) {
    free(pair_table);
    free(parent_open);
    return 0;
  }

  *subseq = (char *)calloc((size_t)sublen + 1, sizeof(char));
  *substruc = (char *)calloc((size_t)sublen + 1, sizeof(char));
  assert(*subseq != NULL);
  assert(*substruc != NULL);
  memcpy(*subseq, sequence + left, (size_t)sublen);
  memcpy(*substruc, structure + left, (size_t)sublen);

  for (i = left; i <= right; i++) {
    int partner = pair_table[i];

    if ((partner >= 0) && ((partner < left) || (partner > right)))
      (*substruc)[i - left] = '.';
  }

  free(pair_table);
  free(parent_open);
  return 1;
}

static int pair_crosses_existing(const int *pair_table, int len, int i, int j) {
  int a, b, k;

  if (i > j) {
    int tmp = i;
    i = j;
    j = tmp;
  }

  for (k = 0; k < len; k++) {
    if (pair_table[k] > k) {
      a = k;
      b = pair_table[k];
      if ((a == i) || (a == j) || (b == i) || (b == j))
        continue;
      if (((a < i) && (i < b) && (b < j)) ||
          ((i < a) && (a < j) && (j < b)))
        return 1;
    }
  }

  return 0;
}

static int position_in_exit_channel(int idx0) {
  return kinfold_bubble_enabled() &&
         (GSV.channel_length > 0) &&
         (idx0 >= GSV.channel_left) &&
         (idx0 < GSV.channel_left + GSV.channel_length);
}

static int exit_channel_left_state(void) {
  return (GSV.channel_length > 0) ? GSV.channel_left : -1;
}

static int exit_channel_right_state(void) {
  return (GSV.channel_length > 0) ? (GSV.channel_left + GSV.channel_length - 1) : -1;
}

static double terminal_component_energy_for_structure(const char *structure) {
  int *pair_table, *parent_open;
  int len, channel_left, channel_right;
  int left = -1, right = -1, i;
  char *subseq, *substruc;
  double energy = 0.0;

  if (!structure)
    return 0.0;

  update_exit_channel_state();
  len = GSV.len;
  channel_left = exit_channel_left_state();
  channel_right = exit_channel_right_state();
  if ((len <= 0) || (channel_left < 0) || (channel_right < channel_left))
    return 0.0;

  pair_table = (int *)calloc((size_t)len, sizeof(int));
  parent_open = (int *)calloc((size_t)len, sizeof(int));
  assert(pair_table != NULL);
  assert(parent_open != NULL);

  if (!build_pair_table(structure, len, pair_table, parent_open)) {
    free(pair_table);
    free(parent_open);
    return 0.0;
  }

  for (i = 0; i < len; i++) {
    if ((pair_table[i] > i) &&
        (parent_open[i] == -1) &&
        (i <= channel_right) &&
        (pair_table[i] >= channel_left) &&
        (pair_table[i] > right)) {
      left = i;
      right = pair_table[i];
    }
  }

  if ((left >= 0) && (right >= left)) {
    subseq = (char *)calloc((size_t)(right - left + 2), sizeof(char));
    substruc = (char *)calloc((size_t)(right - left + 2), sizeof(char));
    assert(subseq != NULL);
    assert(substruc != NULL);
    memcpy(subseq, GAV.farbe + left, (size_t)(right - left + 1));
    memcpy(substruc, structure + left, (size_t)(right - left + 1));
    energy = (double)eval_substructure_dcal(subseq, substruc) / 100.0;
    free(subseq);
    free(substruc);
  }

  free(pair_table);
  free(parent_open);
  return energy;
}

static int rewound_dna_duplex_energy_dcal(int released_hybrid_nts) {
  int rewound_len = 0;

  if (!kinfold_bubble_enabled() || (released_hybrid_nts < 0) || (GSV.bubble_width <= 0))
    return 0;

  if ((GSV.hybrid_len > 0) &&
      (GSV.hybrid_left >= GSV.bubble_left) &&
      (GSV.hybrid_left <= GSV.bubble_left + GSV.bubble_width))
    rewound_len = GSV.hybrid_left - GSV.bubble_left;

  rewound_len += released_hybrid_nts;
  if (rewound_len <= 1)
    return 0;
  if (rewound_len > GSV.bubble_width)
    rewound_len = GSV.bubble_width;

  return dna_duplex_segment_energy_dcal(GSV.bubble_left, rewound_len);
}

static double evaluate_rloop_terminal_plus_energy(int pair_count, int *found_out) {
  int *pair_table, *parent_open, *temp_pairs;
  char *candidate;
  int len, channel_left, channel_right;
  int i, j, target0, target1;
  int found = 0;
  double best = 0.0;

  if (found_out)
    *found_out = 0;

  if (!kinfold_bubble_enabled() || (pair_count <= 0) || (GSV.hybrid_len < pair_count))
    return 0.0;

  update_exit_channel_state();
  len = GSV.len;
  channel_left = exit_channel_left_state();
  channel_right = exit_channel_right_state();
  if ((len <= 0) || (channel_left < 0) || (channel_right < channel_left))
    return 0.0;

  target0 = GSV.hybrid_left;
  target1 = GSV.hybrid_left + 1;
  if ((target0 < 0) || (target0 >= len))
    return 0.0;
  if ((pair_count >= 2) && ((target1 < 0) || (target1 >= len)))
    return 0.0;
  if (GAV.currform[target0] != '.')
    return 0.0;
  if ((pair_count >= 2) && (GAV.currform[target1] != '.'))
    return 0.0;

  pair_table = (int *)calloc((size_t)len, sizeof(int));
  parent_open = (int *)calloc((size_t)len, sizeof(int));
  assert(pair_table != NULL);
  assert(parent_open != NULL);
  if (!build_pair_table(GAV.currform, len, pair_table, parent_open)) {
    free(pair_table);
    free(parent_open);
    return 0.0;
  }
  free(parent_open);

  candidate = strdup(GAV.currform);
  temp_pairs = (int *)calloc((size_t)len, sizeof(int));
  assert(candidate != NULL);
  assert(temp_pairs != NULL);

  for (i = channel_left; i <= channel_right; i++) {
    if (pair_table[i] >= 0)
      continue;
    if (!kinfold_can_pair_positions_raw(i, target0))
      continue;
    if (pair_crosses_existing(pair_table, len, i, target0))
      continue;

    memcpy(temp_pairs, pair_table, sizeof(int) * (size_t)len);
    temp_pairs[i] = target0;
    temp_pairs[target0] = i;
    if (i < target0) {
      candidate[i] = '(';
      candidate[target0] = ')';
    } else {
      candidate[target0] = '(';
      candidate[i] = ')';
    }

    if (pair_count == 1) {
      double energy = terminal_component_energy_for_structure(candidate);
      if (!found || (energy < best)) {
        best = energy;
        found = 1;
      }
    } else {
      for (j = channel_left; j <= channel_right; j++) {
        double energy;

        if ((j == i) || (pair_table[j] >= 0))
          continue;
        if (!kinfold_can_pair_positions_raw(j, target1))
          continue;
        if (pair_crosses_existing(temp_pairs, len, j, target1))
          continue;

        temp_pairs[j] = target1;
        temp_pairs[target1] = j;
        if (j < target1) {
          candidate[j] = '(';
          candidate[target1] = ')';
        } else {
          candidate[target1] = '(';
          candidate[j] = ')';
        }

        energy = terminal_component_energy_for_structure(candidate);
        if (!found || (energy < best)) {
          best = energy;
          found = 1;
        }

        temp_pairs[j] = -1;
        temp_pairs[target1] = -1;
        candidate[j] = '.';
        candidate[target1] = '.';
      }
    }

    temp_pairs[i] = -1;
    temp_pairs[target0] = -1;
    candidate[i] = '.';
    candidate[target0] = '.';
  }

  free(candidate);
  free(temp_pairs);
  free(pair_table);
  if (found_out)
    *found_out = found;
  return found ? best : 0.0;
}

static double pause_rate_scale_for_state(void) {
  static const int offsets[16] = {-12, -11, -10, -9, -8, -7, -6, -5,
                                  -4, -3, -2, -1, 1, 2, 3, 4};
  static const char consensus[16] = {'G','G','c','a','t','a','a','t',
                                     'T','G','?','G','G','C','c','g'};
  double scale = 1.0;
  int site_idx, i, idx;
  char dna_nt;
  char expect;

  if (!kinfold_transcription_enabled() || (GSV.transcription_elongation_rate_nt_per_s <= 0.0))
    return 0.0;

  site_idx = GSV.tx_len - 1;
  if (site_idx < 0)
    return 1.0;

  for (i = 0; i < 16; i++) {
    idx = site_idx + offsets[i];
    if ((idx < 0) || (idx >= GSV.full_len))
      continue;

    dna_nt = coding_dna_nt(idx);
    expect = consensus[i];

    if (expect == '?') {
      if ((dna_nt == 'C') || (dna_nt == 'T'))
        scale -= 0.075;
      continue;
    }

    if (toupper((unsigned char)expect) == dna_nt) {
      if (isupper((unsigned char)expect))
        scale -= 0.075;
      else
        scale -= 0.0375;
    }
  }

  if (scale < 0.1)
    scale = 0.1;

  return scale;
}

static int evaluate_invasion_energy(int invade_count) {
  int *pair_table, *parent_open;
  int hybrid_end, idx0, idx1, best = 0, found = 0;
  char *candidate;
  int i, j;

  if (!kinfold_bubble_enabled() || (invade_count <= 0) || (GSV.hybrid_len <= 0))
    return 0;

  hybrid_end = GSV.hybrid_left + GSV.hybrid_len - 1;
  if (hybrid_end < 0 || hybrid_end >= GSV.len)
    return 0;

  pair_table = (int *)calloc((size_t)GSV.len, sizeof(int));
  parent_open = (int *)calloc((size_t)GSV.len, sizeof(int));
  assert(pair_table != NULL);
  assert(parent_open != NULL);
  if (!build_pair_table(GAV.currform, GSV.len, pair_table, parent_open)) {
    free(pair_table);
    free(parent_open);
    return 0;
  }

  candidate = strdup(GAV.currform);
  assert(candidate != NULL);

  idx0 = hybrid_end;
  for (i = 0; i < GSV.hybrid_left; i++) {
    int candidate_found = 0;
    int eval_dcal;

    if (!kinfold_is_unpaired_position(i))
      continue;
    if (!kinfold_can_pair_positions_raw(i, idx0))
      continue;
    if (pair_crosses_existing(pair_table, GSV.len, i, idx0))
      continue;

    if (i < idx0) {
      candidate[i] = '(';
      candidate[idx0] = ')';
    } else {
      candidate[idx0] = '(';
      candidate[i] = ')';
    }

    if (invade_count >= 2 && (GSV.hybrid_len >= 2)) {
      idx1 = hybrid_end - 1;
      for (j = 0; j < GSV.hybrid_left; j++) {
        if (j == i)
          continue;
        if (!(i < j && j < idx1))
          continue;
        if (!kinfold_is_unpaired_position(j))
          continue;
        if (!kinfold_can_pair_positions_raw(j, idx1))
          continue;
        if (pair_crosses_existing(pair_table, GSV.len, j, idx1))
          continue;

        if (j < idx1) {
          candidate[j] = '(';
          candidate[idx1] = ')';
        } else {
          candidate[idx1] = '(';
          candidate[j] = ')';
        }
        eval_dcal = kinfold_eval_structure_dcal(candidate);
        candidate_found = 1;
        found = 1;
        if (!best || (eval_dcal < best))
          best = eval_dcal;
        candidate[j] = '.';
        candidate[idx1] = '.';
      }
    }

    if (invade_count == 1) {
      eval_dcal = kinfold_eval_structure_dcal(candidate);
      candidate_found = 1;
      found = 1;
      if (!best || (eval_dcal < best))
        best = eval_dcal;
    } else if (invade_count >= 2) {
      /* also allow only one of the two nts to invade if that is the minimum */
      eval_dcal = kinfold_eval_structure_dcal(candidate);
      candidate_found = 1;
      found = 1;
      if (!best || (eval_dcal < best))
        best = eval_dcal;
    }

    (void)candidate_found;
    candidate[i] = '.';
    candidate[idx0] = '.';
  }

  free(candidate);
  free(pair_table);
  free(parent_open);

  return found ? best : 0;
}

static int compute_terminal_hairpin_metrics_raw(KinfoldInformationalMetrics *metrics) {
  int *pair_table, *parent_open;
  int len, bubble_left, left, right, i;
  int channel_left, channel_right;
  double current_hybrid_energy, minus1_hybrid_energy, minus2_hybrid_energy;
  double current_rewound_dna_energy, plus1_rewound_dna_energy, plus2_rewound_dna_energy;
  double terminal_plus1_energy, terminal_plus2_energy;
  int terminal_plus1_found, terminal_plus2_found;
  char *subseq, *substruc, *plus1_structure, *plus1_seq, *plus1_struc;
  int supported_shift;
  int shifted_dcal;

  if (!metrics)
    return 0;

  memset(metrics, 0, sizeof(*metrics));

  if (!kinfold_bubble_enabled())
    return 0;

  update_exit_channel_state();
  metrics->elongation_rate_current = kinfold_current_elongation_rate();
  current_hybrid_energy = kinfold_current_hybrid_energy();
  minus1_hybrid_energy = 0.0;
  minus2_hybrid_energy = 0.0;
  if (GSV.hybrid_len > 1)
    minus1_hybrid_energy =
      (double)hybrid_segment_energy_dcal(GSV.hybrid_left + 1, GSV.hybrid_len - 1) / 100.0;
  if (GSV.hybrid_len > 2)
    minus2_hybrid_energy =
      (double)hybrid_segment_energy_dcal(GSV.hybrid_left + 2, GSV.hybrid_len - 2) / 100.0;
  metrics->rna_dna_minus1_ddg = (GSV.hybrid_len > 0) ?
                                (minus1_hybrid_energy - current_hybrid_energy) :
                                0.0;
  shifted_dcal = hybrid_shifted_upstream_energy_dcal(GSV.hybrid_left, GSV.hybrid_len, &supported_shift);
  metrics->rna_dna_1nt_shift_energy = supported_shift ? ((double)shifted_dcal / 100.0) : 0.0;
  metrics->invasion_1nt_energy = (double)evaluate_invasion_energy(1) / 100.0;
  metrics->invasion_2nt_energy = (double)evaluate_invasion_energy(2) / 100.0;
  current_rewound_dna_energy = (double)rewound_dna_duplex_energy_dcal(0) / 100.0;
  plus1_rewound_dna_energy = (GSV.hybrid_len >= 1) ?
                             ((double)rewound_dna_duplex_energy_dcal(1) / 100.0) :
                             0.0;
  plus2_rewound_dna_energy = (GSV.hybrid_len >= 2) ?
                             ((double)rewound_dna_duplex_energy_dcal(2) / 100.0) :
                             0.0;
  terminal_plus1_energy = evaluate_rloop_terminal_plus_energy(1, &terminal_plus1_found);
  terminal_plus2_energy = evaluate_rloop_terminal_plus_energy(2, &terminal_plus2_found);

  len = GSV.len;
  bubble_left = GSV.bubble_left;
  channel_left = exit_channel_left_state();
  channel_right = exit_channel_right_state();
  if ((len <= 0) || (bubble_left <= 0))
    return 0;

  pair_table = (int *)calloc((size_t)len, sizeof(int));
  parent_open = (int *)calloc((size_t)len, sizeof(int));
  assert(pair_table != NULL);
  assert(parent_open != NULL);

  if (!build_pair_table(GAV.currform, len, pair_table, parent_open)) {
    free(pair_table);
    free(parent_open);
    return 0;
  }

  if (channel_left >= 0) {
    for (i = 0; i < len; i++) {
      if (pair_table[i] > i) {
        if (position_in_exit_channel(i) || position_in_exit_channel(pair_table[i]))
          metrics->exit_channel_bp_count++;
      }
    }
  }

  left = -1;
  right = -1;
  for (i = 0; i < bubble_left; i++) {
    if ((pair_table[i] > i) &&
        (pair_table[i] < bubble_left) &&
        (parent_open[i] == -1) &&
        (pair_table[i] > right)) {
      left = i;
      right = pair_table[i];
    }
  }

  if ((left < 0) || (right < left)) {
    goto exit_channel_metrics;
  }

  subseq = (char *)calloc((size_t)(right - left + 2), sizeof(char));
  substruc = (char *)calloc((size_t)(right - left + 2), sizeof(char));
  assert(subseq != NULL);
  assert(substruc != NULL);
  memcpy(subseq, GAV.farbe + left, (size_t)(right - left + 1));
  memcpy(substruc, GAV.currform + left, (size_t)(right - left + 1));
  metrics->terminal_hairpin_energy = (double)eval_substructure_dcal(subseq, substruc) / 100.0;
  free(subseq);
  free(substruc);

  if ((GSV.hybrid_len > 0) &&
      (left > 0) &&
      (GSV.hybrid_left >= 0) &&
      (GSV.hybrid_left < len) &&
      ((left - 1) < GSV.bubble_left) &&
      (pair_table[left - 1] < 0) &&
      (GAV.currform[GSV.hybrid_left] == '.') &&
      kinfold_can_pair_positions_raw(left - 1, GSV.hybrid_left)) {
    size_t plus1_len;

    plus1_structure = strdup(GAV.currform);
    assert(plus1_structure != NULL);
    plus1_structure[left - 1] = '(';
    plus1_structure[GSV.hybrid_left] = ')';

    plus1_len = (size_t)(GSV.hybrid_left - (left - 1) + 1);
    plus1_seq = (char *)calloc(plus1_len + 1, sizeof(char));
    plus1_struc = (char *)calloc(plus1_len + 1, sizeof(char));
    assert(plus1_seq != NULL);
    assert(plus1_struc != NULL);
    memcpy(plus1_seq, GAV.farbe + left - 1, plus1_len);
    memcpy(plus1_struc, plus1_structure + left - 1, plus1_len);
    metrics->terminal_hairpin_plus1_energy =
      (double)eval_substructure_dcal(plus1_seq, plus1_struc) / 100.0;
    metrics->plus1_ddg =
      metrics->terminal_hairpin_plus1_energy -
      metrics->terminal_hairpin_energy -
      metrics->rna_dna_minus1_ddg;
    free(plus1_seq);
    free(plus1_struc);
    free(plus1_structure);
  }

exit_channel_metrics:
  left = -1;
  right = -1;
  if (channel_left >= 0) {
    for (i = 0; i < len; i++) {
      if ((pair_table[i] > i) &&
          (parent_open[i] == -1) &&
          (i <= channel_right) &&
          (pair_table[i] >= channel_left) &&
          (pair_table[i] > right)) {
        left = i;
        right = pair_table[i];
      }
    }
  }

  if ((left >= 0) && (right >= left)) {
    subseq = (char *)calloc((size_t)(right - left + 2), sizeof(char));
    substruc = (char *)calloc((size_t)(right - left + 2), sizeof(char));
    assert(subseq != NULL);
    assert(substruc != NULL);
    memcpy(subseq, GAV.farbe + left, (size_t)(right - left + 1));
    memcpy(substruc, GAV.currform + left, (size_t)(right - left + 1));
    metrics->terminal_energy = (double)eval_substructure_dcal(subseq, substruc) / 100.0;
    free(subseq);
    free(substruc);
  }

  if ((channel_left >= 0) && (channel_right >= channel_left) &&
      build_substructure_with_outside_pairs_dotted(GAV.farbe,
                                                   GAV.currform,
                                                   channel_left,
                                                   channel_right,
                                                   &subseq,
                                                   &substruc)) {
    metrics->folding_inside_channel_energy =
      (double)eval_substructure_dcal(subseq, substruc) / 100.0;
    free(subseq);
    free(substruc);
  }

  if (GSV.hybrid_len >= 1 && terminal_plus1_found) {
    metrics->rloop_1nt_energy =
      terminal_plus1_energy +
      minus1_hybrid_energy +
      plus1_rewound_dna_energy -
      metrics->terminal_energy -
      current_hybrid_energy -
      current_rewound_dna_energy;
  }

  if (GSV.hybrid_len >= 2 && terminal_plus2_found) {
    metrics->rloop_2nt_energy =
      terminal_plus2_energy +
      minus2_hybrid_energy +
      plus2_rewound_dna_energy -
      metrics->terminal_energy -
      current_hybrid_energy -
      current_rewound_dna_energy;
  }

  free(pair_table);
  free(parent_open);

  return 1;
}

static const char *full_sequence_state(void) {
  static char *buffer = NULL;
  static int size = 0;
  int i;

  if (GSV.full_len + 1 > size) {
    size = GSV.full_len + 1;
    buffer = (char *)realloc(buffer, (size_t)size);
  }

  for (i = 0; i < GSV.full_len; i++)
    buffer[i] = is_transcribed_position(i) ? GAV.farbe_full[i] : 'X';
  buffer[GSV.full_len] = '\0';

  return buffer;
}

static const char *full_structure_state(void) {
  return full_structure_state_from(GAV.currform);
}

static const char *full_structure_state_from(const char *structure) {
  static char *buffer = NULL;
  static int size = 0;
  int i;

  if (GSV.full_len + 1 > size) {
    size = GSV.full_len + 1;
    buffer = (char *)realloc(buffer, (size_t)size);
  }

  for (i = 0; i < GSV.full_len; i++)
    buffer[i] = is_transcribed_position(i) ? structure[i] : 'x';
  buffer[GSV.full_len] = '\0';

  return buffer;
}

static const char *state_cache_key(void) {
  static char *buffer = NULL;
  static int size = 0;
  int needed;

  needed = (GSV.full_len * 2) + 128;
  if (needed > size) {
    size = needed;
    buffer = (char *)realloc(buffer, (size_t)size);
  }

  snprintf(buffer,
           (size_t)size,
           "%d|%d|%d|%d|%d|%s|%s",
           GSV.tx_len,
           GSV.bubble_left,
           GSV.bubble_width,
           GSV.hybrid_left,
           GSV.hybrid_len,
           full_sequence_state(),
           full_structure_state());

  return buffer;
}
/* End of file */
