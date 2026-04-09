#include <Python.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

extern "C" {
#include "../../src/Kinfold/kinfold_api.h"
}

namespace {

struct KinfoldPyOptions {
  std::string               sequence;
  bool                      has_start;
  std::string               start;
  std::vector<std::string>  stops;
  int                       dangle;
  double                    temp;
  bool                      has_par;
  std::string               par;
  bool                      logML;
  bool                      noShift;
  bool                      noLP;
  bool                      has_seed;
  std::string               seed;
  double                    time;
  double                    time_step;
  int                       num;
  bool                      met;
  bool                      fpt;
  bool                      rect;
  bool                      has_glen;
  int                       glen;
  double                    transcription_elongation_rate;
  double                    kinfold_seconds_per_time_unit;
  int                       max_bubble_width;
  bool                      has_phi;
  double                    phi;
  bool                      has_pbounds;
  std::string               pbounds;
  std::string               log;
  bool                      silent;
  bool                      verbose;
  bool                      lmin;
  double                    cut;
};


static bool
pyobject_to_string(PyObject     *obj,
                   std::string  &out,
                   std::string  &error)
{
  PyObject *utf8;

  if (!PyUnicode_Check(obj)) {
    error = "Expected a string value";
    return false;
  }

  utf8 = PyUnicode_AsUTF8String(obj);
  if (!utf8) {
    error = "Failed to encode string argument";
    return false;
  }

  out.assign(PyBytes_AS_STRING(utf8), (size_t)PyBytes_GET_SIZE(utf8));
  Py_DECREF(utf8);
  return true;
}


static bool
parse_optional_string(PyObject     *obj,
                      bool         &present,
                      std::string  &out,
                      std::string  &error)
{
  present = false;
  out.clear();

  if ((obj == NULL) || (obj == Py_None))
    return true;

  present = true;
  return pyobject_to_string(obj, out, error);
}


static bool
parse_optional_double(PyObject *obj,
                      bool     &present,
                      double   &value,
                      std::string &error)
{
  if ((obj == NULL) || (obj == Py_None)) {
    present = false;
    value = 0.0;
    return true;
  }

  present = true;
  value = PyFloat_AsDouble(obj);
  if (PyErr_Occurred()) {
    PyErr_Clear();
    error = "Expected a float-compatible value";
    return false;
  }

  return true;
}


static bool
parse_seed_string(PyObject     *obj,
                  bool         &present,
                  std::string  &out,
                  std::string  &error)
{
  long values[3];
  Py_ssize_t i;

  present = false;
  out.clear();

  if ((obj == NULL) || (obj == Py_None))
    return true;

  if (PyUnicode_Check(obj))
    return parse_optional_string(obj, present, out, error);

  if ((!PyTuple_Check(obj)) && (!PyList_Check(obj))) {
    error = "seed must be None, a string, or a 3-item tuple/list";
    return false;
  }

  if (PySequence_Size(obj) != 3) {
    error = "seed must contain exactly 3 integers";
    return false;
  }

  for (i = 0; i < 3; i++) {
    PyObject *item = PySequence_GetItem(obj, i);
    if (!item) {
      error = "Failed to read seed item";
      return false;
    }
    values[i] = PyLong_AsLong(item);
    Py_DECREF(item);
    if (PyErr_Occurred()) {
      PyErr_Clear();
      error = "seed items must be integers";
      return false;
    }
  }

  std::ostringstream ss;
  ss << values[0] << "=" << values[1] << "=" << values[2];
  out = ss.str();
  present = true;
  return true;
}


static bool
parse_pbounds_string(PyObject     *obj,
                     bool         &present,
                     std::string  &out,
                     std::string  &error)
{
  double values[3];
  Py_ssize_t i;

  present = false;
  out.clear();

  if ((obj == NULL) || (obj == Py_None))
    return true;

  if (PyUnicode_Check(obj))
    return parse_optional_string(obj, present, out, error);

  if ((!PyTuple_Check(obj)) && (!PyList_Check(obj))) {
    error = "pbounds must be None, a string, or a 3-item tuple/list";
    return false;
  }

  if (PySequence_Size(obj) != 3) {
    error = "pbounds must contain exactly 3 numbers";
    return false;
  }

  for (i = 0; i < 3; i++) {
    PyObject *item = PySequence_GetItem(obj, i);
    if (!item) {
      error = "Failed to read pbounds item";
      return false;
    }
    values[i] = PyFloat_AsDouble(item);
    Py_DECREF(item);
    if (PyErr_Occurred()) {
      PyErr_Clear();
      error = "pbounds items must be numeric";
      return false;
    }
  }

  std::ostringstream ss;
  ss << values[0] << "=" << values[1] << "=" << values[2];
  out = ss.str();
  present = true;
  return true;
}


static bool
parse_stop_list(PyObject                      *obj,
                std::vector<std::string>      &out,
                std::string                   &error)
{
  Py_ssize_t i, n;

  out.clear();

  if ((obj == NULL) || (obj == Py_None))
    return true;

  if (PyUnicode_Check(obj)) {
    std::string tmp;
    if (!pyobject_to_string(obj, tmp, error))
      return false;
    out.push_back(tmp);
    return true;
  }

  if (!PySequence_Check(obj)) {
    error = "stop must be None, a string, or a sequence of strings";
    return false;
  }

  n = PySequence_Size(obj);
  for (i = 0; i < n; i++) {
    PyObject    *item = PySequence_GetItem(obj, i);
    std::string value;

    if (!item) {
      error = "Failed to read stop item";
      return false;
    }
    if (!pyobject_to_string(item, value, error)) {
      Py_DECREF(item);
      return false;
    }
    Py_DECREF(item);
    out.push_back(value);
  }

  return true;
}


static bool
validate_options(const KinfoldPyOptions  &opts,
                 std::string             &error)
{
  size_t seq_len, i;
  int effective_glen;

  seq_len = opts.sequence.size();
  if (seq_len == 0) {
    error = "sequence must not be empty";
    return false;
  }

  if ((opts.transcription_elongation_rate > 0.0) &&
      (opts.sequence.find('&') != std::string::npos)) {
    error = "transcription-aware mode does not support cut-point/cofold sequences";
    return false;
  }

  if (opts.num <= 0) {
    error = "num must be > 0";
    return false;
  }

  if (opts.time_step <= 0.0) {
    error = "time_step must be > 0";
    return false;
  }

  effective_glen = opts.has_glen ? opts.glen : ((opts.transcription_elongation_rate > 0.0) ? 0 : 15);
  if (effective_glen < 0) {
    error = "glen must be >= 0";
    return false;
  }

  if (opts.has_start) {
    if (opts.transcription_elongation_rate > 0.0) {
      if ((opts.start.size() != seq_len) &&
          ((int)opts.start.size() != effective_glen)) {
        error = "in transcription mode, start must match the full sequence length or the initial transcribed length";
        return false;
      }
    } else if (opts.start.size() != seq_len) {
      error = "start must match the sequence length";
      return false;
    }
  }

  for (i = 0; i < opts.stops.size(); i++) {
    if (opts.stops[i].size() != seq_len) {
      error = "every stop structure must match the sequence length";
      return false;
    }
  }

  if (opts.kinfold_seconds_per_time_unit <= 0.0) {
    error = "kinfold_seconds_per_time_unit must be > 0";
    return false;
  }

  if (opts.max_bubble_width < 0) {
    error = "max_bubble_width must be >= 0";
    return false;
  }

  return true;
}


static std::string
format_double(double value)
{
  std::ostringstream ss;
  ss.precision(17);
  ss << value;
  return ss.str();
}


static void
append_flag_if(std::vector<std::string>  &argv,
               bool                      condition,
               const std::string         &flag)
{
  if (condition)
    argv.push_back(flag);
}


static std::vector<std::string>
build_cli_args(const KinfoldPyOptions &opts)
{
  std::vector<std::string> argv;

  argv.push_back("kinfold");
  argv.push_back("--dangle=" + std::to_string(opts.dangle));
  argv.push_back("--Temp=" + format_double(opts.temp));
  append_flag_if(argv, opts.has_par, "--Par=" + opts.par);
  append_flag_if(argv, !opts.logML, "--logML");
  append_flag_if(argv, opts.noShift, "--noShift");
  append_flag_if(argv, opts.noLP, "--noLP");
  if (opts.has_seed)
    argv.push_back("--seed=" + opts.seed);
  argv.push_back("--time=" + format_double(opts.time));
  argv.push_back("--time-step=" + format_double(opts.time_step));
  argv.push_back("--num=" + std::to_string(opts.num));
  append_flag_if(argv, opts.has_start, "--start");
  append_flag_if(argv, !opts.stops.empty(), "--stop");
  append_flag_if(argv, opts.met, "--met");
  append_flag_if(argv, !opts.fpt, "--fpt");
  append_flag_if(argv, opts.rect, "--rect");
  if (opts.has_glen)
    argv.push_back("--glen=" + std::to_string(opts.glen));
  if (opts.transcription_elongation_rate > 0.0)
    argv.push_back("--transcription_elongation_rate=" + format_double(opts.transcription_elongation_rate));
  if ((opts.transcription_elongation_rate > 0.0) || (opts.kinfold_seconds_per_time_unit != 1e-5))
    argv.push_back("--kinfold_seconds_per_time_unit=" + format_double(opts.kinfold_seconds_per_time_unit));
  if (opts.max_bubble_width > 0)
    argv.push_back("--max_bubble_width=" + std::to_string(opts.max_bubble_width));
  if (opts.has_phi)
    argv.push_back("--phi=" + format_double(opts.phi));
  if (opts.has_pbounds)
    argv.push_back("--pbounds=" + opts.pbounds);
  if (!opts.log.empty())
    argv.push_back("--log=" + opts.log);
  append_flag_if(argv, opts.silent, "--silent");
  append_flag_if(argv, opts.verbose, "--verbose");
  append_flag_if(argv, opts.lmin, "--lmin");
  if (std::isfinite(opts.cut))
    argv.push_back("--cut=" + format_double(opts.cut));

  return argv;
}


static std::string
build_input_blob(const KinfoldPyOptions &opts)
{
  std::string blob;
  size_t      i;

  blob.append(opts.sequence);
  blob.push_back('\n');
  if (opts.has_start) {
    blob.append(opts.start);
    blob.push_back('\n');
  }
  for (i = 0; i < opts.stops.size(); i++) {
    blob.append(opts.stops[i]);
    blob.push_back('\n');
  }

  return blob;
}


static bool
run_kinfold_capture(const KinfoldPyOptions  &opts,
                    kinfold_result_t        &result,
                    std::string             &error)
{
  std::vector<std::string>  args;
  std::vector<char *>       argv;
  std::string               input_blob;
  FILE                      *fp;
  int                       status;
  size_t                    i;

  if (!validate_options(opts, error))
    return false;

  args = build_cli_args(opts);
  argv.reserve(args.size());
  for (i = 0; i < args.size(); i++)
    argv.push_back(const_cast<char *>(args[i].c_str()));

  input_blob = build_input_blob(opts);
  fp = tmpfile();
  if (!fp) {
    error = "Failed to create a temporary input stream for Kinfold";
    return false;
  }

  if (!input_blob.empty()) {
    if (fwrite(input_blob.data(), 1, input_blob.size(), fp) != input_blob.size()) {
      fclose(fp);
      error = "Failed to populate the temporary Kinfold input stream";
      return false;
    }
  }
  rewind(fp);

  kinfold_result_init(&result);
  Py_BEGIN_ALLOW_THREADS
  status = kinfold_run_with_args((int)argv.size(), argv.data(), fp, &result, 1);
  Py_END_ALLOW_THREADS
  fclose(fp);

  if (status != 0) {
    kinfold_result_free(&result);
    error = "Kinfold execution failed";
    return false;
  }

  return true;
}


static PyObject *
build_parameters_dict(const KinfoldPyOptions &opts)
{
  PyObject *params = PyDict_New();
  PyObject *stop_list = PyList_New((Py_ssize_t)opts.stops.size());
  size_t i;

  PyDict_SetItemString(params, "dangle", PyLong_FromLong(opts.dangle));
  PyDict_SetItemString(params, "Temp", PyFloat_FromDouble(opts.temp));
  if (opts.has_par)
    PyDict_SetItemString(params, "Par", PyUnicode_FromString(opts.par.c_str()));
  else {
    Py_INCREF(Py_None);
    PyDict_SetItemString(params, "Par", Py_None);
  }
  PyDict_SetItemString(params, "logML", PyBool_FromLong(opts.logML ? 1 : 0));
  PyDict_SetItemString(params, "noShift", PyBool_FromLong(opts.noShift ? 1 : 0));
  PyDict_SetItemString(params, "noLP", PyBool_FromLong(opts.noLP ? 1 : 0));
  if (opts.has_seed)
    PyDict_SetItemString(params, "seed", PyUnicode_FromString(opts.seed.c_str()));
  else {
    Py_INCREF(Py_None);
    PyDict_SetItemString(params, "seed", Py_None);
  }
  PyDict_SetItemString(params, "time", PyFloat_FromDouble(opts.time));
  PyDict_SetItemString(params, "time_step", PyFloat_FromDouble(opts.time_step));
  PyDict_SetItemString(params, "num", PyLong_FromLong(opts.num));
  PyDict_SetItemString(params, "met", PyBool_FromLong(opts.met ? 1 : 0));
  PyDict_SetItemString(params, "fpt", PyBool_FromLong(opts.fpt ? 1 : 0));
  PyDict_SetItemString(params, "rect", PyBool_FromLong(opts.rect ? 1 : 0));
  if (opts.has_glen)
    PyDict_SetItemString(params, "glen", PyLong_FromLong(opts.glen));
  else {
    Py_INCREF(Py_None);
    PyDict_SetItemString(params, "glen", Py_None);
  }
  PyDict_SetItemString(params, "transcription_elongation_rate", PyFloat_FromDouble(opts.transcription_elongation_rate));
  PyDict_SetItemString(params, "kinfold_seconds_per_time_unit", PyFloat_FromDouble(opts.kinfold_seconds_per_time_unit));
  PyDict_SetItemString(params, "max_bubble_width", PyLong_FromLong(opts.max_bubble_width));
  if (opts.has_phi)
    PyDict_SetItemString(params, "phi", PyFloat_FromDouble(opts.phi));
  else {
    Py_INCREF(Py_None);
    PyDict_SetItemString(params, "phi", Py_None);
  }
  if (opts.has_pbounds)
    PyDict_SetItemString(params, "pbounds", PyUnicode_FromString(opts.pbounds.c_str()));
  else {
    Py_INCREF(Py_None);
    PyDict_SetItemString(params, "pbounds", Py_None);
  }
  PyDict_SetItemString(params, "log", PyUnicode_FromString(opts.log.c_str()));
  PyDict_SetItemString(params, "silent", PyBool_FromLong(opts.silent ? 1 : 0));
  PyDict_SetItemString(params, "verbose", PyBool_FromLong(opts.verbose ? 1 : 0));
  PyDict_SetItemString(params, "lmin", PyBool_FromLong(opts.lmin ? 1 : 0));
  PyDict_SetItemString(params, "cut", PyFloat_FromDouble(opts.cut));

  for (i = 0; i < opts.stops.size(); i++)
    PyList_SET_ITEM(stop_list, (Py_ssize_t)i, PyUnicode_FromString(opts.stops[i].c_str()));
  PyDict_SetItemString(params, "stop", stop_list);
  Py_DECREF(stop_list);

  if (opts.has_start)
    PyDict_SetItemString(params, "start", PyUnicode_FromString(opts.start.c_str()));
  else {
    Py_INCREF(Py_None);
    PyDict_SetItemString(params, "start", Py_None);
  }

  return params;
}


static PyObject *
build_trajectory_dict(const kinfold_result_t      &result,
                      const kinfold_trajectory_t  &trajectory,
                      PyObject                    *params)
{
  PyObject *item;
  PyObject *times = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *transcribed = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *structures = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *energies = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *total_energies = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *dna_duplex_energies = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *terminal_hairpin_energies = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *terminal_hairpin_dg_dts = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *terminal_hairpin_plus1_energies = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *rna_dna_minus1_ddgs = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *plus1_ddgs = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *bubble_lengths = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *bubble_lefts = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *hybrid_lengths = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *hybrid_lefts = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *sequence_states = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *moves = PyList_New((Py_ssize_t)trajectory.step_count);
  PyObject *seed_initial = Py_BuildValue("(HHH)",
                                         trajectory.seed_initial[0],
                                         trajectory.seed_initial[1],
                                         trajectory.seed_initial[2]);
  PyObject *seed_final = Py_BuildValue("(HHH)",
                                       trajectory.seed_final[0],
                                       trajectory.seed_final[1],
                                       trajectory.seed_final[2]);
  Py_ssize_t i;

  for (i = 0; i < trajectory.step_count; i++) {
    const kinfold_step_t &step = trajectory.steps[i];

    PyList_SET_ITEM(times, i, PyFloat_FromDouble(step.time));
    PyList_SET_ITEM(transcribed, i, PyLong_FromLong(step.transcribed_length));
    PyList_SET_ITEM(structures, i, PyUnicode_FromString(step.structure ? step.structure : ""));
    PyList_SET_ITEM(energies, i, PyFloat_FromDouble(step.energy));
    PyList_SET_ITEM(total_energies, i, PyFloat_FromDouble(step.total_energy));
    PyList_SET_ITEM(dna_duplex_energies, i, PyFloat_FromDouble(step.dna_duplex_energy));
    PyList_SET_ITEM(terminal_hairpin_energies, i, PyFloat_FromDouble(step.terminal_hairpin_energy));
    PyList_SET_ITEM(terminal_hairpin_dg_dts, i, PyFloat_FromDouble(step.terminal_hairpin_dg_dt));
    PyList_SET_ITEM(terminal_hairpin_plus1_energies, i, PyFloat_FromDouble(step.terminal_hairpin_plus1_energy));
    PyList_SET_ITEM(rna_dna_minus1_ddgs, i, PyFloat_FromDouble(step.rna_dna_minus1_ddg));
    PyList_SET_ITEM(plus1_ddgs, i, PyFloat_FromDouble(step.plus1_ddg));
    PyList_SET_ITEM(bubble_lengths, i, PyLong_FromLong(step.bubble_length));
    PyList_SET_ITEM(bubble_lefts, i, PyLong_FromLong(step.bubble_left));
    PyList_SET_ITEM(hybrid_lengths, i, PyLong_FromLong(step.hybrid_length));
    PyList_SET_ITEM(hybrid_lefts, i, PyLong_FromLong(step.hybrid_left));
    PyList_SET_ITEM(sequence_states, i, PyUnicode_FromString(step.sequence_state ? step.sequence_state : ""));
    if (step.move != '\0')
      PyList_SET_ITEM(moves, i, PyUnicode_FromStringAndSize(&(step.move), 1));
    else {
      Py_INCREF(Py_None);
      PyList_SET_ITEM(moves, i, Py_None);
    }
  }

  item = PyDict_New();
  PyDict_SetItemString(item, "sequence", PyUnicode_FromString(result.sequence ? result.sequence : ""));
  PyDict_SetItemString(item, "parameters", PyDict_Copy(params));
  PyDict_SetItemString(item, "start_structure", PyUnicode_FromString(result.start_structure ? result.start_structure : ""));
  PyDict_SetItemString(item, "time", times);
  PyDict_SetItemString(item, "transcribed_length", transcribed);
  PyDict_SetItemString(item, "structure", structures);
  PyDict_SetItemString(item, "energy", energies);
  PyDict_SetItemString(item, "total_energy", total_energies);
  PyDict_SetItemString(item, "dna_duplex_energy", dna_duplex_energies);
  PyDict_SetItemString(item, "terminal_hairpin_energy", terminal_hairpin_energies);
  PyDict_SetItemString(item, "terminal_hairpin_dg_dt", terminal_hairpin_dg_dts);
  PyDict_SetItemString(item, "terminal_hairpin_plus1_energy", terminal_hairpin_plus1_energies);
  PyDict_SetItemString(item, "rna_dna_minus1_ddg", rna_dna_minus1_ddgs);
  PyDict_SetItemString(item, "plus1_ddg", plus1_ddgs);
  PyDict_SetItemString(item, "bubble_length", bubble_lengths);
  PyDict_SetItemString(item, "bubble_left", bubble_lefts);
  PyDict_SetItemString(item, "hybrid_length", hybrid_lengths);
  PyDict_SetItemString(item, "hybrid_left", hybrid_lefts);
  PyDict_SetItemString(item, "sequence_state", sequence_states);
  PyDict_SetItemString(item, "move", moves);
  PyDict_SetItemString(item, "trajectory_index", PyLong_FromLong(trajectory.trajectory_index));
  PyDict_SetItemString(item, "termination", PyUnicode_FromString(trajectory.termination ? trajectory.termination : ""));
  if (trajectory.stop_index > 0)
    PyDict_SetItemString(item, "stop_index", PyLong_FromLong(trajectory.stop_index));
  else {
    Py_INCREF(Py_None);
    PyDict_SetItemString(item, "stop_index", Py_None);
  }
  PyDict_SetItemString(item, "seed_initial", seed_initial);
  PyDict_SetItemString(item, "seed_final", seed_final);

  Py_DECREF(times);
  Py_DECREF(transcribed);
  Py_DECREF(structures);
  Py_DECREF(energies);
  Py_DECREF(total_energies);
  Py_DECREF(dna_duplex_energies);
  Py_DECREF(terminal_hairpin_energies);
  Py_DECREF(terminal_hairpin_dg_dts);
  Py_DECREF(terminal_hairpin_plus1_energies);
  Py_DECREF(rna_dna_minus1_ddgs);
  Py_DECREF(plus1_ddgs);
  Py_DECREF(bubble_lengths);
  Py_DECREF(bubble_lefts);
  Py_DECREF(hybrid_lengths);
  Py_DECREF(hybrid_lefts);
  Py_DECREF(sequence_states);
  Py_DECREF(moves);
  Py_DECREF(seed_initial);
  Py_DECREF(seed_final);

  return item;
}


static PyObject *
result_to_python(const KinfoldPyOptions &opts,
                 const kinfold_result_t &result)
{
  PyObject *params = build_parameters_dict(opts);
  PyObject *trajectory_list = PyList_New((Py_ssize_t)result.trajectory_count);
  int i;

  for (i = 0; i < result.trajectory_count; i++) {
    PyObject *trajectory = build_trajectory_dict(result, result.trajectories[i], params);
    PyList_SET_ITEM(trajectory_list, i, trajectory);
  }

  Py_DECREF(params);
  return trajectory_list;
}


static bool
parse_sequence_list(PyObject                   *obj,
                    std::vector<std::string>   &out,
                    std::string                &error)
{
  Py_ssize_t i, n;

  out.clear();

  if ((!obj) || (!PySequence_Check(obj)) || PyUnicode_Check(obj)) {
    error = "sequences must be a sequence of strings";
    return false;
  }

  n = PySequence_Size(obj);
  for (i = 0; i < n; i++) {
    PyObject *item = PySequence_GetItem(obj, i);
    std::string value;
    if (!item) {
      error = "Failed to read sequence item";
      return false;
    }
    if (!pyobject_to_string(item, value, error)) {
      Py_DECREF(item);
      return false;
    }
    Py_DECREF(item);
    out.push_back(value);
  }

  return true;
}


static bool
parse_broadcast_or_list_strings(PyObject                   *obj,
                                size_t                     n,
                                std::vector<bool>          &present,
                                std::vector<std::string>   &values,
                                std::string                &error)
{
  size_t i;

  present.assign(n, false);
  values.assign(n, std::string());

  if ((obj == NULL) || (obj == Py_None))
    return true;

  if (PyUnicode_Check(obj)) {
    std::string s;
    if (!pyobject_to_string(obj, s, error))
      return false;
    for (i = 0; i < n; i++) {
      present[i] = true;
      values[i] = s;
    }
    return true;
  }

  if (!PySequence_Check(obj)) {
    error = "Expected a string, None, or a sequence of strings";
    return false;
  }

  if ((size_t)PySequence_Size(obj) != n) {
    error = "Per-sequence string options must match the number of sequences";
    return false;
  }

  for (i = 0; i < n; i++) {
    PyObject *item = PySequence_GetItem(obj, (Py_ssize_t)i);
    if (!item) {
      error = "Failed to read per-sequence string option";
      return false;
    }
    if (item != Py_None) {
      present[i] = true;
      if (!pyobject_to_string(item, values[i], error)) {
        Py_DECREF(item);
        return false;
      }
    }
    Py_DECREF(item);
  }

  return true;
}


static bool
parse_broadcast_or_list_stops(PyObject                              *obj,
                              size_t                                n,
                              std::vector<std::vector<std::string>> &values,
                              std::string                           &error)
{
  size_t i;

  values.assign(n, std::vector<std::string>());

  if ((obj == NULL) || (obj == Py_None))
    return true;

  if (PyUnicode_Check(obj)) {
    std::string stop;
    if (!pyobject_to_string(obj, stop, error))
      return false;
    for (i = 0; i < n; i++)
      values[i].push_back(stop);
    return true;
  }

  if (!PySequence_Check(obj)) {
    error = "stops must be None, a string, or a per-sequence sequence";
    return false;
  }

  if ((size_t)PySequence_Size(obj) != n) {
    error = "Per-sequence stops must match the number of sequences";
    return false;
  }

  for (i = 0; i < n; i++) {
    PyObject *item = PySequence_GetItem(obj, (Py_ssize_t)i);
    if (!item) {
      error = "Failed to read stop collection";
      return false;
    }
    if (!parse_stop_list(item, values[i], error)) {
      Py_DECREF(item);
      return false;
    }
    Py_DECREF(item);
  }

  return true;
}


static bool
parse_broadcast_or_list_seeds(PyObject                 *obj,
                              size_t                   n,
                              std::vector<bool>        &present,
                              std::vector<std::string> &values,
                              std::string              &error)
{
  size_t i;

  present.assign(n, false);
  values.assign(n, std::string());

  if ((obj == NULL) || (obj == Py_None))
    return true;

  if (PyUnicode_Check(obj) || PyTuple_Check(obj) || PyList_Check(obj)) {
    if ((PySequence_Check(obj)) && (!PyUnicode_Check(obj)) && ((size_t)PySequence_Size(obj) == n)) {
      for (i = 0; i < n; i++) {
        PyObject *item = PySequence_GetItem(obj, (Py_ssize_t)i);
        bool     present_value = false;
        if (!item) {
          error = "Failed to read per-sequence seed";
          return false;
        }
        if (!parse_seed_string(item, present_value, values[i], error)) {
          Py_DECREF(item);
          return false;
        }
        present[i] = present_value;
        Py_DECREF(item);
      }
      return true;
    }

    {
      bool present_one;
      std::string one;
      if (!parse_seed_string(obj, present_one, one, error))
        return false;
      for (i = 0; i < n; i++) {
        present[i] = present_one;
        values[i] = one;
      }
    }
    return true;
  }

  error = "seed/seeds must be None, a string/3-tuple, or a per-sequence sequence";
  return false;
}

} /* namespace */


PyObject *
kinfold_python_run(std::string sequence,
                   PyObject    *start_obj,
                   PyObject    *stop_obj,
                   int         dangle,
                   double      Temp,
                   std::string Par,
                   bool        logML,
                   bool        noShift,
                   bool        noLP,
                   PyObject    *seed_obj,
                   double      time,
                   double      time_step,
                   int         num,
                   bool        met,
                   bool        fpt,
                   bool        rect,
                   int         glen,
                   double      transcription_elongation_rate,
                   double      kinfold_seconds_per_time_unit,
                   int         max_bubble_width,
                   PyObject    *phi_obj,
                   PyObject    *pbounds_obj,
                   std::string log,
                   bool        silent,
                   bool        verbose,
                   bool        lmin,
                   double      cut)
{
  KinfoldPyOptions opts;
  kinfold_result_t result;
  PyObject         *items;
  std::string      error;

  opts.sequence = sequence;
  opts.dangle = dangle;
  opts.temp = Temp;
  opts.has_par = !Par.empty();
  opts.par = Par;
  opts.logML = logML;
  opts.noShift = noShift;
  opts.noLP = noLP;
  opts.time = time;
  opts.time_step = time_step;
  opts.num = num;
  opts.met = met;
  opts.fpt = fpt;
  opts.rect = rect;
  opts.has_glen = (glen >= 0);
  opts.glen = glen;
  opts.transcription_elongation_rate = transcription_elongation_rate;
  opts.kinfold_seconds_per_time_unit = kinfold_seconds_per_time_unit;
  opts.max_bubble_width = max_bubble_width;
  opts.log = log;
  opts.silent = silent;
  opts.verbose = verbose;
  opts.lmin = lmin;
  opts.cut = (cut >= 1e299) ? std::numeric_limits<double>::infinity() : cut;

  if (!parse_optional_string(start_obj, opts.has_start, opts.start, error) ||
      !parse_stop_list(stop_obj, opts.stops, error) ||
      !parse_seed_string(seed_obj, opts.has_seed, opts.seed, error) ||
      !parse_optional_double(phi_obj, opts.has_phi, opts.phi, error) ||
      !parse_pbounds_string(pbounds_obj, opts.has_pbounds, opts.pbounds, error)) {
    PyErr_SetString(PyExc_ValueError, error.c_str());
    return NULL;
  }

  if (!run_kinfold_capture(opts, result, error)) {
    PyErr_SetString(PyExc_ValueError, error.c_str());
    return NULL;
  }

  items = result_to_python(opts, result);
  kinfold_result_free(&result);

  if (num == 1) {
    PyObject *first = PySequence_GetItem(items, 0);
    Py_DECREF(items);
    return first;
  }

  return items;
}


PyObject *
kinfold_python_batch(PyObject    *sequences_obj,
                     PyObject    *starts_obj,
                     PyObject    *stops_obj,
                     int         trajectories_per_sequence,
                     int         dangle,
                     double      Temp,
                     std::string Par,
                     bool        logML,
                     bool        noShift,
                     bool        noLP,
                     PyObject    *seeds_obj,
                     double      time,
                     double      time_step,
                     bool        met,
                     bool        fpt,
                     bool        rect,
                     int         glen,
                     double      transcription_elongation_rate,
                     double      kinfold_seconds_per_time_unit,
                     int         max_bubble_width,
                     PyObject    *phi_obj,
                     PyObject    *pbounds_obj,
                     std::string log,
                     bool        silent,
                     bool        verbose,
                     bool        lmin,
                     double      cut)
{
  std::vector<std::string>              sequences;
  std::vector<bool>                     has_starts;
  std::vector<std::string>              starts;
  std::vector<std::vector<std::string>> stops;
  std::vector<bool>                     has_seeds;
  std::vector<std::string>              seeds;
  KinfoldPyOptions                      opts;
  PyObject                              *grouped;
  std::string                           error;
  size_t                                i;

  if (!parse_sequence_list(sequences_obj, sequences, error) ||
      !parse_broadcast_or_list_strings(starts_obj, sequences.size(), has_starts, starts, error) ||
      !parse_broadcast_or_list_stops(stops_obj, sequences.size(), stops, error) ||
      !parse_broadcast_or_list_seeds(seeds_obj, sequences.size(), has_seeds, seeds, error) ||
      !parse_optional_double(phi_obj, opts.has_phi, opts.phi, error) ||
      !parse_pbounds_string(pbounds_obj, opts.has_pbounds, opts.pbounds, error)) {
    PyErr_SetString(PyExc_ValueError, error.c_str());
    return NULL;
  }

  grouped = PyList_New((Py_ssize_t)sequences.size());

  for (i = 0; i < sequences.size(); i++) {
    kinfold_result_t  result;
    PyObject          *group;
    PyObject          *traj_list;

    opts.sequence = sequences[i];
    opts.has_start = has_starts[i];
    opts.start = starts[i];
    opts.stops = stops[i];
    opts.dangle = dangle;
    opts.temp = Temp;
    opts.has_par = !Par.empty();
    opts.par = Par;
    opts.logML = logML;
    opts.noShift = noShift;
    opts.noLP = noLP;
    opts.has_seed = has_seeds[i];
    opts.seed = seeds[i];
    opts.time = time;
    opts.time_step = time_step;
    opts.num = trajectories_per_sequence;
    opts.met = met;
    opts.fpt = fpt;
    opts.rect = rect;
    opts.has_glen = (glen >= 0);
    opts.glen = glen;
    opts.transcription_elongation_rate = transcription_elongation_rate;
    opts.kinfold_seconds_per_time_unit = kinfold_seconds_per_time_unit;
    opts.max_bubble_width = max_bubble_width;
    opts.log = log;
    opts.silent = silent;
    opts.verbose = verbose;
    opts.lmin = lmin;
    opts.cut = (cut >= 1e299) ? std::numeric_limits<double>::infinity() : cut;

    if (!run_kinfold_capture(opts, result, error)) {
      Py_DECREF(grouped);
      PyErr_SetString(PyExc_ValueError, error.c_str());
      return NULL;
    }

    traj_list = result_to_python(opts, result);
    group = PyDict_New();
    PyDict_SetItemString(group, "sequence", PyUnicode_FromString(result.sequence ? result.sequence : ""));
    PyDict_SetItemString(group, "parameters", build_parameters_dict(opts));
    PyDict_SetItemString(group, "start_structure", PyUnicode_FromString(result.start_structure ? result.start_structure : ""));
    PyDict_SetItemString(group, "trajectories", traj_list);
    PyList_SET_ITEM(grouped, (Py_ssize_t)i, group);
    Py_DECREF(traj_list);
    kinfold_result_free(&result);
  }

  return grouped;
}
