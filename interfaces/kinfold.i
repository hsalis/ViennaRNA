/**********************************************/
/* BEGIN Python interface for Kinfold         */
/**********************************************/

#ifdef SWIGPYTHON
%{
PyObject *
kinfold_python_run(std::string sequence,
                   PyObject    *start,
                   PyObject    *stop,
                   int         dangle,
                   double      Temp,
                   std::string Par,
                   bool        logML,
                   bool        noShift,
                   bool        noLP,
                   PyObject    *seed,
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
                   PyObject    *phi,
                   PyObject    *pbounds,
                   std::string log,
                   bool        silent,
                   bool        verbose,
                   bool        lmin,
                   double      cut);

PyObject *
kinfold_python_batch(PyObject    *sequences,
                     PyObject    *starts,
                     PyObject    *stops,
                     int         trajectories_per_sequence,
                     int         dangle,
                     double      Temp,
                     std::string Par,
                     bool        logML,
                     bool        noShift,
                     bool        noLP,
                     PyObject    *seeds,
                     double      time,
                     double      time_step,
                     bool        met,
                     bool        fpt,
                     bool        rect,
                     int         glen,
                     double      transcription_elongation_rate,
                     double      kinfold_seconds_per_time_unit,
                     int         max_bubble_width,
                     PyObject    *phi,
                     PyObject    *pbounds,
                     std::string log,
                     bool        silent,
                     bool        verbose,
                     bool        lmin,
                     double      cut);
%}

%rename(_kinfold_run_raw) kinfold_python_run;
%rename(_kinfold_batch_raw) kinfold_python_batch;

%feature("autodoc") kinfold_python_run;
%feature("kwargs") kinfold_python_run;
%feature("autodoc") kinfold_python_batch;
%feature("kwargs") kinfold_python_batch;

PyObject *
kinfold_python_run(std::string sequence,
                   PyObject    *start = Py_None,
                   PyObject    *stop = Py_None,
                   int         dangle = 2,
                   double      Temp = 37.0,
                   std::string Par = "",
                   bool        logML = true,
                   bool        noShift = false,
                   bool        noLP = false,
                   PyObject    *seed = Py_None,
                   double      time = 500.0,
                   double      time_step = 1.0,
                   int         num = 1,
                   bool        met = false,
                   bool        fpt = true,
                   bool        rect = false,
                   int         glen = -1,
                   double      transcription_elongation_rate = 0.0,
                   double      kinfold_seconds_per_time_unit = 1e-5,
                   int         max_bubble_width = 0,
                   PyObject    *phi = Py_None,
                   PyObject    *pbounds = Py_None,
                   std::string log = "kinout",
                   bool        silent = false,
                   bool        verbose = false,
                   bool        lmin = false,
                   double      cut = 1e300);

PyObject *
kinfold_python_batch(PyObject    *sequences,
                     PyObject    *starts = Py_None,
                     PyObject    *stops = Py_None,
                     int         trajectories_per_sequence = 1,
                     int         dangle = 2,
                     double      Temp = 37.0,
                     std::string Par = "",
                     bool        logML = true,
                     bool        noShift = false,
                     bool        noLP = false,
                     PyObject    *seeds = Py_None,
                     double      time = 500.0,
                     double      time_step = 1.0,
                     bool        met = false,
                     bool        fpt = true,
                     bool        rect = false,
                     int         glen = -1,
                     double      transcription_elongation_rate = 0.0,
                     double      kinfold_seconds_per_time_unit = 1e-5,
                     int         max_bubble_width = 0,
                     PyObject    *phi = Py_None,
                     PyObject    *pbounds = Py_None,
                     std::string log = "kinout",
                     bool        silent = false,
                     bool        verbose = false,
                     bool        lmin = false,
                     double      cut = 1e300);

%pythoncode %{
def kinfold(sequence, **kwargs):
    return _kinfold_run_raw(sequence, **kwargs)


def kinfold_batch(sequences, **kwargs):
    return _kinfold_batch_raw(sequences, **kwargs)
%}
#endif
