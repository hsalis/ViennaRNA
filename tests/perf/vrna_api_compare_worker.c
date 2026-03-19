#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/subopt/wuchty.h>
#include <ViennaRNA/utils/cpu.h>


typedef struct {
  char  *id;
  char  *sequence;
  int   delta;
} workload_item_t;


static long long
now_ns(void)
{
  struct timespec ts;

  clock_gettime(CLOCK_MONOTONIC, &ts);

  return (long long)ts.tv_sec * 1000000000LL + (long long)ts.tv_nsec;
}


static void
write_json_string(FILE       *fp,
                  const char *s)
{
  const unsigned char *p = (const unsigned char *)s;

  fputc('"', fp);

  while (*p) {
    switch (*p) {
      case '\\':
      case '"':
        fputc('\\', fp);
        fputc(*p, fp);
        break;

      case '\n':
        fputs("\\n", fp);
        break;

      case '\r':
        fputs("\\r", fp);
        break;

      case '\t':
        fputs("\\t", fp);
        break;

      default:
        fputc(*p, fp);
        break;
    }

    p++;
  }

  fputc('"', fp);
}


static workload_item_t *
read_workload(const char  *path,
              const char  *task,
              size_t      *count)
{
  FILE            *fp;
  char            line[16384];
  size_t          used = 0;
  size_t          cap = 32;
  workload_item_t *items;

  items = (workload_item_t *)calloc(cap, sizeof(*items));
  if (!items)
    return NULL;

  fp = fopen(path, "r");
  if (!fp) {
    free(items);
    return NULL;
  }

  while (fgets(line, sizeof(line), fp)) {
    char *id;
    char *sequence;
    char *delta_s;
    char *nl;

    nl = strchr(line, '\n');
    if (nl)
      *nl = '\0';

    if (line[0] == '\0')
      continue;

    if (used == cap) {
      workload_item_t *tmp;

      cap *= 2;
      tmp = (workload_item_t *)realloc(items, cap * sizeof(*items));
      if (!tmp) {
        fclose(fp);
        free(items);
        return NULL;
      }
      items = tmp;
    }

    id = strtok(line, "\t");
    if (!id)
      continue;

    if (strcmp(task, "subopt") == 0) {
      delta_s   = strtok(NULL, "\t");
      sequence  = strtok(NULL, "\t");

      if ((!delta_s) || (!sequence))
        continue;

      items[used].delta = atoi(delta_s);
    } else {
      sequence = strtok(NULL, "\t");
      if (!sequence)
        continue;
      items[used].delta = 0;
    }

    items[used].id        = strdup(id);
    items[used].sequence  = strdup(sequence);
    used++;
  }

  fclose(fp);
  *count = used;

  return items;
}


static void
free_workload(workload_item_t *items,
              size_t          count)
{
  size_t i;

  if (!items)
    return;

  for (i = 0; i < count; i++) {
    free(items[i].id);
    free(items[i].sequence);
  }

  free(items);
}


static void
write_feature_names(FILE         *fp,
                    unsigned int caps)
{
  int first = 1;

  fputc('[', fp);

  if (caps & VRNA_CPU_SIMD_SSE2) {
    write_json_string(fp, "SSE2");
    first = 0;
  }

  if (caps & VRNA_CPU_SIMD_SSE3) {
    if (!first)
      fputc(',', fp);
    write_json_string(fp, "SSE3");
    first = 0;
  }

  if (caps & VRNA_CPU_SIMD_SSE41) {
    if (!first)
      fputc(',', fp);
    write_json_string(fp, "SSE4.1");
    first = 0;
  }

  if (caps & VRNA_CPU_SIMD_SSE42) {
    if (!first)
      fputc(',', fp);
    write_json_string(fp, "SSE4.2");
    first = 0;
  }

  if (caps & VRNA_CPU_SIMD_AVX) {
    if (!first)
      fputc(',', fp);
    write_json_string(fp, "AVX");
    first = 0;
  }

  if (caps & VRNA_CPU_SIMD_AVX2) {
    if (!first)
      fputc(',', fp);
    write_json_string(fp, "AVX2");
    first = 0;
  }

  if (caps & VRNA_CPU_SIMD_AVX512F) {
    if (!first)
      fputc(',', fp);
    write_json_string(fp, "AVX512F");
    first = 0;
  }

  if (caps & VRNA_CPU_SIMD_AVX512BW) {
    if (!first)
      fputc(',', fp);
    write_json_string(fp, "AVX512BW");
    first = 0;
  }

  if (caps & VRNA_CPU_SIMD_AVX512VL) {
    if (!first)
      fputc(',', fp);
    write_json_string(fp, "AVX512VL");
    first = 0;
  }

  if (caps & VRNA_CPU_SIMD_NEON) {
    if (!first)
      fputc(',', fp);
    write_json_string(fp, "NEON");
  }

  fputc(']', fp);
}


static void
write_compiled_slices(FILE *fp)
{
  int first = 1;

  fputc('[', fp);
  write_json_string(fp, "scalar");
  first = 0;

#if VRNA_WITH_SIMD_SSE41
  if (!first)
    fputc(',', fp);
  write_json_string(fp, "SSE4.1");
  first = 0;
#endif

#if VRNA_WITH_SIMD_AVX2
  if (!first)
    fputc(',', fp);
  write_json_string(fp, "AVX2");
  first = 0;
#endif

#if VRNA_WITH_SIMD_AVX512
  if (!first)
    fputc(',', fp);
  write_json_string(fp, "AVX512F");
  first = 0;
#endif

#if VRNA_WITH_SIMD_NEON
  if (!first)
    fputc(',', fp);
  write_json_string(fp, "NEON");
#endif

  fputc(']', fp);
}


static int
write_probe(const char *output_path)
{
  FILE         *fp;
  unsigned int caps;

  fp = fopen(output_path, "w");
  if (!fp)
    return 1;

  caps = vrna_cpu_simd_capabilities();

  fputc('{', fp);
  fputs("\"cpu_vendor\":", fp);
  write_json_string(fp, vrna_cpu_vendor_string());
  fputs(",\"runtime_capabilities\":", fp);
  fprintf(fp, "%u", caps);
  fputs(",\"runtime_features\":", fp);
  write_feature_names(fp, caps);
  fputs(",\"compiled_slices\":", fp);
  write_compiled_slices(fp);
  fputs("}\n", fp);

  fclose(fp);
  return 0;
}


static int
run_fold_like(FILE                  *fp,
              const char            *task,
              const workload_item_t *item,
              int                   repeats,
              int                   warmups)
{
  int                 r;
  vrna_md_t           md;
  float               energy = 0.;
  char                *structure = NULL;
  vrna_fold_compound_t *fc;

  vrna_md_set_default(&md);

  fprintf(fp, "{\"id\":");
  write_json_string(fp, item->id);
  fputs(",\"meta\":{\"sequence\":", fp);
  write_json_string(fp, item->sequence);
  if (strcmp(task, "cofold") == 0) {
    const char *cut = strchr(item->sequence, '&');
    int        len1 = cut ? (int)(cut - item->sequence) : (int)strlen(item->sequence);
    int        len2 = cut ? (int)strlen(cut + 1) : 0;
    fprintf(fp, ",\"len1\":%d,\"len2\":%d", len1, len2);
  } else {
    fprintf(fp, ",\"length\":%d", (int)strlen(item->sequence));
  }
  fputs("},\"durations_ns\":[", fp);

  for (r = 0; r < warmups; r++) {
    fc = vrna_fold_compound(item->sequence, &md, VRNA_OPTION_MFE);
    if (!fc)
      return 1;

    structure = (char *)malloc(strlen(item->sequence) + 1);
    if (!structure)
      return 1;

    (void)vrna_mfe(fc, structure);
    free(structure);
    vrna_fold_compound_free(fc);
  }

  for (r = 0; r < repeats; r++) {
    long long start;
    long long elapsed;

    fc = vrna_fold_compound(item->sequence, &md, VRNA_OPTION_MFE);
    if (!fc)
      return 1;

    structure = (char *)malloc(strlen(item->sequence) + 1);
    if (!structure)
      return 1;

    start   = now_ns();
    energy   = vrna_mfe(fc, structure);
    elapsed  = now_ns() - start;

    if (r > 0)
      fputc(',', fp);
    fprintf(fp, "%lld", elapsed);

    if (r + 1 == repeats) {
      fputs("],\"output\":{\"structure\":", fp);
      write_json_string(fp, structure);
      fputs(",\"energy\":", fp);
      fprintf(fp, "%.8f", energy);
      fputs("}}", fp);
    }

    free(structure);
    vrna_fold_compound_free(fc);
  }

  return 0;
}


static void
write_subopt_solutions(FILE                  *fp,
                       vrna_subopt_solution_t *sol)
{
  int i;

  fputc('[', fp);

  for (i = 0; sol[i].structure; i++) {
    if (i > 0)
      fputc(',', fp);

    fputs("{\"structure\":", fp);
    write_json_string(fp, sol[i].structure);
    fputs(",\"energy\":", fp);
    fprintf(fp, "%.8f}", sol[i].energy);
  }

  fputc(']', fp);
}


static void
free_subopt_solutions(vrna_subopt_solution_t *sol)
{
  int i;

  if (!sol)
    return;

  for (i = 0; sol[i].structure; i++)
    free(sol[i].structure);

  free(sol);
}


static int
run_subopt(FILE                  *fp,
           const workload_item_t *item,
           int                   repeats,
           int                   warmups)
{
  int                   r;
  vrna_md_t             md;
  vrna_fold_compound_t  *fc;
  vrna_subopt_solution_t *sol = NULL;

  vrna_md_set_default(&md);
  md.uniq_ML = 1;

  fprintf(fp, "{\"id\":");
  write_json_string(fp, item->id);
  fprintf(fp, ",\"meta\":{\"length\":%d,\"delta\":%d,\"sequence\":",
          (int)strlen(item->sequence),
          item->delta);
  write_json_string(fp, item->sequence);
  fputs("},\"durations_ns\":[", fp);

  for (r = 0; r < warmups; r++) {
    fc = vrna_fold_compound(item->sequence, &md, VRNA_OPTION_MFE);
    if (!fc)
      return 1;
    sol = vrna_subopt(fc, item->delta, VRNA_SORT_BY_ENERGY_LEXICOGRAPHIC_ASC, NULL);
    free_subopt_solutions(sol);
    vrna_fold_compound_free(fc);
  }

  for (r = 0; r < repeats; r++) {
    long long start;
    long long elapsed;

    fc = vrna_fold_compound(item->sequence, &md, VRNA_OPTION_MFE);
    if (!fc)
      return 1;

    start   = now_ns();
    sol     = vrna_subopt(fc, item->delta, VRNA_SORT_BY_ENERGY_LEXICOGRAPHIC_ASC, NULL);
    elapsed = now_ns() - start;

    if (r > 0)
      fputc(',', fp);
    fprintf(fp, "%lld", elapsed);

    if (r + 1 == repeats) {
      fputs("],\"output\":{\"solutions\":", fp);
      write_subopt_solutions(fp, sol);
      fputs("}}", fp);
    }

    free_subopt_solutions(sol);
    vrna_fold_compound_free(fc);
  }

  return 0;
}


int
main(int argc, char *argv[])
{
  const char      *task;
  const char      *input_path;
  const char      *output_path;
  int             repeats;
  int             warmups;
  size_t          count;
  size_t          i;
  FILE            *fp;
  workload_item_t *items;

  if (argc < 2)
    return 1;

  task = argv[1];

  if (strcmp(task, "probe") == 0) {
    if (argc != 3)
      return 1;
    return write_probe(argv[2]);
  }

  if (argc != 6)
    return 1;

  input_path   = argv[2];
  output_path  = argv[3];
  repeats      = atoi(argv[4]);
  warmups      = atoi(argv[5]);

  items = read_workload(input_path, task, &count);
  if (!items)
    return 1;

  fp = fopen(output_path, "w");
  if (!fp) {
    free_workload(items, count);
    return 1;
  }

  fputs("{\"task\":", fp);
  write_json_string(fp, task);
  fputs(",\"results\":[", fp);

  for (i = 0; i < count; i++) {
    int ret;

    if (i > 0)
      fputc(',', fp);

    if (strcmp(task, "subopt") == 0)
      ret = run_subopt(fp, &items[i], repeats, warmups);
    else
      ret = run_fold_like(fp, task, &items[i], repeats, warmups);

    if (ret) {
      fclose(fp);
      free_workload(items, count);
      return 1;
    }
  }

  fputs("]}\n", fp);
  fclose(fp);
  free_workload(items, count);

  return 0;
}
