#include <glib.h>
#include "sbm_dirs.h"

void initialize (void);
void run (double);
void finalize (void);

gchar*
read_data_file (gchar* file)
{
  gchar* contents = NULL;
  GError* error = NULL;
  gsize len = 0;

  g_assert (file!=NULL);

  g_file_get_contents (file, &contents, &len, &error);

  g_assert (error==NULL);
  g_assert (len>0);
  g_assert (contents!=NULL);

  return contents;
}

gchar*
read_benchmark_file (gchar* file)
{
  gchar* contents = NULL;
  char* path = g_build_filename (SBM_DATADIR, "benchmark_0", file, NULL);

  g_assert (path);
  contents = read_data_file (path);
  g_assert (contents!=NULL);

  return contents;
}

void
test_murray_1 (void)
{
  gchar* contents_bench = NULL;
  gchar* contents_test = NULL;

  initialize ();
  run (1);
  finalize ();

  contents_bench = read_benchmark_file ("ripp.0001");
  g_assert (contents_bench);

  contents_test = read_data_file ("ripp.0001");
  g_assert (contents_test);

  g_assert_cmpstr (contents_bench, ==, contents_test);

  contents_bench = read_benchmark_file ("save.0001");
  g_assert (contents_bench);

  contents_test = read_data_file ("save.0001");
  g_assert (contents_test);

  g_assert_cmpstr (contents_bench, ==, contents_test);
}

void
test_murray_0 (void)
{
  gchar* contents_bench = NULL;
  gchar* contents_test = NULL;

  initialize ();
  run (.002);
  finalize ();

  contents_bench = read_benchmark_file ("ripp.0001");
  g_assert (contents_bench);

  contents_test = read_data_file ("ripp.0001");
  g_assert (contents_test);

  g_assert_cmpstr (contents_bench, ==, contents_test);

  contents_bench = read_benchmark_file ("save.0001");
  g_assert (contents_bench);

  contents_test = read_data_file ("save.0001");
  g_assert (contents_test);

  g_assert_cmpstr (contents_bench, ==, contents_test);
}

int
main (int argc, char* argv[])
{
  g_test_init (&argc, &argv, NULL);

  if (g_test_slow ())
  {
    g_test_add_func ("/murray/main/1", &test_murray_1);
  }
  g_test_add_func ("/murray/main/0", &test_murray_0);

  g_test_run ();
  return 0;
}

