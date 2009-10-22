#include <glib.h>

void
test_murray_0 (void)
{
}

int
main (int argc, char* argv[])
{
  g_test_init (&argc, &argv, NULL);

  if (g_test_slow ())
  {
    g_test_add_func ("murray/main/0", &test_murray_0);
  }

  g_test_run ();
  return 0;
}

