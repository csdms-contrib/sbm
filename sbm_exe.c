
void initialize (void);
void run (double);
void finalize (void);

int
main (void)
{
  initialize ();
  run (-1);
  finalize ();

  return 0;
}

