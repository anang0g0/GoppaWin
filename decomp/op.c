#include <stdio.h>
#include <stdlib.h>

//#include "global.h"
//#include "struct.h"
#include "struct.h"
#include "oplib.h"


//多項式を表示する(default)
void printpol (vec a){
  int i, n;

  n = deg (a);
  if (n < 0)
    {
      printf ("baka\n");
      exit (1);
    }


  for (i = n; i > -1; i--)
    {
      if (a.x[i] > 0)
	{
	  printf ("%u", a.x[i]);
	  if (i > 0)
	    printf ("x^%d", i);
	  if (i > 0)
	    printf ("+");
	}
    }
  //  printf("\n");

  return;
}
