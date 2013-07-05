#include <stdlib.h>
#include "timer.h"


/* For experimental measurements!!! */
static struct timeval tp1, tp2;

void _tstart(void)
{
  gettimeofday(&tp1, NULL);
}

void _tend(void)
{
  gettimeofday(&tp2, NULL);
}

double tdiff(struct timeval tpb, struct timeval tpe)
{
  struct timeval tp;
  double t;
  tp.tv_sec = tpe.tv_sec - tpb.tv_sec;
  tp.tv_usec = tpe.tv_usec - tpb.tv_usec;
  if (tp.tv_usec < 0) {
    tp.tv_sec--;
    tp.tv_usec += 1000000;
  }
  t = tp.tv_sec + ((double)tp.tv_usec)/1000000.0;
  return t;
}

void _tdiff(struct timeval *tp)
{
  tp->tv_sec = tp2.tv_sec - tp1.tv_sec;
  tp->tv_usec = tp2.tv_usec - tp1.tv_usec;
  if (tp->tv_usec < 0) {
    tp->tv_sec--;
    tp->tv_usec += 1000000;
  }
}

void _tadd(struct timeval tp, struct timeval *tpsum)
{
  tpsum->tv_sec  += tp.tv_sec;
  tpsum->tv_usec += tp.tv_usec;
  if (tpsum->tv_usec > 999999) {
    tpsum->tv_sec++;
    tpsum->tv_usec -= 1000000;
  }
}
