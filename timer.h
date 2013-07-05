#include <sys/time.h>

void _tstart(void);
void _tend(void);
void _tdiff(struct timeval *tp);
void _tadd(struct timeval tp, struct timeval *tpsum);


