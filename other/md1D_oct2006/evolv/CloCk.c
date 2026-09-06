#include <time.h>

#define CHRONO(t2,t1) ((double)(t2-t1)/CLOCKS_PER_SEC)

int Hey_Hey_attend(double t_en_sec)
{
	clock_t clock0
	while(CHRONO(clock(), clock0) < t_en_sec){}
	return 0;
}
