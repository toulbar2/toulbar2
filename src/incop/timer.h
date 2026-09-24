#ifndef _timer_h_
#define _timer_h_

typedef double Time;
typedef enum type_timer { REAL,
    VIRTUAL } Timer;

#ifndef __WIN32__
#include "tb2config.hpp"
void start_timers(void);
void stop_timers(Timer type);
extern TB2_THREAD_LOCAL Time real_lapse;
extern TB2_THREAD_LOCAL Time virtual_ulapse, virtual_slapse;
#define REAL_TIMELAPSE (real_lapse)
#define VIRTUAL_TIMELAPSE (virtual_ulapse + virtual_slapse)
#endif

#endif
