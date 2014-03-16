#include <stdio.h>
#include "io.h"

#ifndef IO_OUTBIN

#define VERSION_TAG "SPHView00 "


void write_header(FILE* fp, int n, int framecount, float h)
{
    fprintf(fp, "%s%d %d %g\n", VERSION_TAG, n, framecount, h);
}


void write_frame_data(FILE* fp, int n, sim_state_t* s, int* c)
{
    particle_t* p = s->part;
    for (int i = 0; i < n; ++i, ++p)
        fprintf(fp, "%e %e %e\n", p->x[0], p->x[1], p->x[2]);
}

#endif /* IO_OUTBIN */
