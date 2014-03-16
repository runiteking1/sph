#ifndef IO_H
#define IO_H

#include "state.h"

void write_header(FILE* fp, int n, int framecount, float h);
void write_frame_data(FILE* fp, int n, sim_state_t* state, int* c);

#endif /* IO_H */
