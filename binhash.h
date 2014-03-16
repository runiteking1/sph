#ifndef BINHASH_H
#define BINHASH_H

#include "state.h"

/*@T
 * \section{Spatial hashing}
 *
 * We conceptually partition our computational domain into bins that
 * are at least $h$ on a side, and label each with integer coordinates
 * $(i_x, i_y, i_z)$.  In three dimensions, the bin size doesn't have
 * to be all that small before the number of such bins is quite large,
 * and most bins will be empty.  Rather than represent each bin
 * explicitly, we will map the bins to locations into a hash table,
 * allowing the possibility that different bins can map to the same
 * location in the hash table (though ideally this should not happen
 * too often).  We compute the hash function by mapping $(i_x, i_y, i_z)$ 
 * to a Z-Morton integer index, which we then associate with a hash
 * bucket.  By figuring out the hash buckets in which the neighbor of
 * a given particle could possibly lie, we significantly reduce the
 * cost of checking interactions.
 *
 * In the current implementation, we set the bin size equal to $h$,
 * which implies that particles interacting with a given particle
 * might lie in any of 27 possible neighbors.  If you use a bin of
 * size $2h$, you may have more particles in each bin, but you would
 * only need to check eight bins (at most) for possible interactions.
 * In order to allow for the possibility of more or fewer possible
 * bins containing neighbors, we define [[MAX_NBR_BINS]] to be the
 * maximum number of bins we will ever need to check for interactions,
 * and let [[particle_neighborhood]] return both which bins are
 * involved (as an output argument) and the number of bins needed (via
 * the return value).
 *
 * We also currently use the last few bits of each of the $i_x, i_y,
 * i_z$ indices to form the Z-Morton index.  You may want to change
 * the number of bits used, or change to some other hashing scheme
 * that potentially spreads the bins more uniformly across the table.
 *
 *@c*/

#define HASH_DIM 0x10
#define HASH_SIZE (HASH_DIM*HASH_DIM*HASH_DIM)
#define MAX_NBR_BINS 27

unsigned particle_bucket(particle_t* p, float h);
unsigned particle_neighborhood(unsigned* buckets, particle_t* p, float h);
void hash_particles(sim_state_t* s, float h);

/*@q*/
#endif /* BINHASH_H */
