#ifndef STATE_H
#define STATE_H

#define STATE_HASH_SIZE 4096

/*@T
 * \section{System state}
 * 
 * The [[sim_state_t]] structure holds the information for the current
 * state of the system and of the integration algorithm.
 * 
 * The [[alloc_state]] and [[free_state]] functions take care of storage
 * for the local simulation state.
 *@c*/
typedef struct particle_t {
    float rho;         /* Particle density   */
    float x[3];        /* Particle positions */
    float v[3];        /* Particle velocities (full step) */
    float vh[3];       /* Particle velocities (half step) */
    float a[3];        /* Particle accelerations */
    struct particle_t* next;  /* List link for spatial hashing */
} particle_t;

typedef struct sim_state_t {
    int n;                /* Number of particles    */
    float mass;           /* Particle mass          */
    particle_t* part;     /* Particles              */
    particle_t** hash;    /* Hash table             */
} sim_state_t;

sim_state_t* alloc_state(int n);
void free_state(sim_state_t* s);

/*@q*/
#endif /* STATE_H */
