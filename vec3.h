#ifndef VEC3_H
#define VEC3_H

/*@T
 * \subsection{Vector operations}
 * 
 * A few inline functions save us from typing the same annoying bits
 * of code repeatedly.
 *
 *@c*/

inline void vec3_set(float* result, float x, float y, float z)
{
    result[0] = x;
    result[1] = y;
    result[2] = z;
}

inline void vec3_copy(float* result, float* v)
{
    vec3_set(result, v[0], v[1], v[2]);
}

inline void vec3_diff(float* result, float* a, float* b)
{
    result[0] = a[0]-b[0];
    result[1] = a[1]-b[1];
    result[2] = a[2]-b[2];
}

inline void vec3_scale(float* result, float alpha, float* v)
{
    result[0] = alpha*v[0];
    result[1] = alpha*v[1];
    result[2] = alpha*v[2];
}

inline float vec3_dist2(float* a, float* b)
{
    float dx = a[0]-b[0];
    float dy = a[1]-b[1];
    float dz = a[2]-b[2];
    return dx*dx + dy*dy + dz*dz;
}

inline float vec3_len2(float* a)
{
    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
}

inline void vec3_saxpy(float* result, float alpha, float* v)
{
    result[0] += alpha*v[0];
    result[1] += alpha*v[1];
    result[2] += alpha*v[2];
}

inline void vec3_scalev(float* result, float alpha)
{
    result[0] *= alpha;
    result[1] *= alpha;
    result[2] *= alpha;
}

/*@q*/
#endif /* VEC3_H */
