#ifndef ZMORTON_H
#define ZMORTON_H

/*@T
 * \subsection{Z-Morton encoding}
 *
 * We use a Z-Morton encoding to map a triple of integer indices
 * $(i_x, i_y, i_z)$ to a single integer.  If you want to see a
 * picture of the Z-Morton ordering, I recommend the Wikipedia page!
 * In terms of computation, though, the Z-Morton ordering simply
 * interleaves the bits of the three independent indices.  That is,
 * if $i_x^j$ is the bit in the $2^j$ place for index $i_x$, the
 * bit pattern for the Z-Morton code looks like
 * \[
 *   c = (\ldots i_z^1 ~ i_y^1 ~ i_x^1 ~ i_z^0 ~ i_y^0 ~ i_x^0)_2.
 * \]
 * While this is not as good a space-filling curve as a Hilbert
 * curve, it's very cheap.
 *
 * The concrete code combines three 10-bit (max) indices into a single
 * 32-bit Z-Morton code.  The code is adapted from
 * \begin{center}
 *   \url{http://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/}
 * \end{center}
 *@c*/

inline unsigned zm_part1by2(unsigned x)
{
    x &= 0x000003ff;
    x = (x ^ (x << 16)) & 0xff0000ff;
    x = (x ^ (x <<  8)) & 0x0300f00f;
    x = (x ^ (x <<  4)) & 0x030c30c3;
    x = (x ^ (x <<  2)) & 0x09249249;
    return x;
}

inline unsigned zm_compact1by2(unsigned x)
{
    x &= 0x09249249;
    x = (x ^ (x >>  2)) & 0x030c30c3;
    x = (x ^ (x >>  4)) & 0x0300f00f;
    x = (x ^ (x >>  8)) & 0xff0000ff;
    x = (x ^ (x >> 16)) & 0x000003ff;
    return x;
}

inline unsigned zm_encode(unsigned x, unsigned y, unsigned z)
{
    return (zm_part1by2(z) << 2) + (zm_part1by2(y) << 1) + zm_part1by2(x);
}

inline void zm_decode(unsigned code, unsigned* x, unsigned* y, unsigned* z)
{
    *x = zm_compact1by2(code);
    *y = zm_compact1by2(code >> 1);
    *z = zm_compact1by2(code >> 2);
}

/*@q*/
#endif /* ZMORTON_H */
