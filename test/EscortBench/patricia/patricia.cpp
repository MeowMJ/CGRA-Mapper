#define PI 3.14159265358979323846
#include <math.h>
#include <stdlib.h>
// patricia/insertR
int __attribute__((inline)) bit(int i, unsigned long key)
{
  return key & (1 << (31 - i));
}

// patricia/insertR
int __attribute__((noinline)) insertR_kernel(int h, int n, int d, int p)
{
  if ((h >= d) || (h <= p))
  {
    n = d;
    n = bit(d, n) ? h : n;
    return n;
  }

  if (bit(h, n))
    h += insertR(h, n, d, h);
  else
    h -= insertR(h, n, d, h);
  return h;
}
