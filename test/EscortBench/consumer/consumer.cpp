// consumer/lame/lame3.70/fft.c
# define        SQRT2                   1.41421356237309504880
#define TRI_SIZE (5-1) /* 1024 =  4**5 */
static float costab[TRI_SIZE*2];
void fht(float *fz, short n)
{
    short k4;
    float *fi, *fn, *gi;
    float *tri;

    fn = fz + n;
    tri = &costab[0];
    k4 = 4;
    do {
	float s1, c1;
	short i, k1, k2, k3, kx;
	kx  = k4 >> 1;
	k1  = k4;
	k2  = k4 << 1;
	k3  = k2 + k1;
	k4  = k2 << 1;
	fi  = fz;
	gi  = fi + kx;
	do {
	    float f0,f1,f2,f3;
	    f1      = fi[0]  - fi[k1];
	    f0      = fi[0]  + fi[k1];
	    f3      = fi[k2] - fi[k3];
	    f2      = fi[k2] + fi[k3];
	    fi[k2]  = f0     - f2;
	    fi[0 ]  = f0     + f2;
	    fi[k3]  = f1     - f3;
	    fi[k1]  = f1     + f3;
	    f1      = gi[0]  - gi[k1];
	    f0      = gi[0]  + gi[k1];
	    f3      = SQRT2  * gi[k3];
	    f2      = SQRT2  * gi[k2];
	    gi[k2]  = f0     - f2;
	    gi[0 ]  = f0     + f2;
	    gi[k3]  = f1     - f3;
	    gi[k1]  = f1     + f3;
	    gi     += k4;
	    fi     += k4;
	} while (fi<fn);
	c1 = tri[0];
	s1 = tri[1];
	for (i = 1; i < kx; i++) {
	    float c2,s2;
	    c2 = 1 - (2*s1)*s1;
	    s2 = (2*s1)*c1;
	    fi = fz + i;
	    gi = fz + k1 - i;
	    do {
		float a,b,g0,f0,f1,g1,f2,g2,f3,g3;
		b       = s2*fi[k1] - c2*gi[k1];
		a       = c2*fi[k1] + s2*gi[k1];
		f1      = fi[0 ]    - a;
		f0      = fi[0 ]    + a;
		g1      = gi[0 ]    - b;
		g0      = gi[0 ]    + b;
		b       = s2*fi[k3] - c2*gi[k3];
		a       = c2*fi[k3] + s2*gi[k3];
		f3      = fi[k2]    - a;
		f2      = fi[k2]    + a;
		g3      = gi[k2]    - b;
		g2      = gi[k2]    + b;
		b       = s1*f2     - c1*g3;
		a       = c1*f2     + s1*g3;
		fi[k2]  = f0        - a;
		fi[0 ]  = f0        + a;
		gi[k3]  = g1        - b;
		gi[k1]  = g1        + b;
		b       = c1*g2     - s1*f3;
		a       = s1*g2     + c1*f3;
		gi[k2]  = g0        - a;
		gi[0 ]  = g0        + a;
		fi[k3]  = f1        - b;
		fi[k1]  = f1        + b;
		gi     += k4;
		fi     += k4;
	    } while (fi<fn);
	    c2 = c1;
	    c1 = c2 * tri[0] - s1 * tri[1];
	    s1 = c2 * tri[1] + s1 * tri[0];
        }
	tri += 2;
    } while (k4<n);
}

typedef   double mad_fixed_t;

#define MAD_F_FRACBITS		28
#define mad_f_todouble(x)	((double)  \
				 ((x) / (double) (1L << MAD_F_FRACBITS)))

#define MAD_F(x)		mad_f_todouble(x)
mad_fixed_t const window_s[12] = {
  MAD_F(0x0216a2a2) /* 0.130526192 */, MAD_F(0x061f78aa) /* 0.382683432 */,
  MAD_F(0x09bd7ca0) /* 0.608761429 */, MAD_F(0x0cb19346) /* 0.793353340 */,
  MAD_F(0x0ec835e8) /* 0.923879533 */, MAD_F(0x0fdcf549) /* 0.991444861 */,
  MAD_F(0x0fdcf549) /* 0.991444861 */, MAD_F(0x0ec835e8) /* 0.923879533 */,
  MAD_F(0x0cb19346) /* 0.793353340 */, MAD_F(0x09bd7ca0) /* 0.608761429 */,
  MAD_F(0x061f78aa) /* 0.382683432 */, MAD_F(0x0216a2a2) /* 0.130526192 */,
};


#define MAX_COMPONENTS  10
unsigned char buffer[MAX_COMPONENTS];
#define ms00(f) (window_s[i] * f(i + k))
#define ms10(f) (window_s[0x7f - i] * f(i + k + 0x80))
#define ms20(f) (window_s[i + 0x40] * f(i + k + 0x40))
#define ms30(f) (window_s[0x3f - i] * f(i + k + 0xc0))

#define ms01(f) (window_s[i + 0x01] * f(i + k + 0x01))
#define ms11(f) (window_s[0x7e - i] * f(i + k + 0x81))
#define ms21(f) (window_s[i + 0x41] * f(i + k + 0x41))
#define ms31(f) (window_s[0x3e - i] * f(i + k + 0xc1))
#define ch01(index)  (buffer[chn][index])
#define BLKSIZE_s 128
#define ch2(index)  (((float)(0.5*SQRT2))*(buffer[0][index] + buffer[1][index]))
#define ch2(index)  (((float)(0.5*SQRT2))*(buffer[0][index] + buffer[1][index]))
#define ch3(index)  (((float)(0.5*SQRT2))*(buffer[0][index] - buffer[1][index]))
void fft_short_kernel(
    float x_real[3][BLKSIZE_s], int chn, short *buffer[2])
{
  short i, j, b;

  float *x = &x_real[3][BLKSIZE_s / 2];
  short k = (576 / 3) * (3 + 1);
  j = BLKSIZE_s / 8 - 1;
  if (chn < 2)
  {
    do
    {
      float f0, f1, f2, f3, w;

      // i = rv_tbl[j << 2];

      f0 = ms00(ch01);
      w = ms10(ch01);
      f1 = f0 - w;
      f0 = f0 + w;
      f2 = ms20(ch01);
      w = ms30(ch01);
      f3 = f2 - w;
      f2 = f2 + w;

      x -= 4;
      x[0] = f0 + f2;
      x[2] = f0 - f2;
      x[1] = f1 + f3;
      x[3] = f1 - f3;

      f0 = ms01(ch01);
      w = ms11(ch01);
      f1 = f0 - w;
      f0 = f0 + w;
      f2 = ms21(ch01);
      w = ms31(ch01);
      f3 = f2 - w;
      f2 = f2 + w;

      x[BLKSIZE_s / 2 + 0] = f0 + f2;
      x[BLKSIZE_s / 2 + 2] = f0 - f2;
      x[BLKSIZE_s / 2 + 1] = f1 + f3;
      x[BLKSIZE_s / 2 + 3] = f1 - f3;
    } while (--j >= 0);
  }
  else if (chn == 2)
  {
    do
    {
      float f0, f1, f2, f3, w;

      // i = rv_tbl[j << 2];

      f0 = ms00(ch2);
      w = ms10(ch2);
      f1 = f0 - w;
      f0 = f0 + w;
      f2 = ms20(ch2);
      w = ms30(ch2);
      f3 = f2 - w;
      f2 = f2 + w;

      x -= 4;
      x[0] = f0 + f2;
      x[2] = f0 - f2;
      x[1] = f1 + f3;
      x[3] = f1 - f3;

      f0 = ms01(ch2);
      w = ms11(ch2);
      f1 = f0 - w;
      f0 = f0 + w;
      f2 = ms21(ch2);
      w = ms31(ch2);
      f3 = f2 - w;
      f2 = f2 + w;

      x[BLKSIZE_s / 2 + 0] = f0 + f2;
      x[BLKSIZE_s / 2 + 2] = f0 - f2;
      x[BLKSIZE_s / 2 + 1] = f1 + f3;
      x[BLKSIZE_s / 2 + 3] = f1 - f3;
    } while (--j >= 0);
  }
  else
  {
    do
    {
      float f0, f1, f2, f3, w;

      // i = rv_tbl[j << 2];

      f0 = ms00(ch3);
      w = ms10(ch3);
      f1 = f0 - w;
      f0 = f0 + w;
      f2 = ms20(ch3);
      w = ms30(ch3);
      f3 = f2 - w;
      f2 = f2 + w;

      x -= 4;
      x[0] = f0 + f2;
      x[2] = f0 - f2;
      x[1] = f1 + f3;
      x[3] = f1 - f3;

      f0 = ms01(ch3);
      w = ms11(ch3);
      f1 = f0 - w;
      f0 = f0 + w;
      f2 = ms21(ch3);
      w = ms31(ch3);
      f3 = f2 - w;
      f2 = f2 + w;

      x[BLKSIZE_s / 2 + 0] = f0 + f2;
      x[BLKSIZE_s / 2 + 2] = f0 - f2;
      x[BLKSIZE_s / 2 + 1] = f1 + f3;
      x[BLKSIZE_s / 2 + 3] = f1 - f3;
    } while (--j >= 0);
  }

  // fht(x, BLKSIZE_s);
}


int ulstrcmp(char const *str1, char const *str2)
{
    register char c1,c2;

    for(;;) {
	c1 = *str1++;
	if (c1 <= 'Z')
	    if (c1 >= 'A')
		c1 += 040;
	c2 = *str2++;
	if (c2 <= 'Z')
	    if (c2 >= 'A')
		c2 += 040;
	if (c1 != c2)
	    break;
	if (c1 == '\0')
	    return(0);
    }
    return(c1 - c2);
}