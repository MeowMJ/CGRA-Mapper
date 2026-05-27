// adpcm_coder_kernel
void adpcm_coder_kernel(unsigned char indata[], char outdata[], int len, const int valprev)
{
  unsigned char *inp;  /* Input buffer pointer */
  unsigned char *outp; /* output buffer pointer */
  int val;             /* Current input sample value */
  int sign;            /* Current adpcm sign bit */
  int delta;           /* Current adpcm output value */
  int diff;            /* Difference between val and valprev */
  int step;            /* Stepsize */
  int valpred;         /* Predicted output value */
  int vpdiff;          /* Current change to valpred */
  int outputbuffer;    /* place to keep previous 4-bit value */
  int bufferstep;      /* toggle between outputbuffer/output */

  outp = (unsigned char *)outdata;
  inp = indata;
  valpred = valprev;
  step = 9;

  bufferstep = 1;

  for (; len > 0; len--)
  {
    val = *inp++;

    /* Step 1 - compute difference with previous value */
    diff = val - valpred;
    sign = (diff < 0) ? 8 : 0;
    if (sign)
      diff = (-diff);

    /* Step 2 - Divide and clamp */
    delta = 0;
    vpdiff = (step >> 3);

    if (diff >= step)
    {
      delta = 4;
      diff -= step;
      vpdiff += step;
    }
    step >>= 1;
    if (diff >= step)
    {
      delta |= 2;
      diff -= step;
      vpdiff += step;
    }
    step >>= 1;
    if (diff >= step)
    {
      delta |= 1;
      vpdiff += step;
    }

    /* Step 3 - Update previous value */
    if (sign)
      valpred -= vpdiff;
    else
      valpred += vpdiff;

    /* Step 4 - Clamp previous value to 16 bits */
    if (valpred > 32767)
      valpred = 32767;
    else if (valpred < -32768)
      valpred = -32768;

    /* Step 5 - Assemble value, update index and step values */

    /* Step 6 - Output value */
    if (bufferstep)
    {
      outputbuffer = (delta << 4) & 0xf0;
    }
    else
    {
      *outp++ = (delta & 0x0f) | outputbuffer;
    }
    bufferstep = !bufferstep;
  }
}
