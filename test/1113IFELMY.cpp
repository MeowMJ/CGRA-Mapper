/*
 * Filename: IFELtest.cpp
 * Author: Miaomiao Jiang
 * Date: 2024-10-22
 * Description: This file tests if-else strucure under CGRA-Mapper.
 * Version: 1.0
 */

// nwkernel
#define ALEN 128
#define BLEN 128
#define ALIGN '\\'
#define SKIPA '^'
#define SKIPB '<'
// word2vec: CreateBinaryTree
#define MAX_CODE_LENGTH 40
long long vocab_size = 0;
// dwt::wtree_per
typedef struct wt_set *wt_object;
wt_object wt_init(int wave, const char *method, int siglength, int J);
struct wt_set
{
  int wave;
  int cobj;
  char method[10];
  int siglength;      // Length of the original signal.
  int modwtsiglength; // Modified signal length for MODWT
  int outlength;      // Length of the output DWT vector
  int lenlength;      // Length of the Output Dimension Vector "length"
  int J;              // Number of decomposition Levels
  int MaxIter;        // Maximum Iterations J <= MaxIter
  int even;           // even = 1 if signal is of even length. even = 0 otherwise
  char ext[10];       // Type of Extension used - "per" or "sym"
  char cmethod[10];   // Convolution Method - "direct" or "FFT"

  int N; //
  int cfftset;
  int zpad;
  int length[102];
  double *output;
  double params[0];
};

int main()
{
  return 0;
}

// adpcm_coder
void adpcm_coder(unsigned char indata[], char outdata[], int len, const int valprev)
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

// nw: needwun
void nwkernel(char SEQA[ALEN], char SEQB[BLEN],
              char alignedA[ALEN + BLEN], char alignedB[ALEN + BLEN],
              int M[(ALEN + 1) * (BLEN + 1)], char ptr[(ALEN + 1) * (BLEN + 1)])
{
  int r;
  int a_idx, b_idx;
  int a_str_idx, b_str_idx;

  // TraceBack (n.b. aligned sequences are backwards to avoid string appending)
  a_idx = ALEN;
  b_idx = BLEN;
  a_str_idx = 0;
  b_str_idx = 0;

  while (a_idx > 0 || b_idx > 0)
  {
    r = b_idx * (ALEN + 1);
    if (ptr[r + a_idx] == ALIGN)
    {
      alignedA[a_str_idx++] = SEQA[a_idx - 1];
      alignedB[b_str_idx++] = SEQB[b_idx - 1];
      a_idx--;
      b_idx--;
    }
    else if (ptr[r + a_idx] == SKIPB)
    {
      alignedA[a_str_idx++] = SEQA[a_idx - 1];
      alignedB[b_str_idx++] = '-';
      a_idx--;
    }
    else
    { // SKIPA
      alignedA[a_str_idx++] = '-';
      alignedB[b_str_idx++] = SEQB[b_idx - 1];
      b_idx--;
    }
  }
}

// word2vec: CreateBinaryTree
void CreateBinaryTree()
{
  long long a, min1i, min2i, pos1, pos2;
  long long *count = new long long[vocab_size * 2 + 1];
  long long *binary = new long long[vocab_size * 2 + 1];
  long long *parent_node = new long long[vocab_size * 2 + 1];

  pos1 = vocab_size - 1;
  pos2 = vocab_size;
  // Following algorithm constructs the Huffman tree by adding one node at a time
  for (a = 0; a < vocab_size - 1; a++)
  {
    // First, find two smallest nodes 'min1, min2'
    if (pos1 >= 0)
    {
      if (count[pos1] < count[pos2])
      {
        min1i = pos1;
        pos1--;
      }
      else
      {
        min1i = pos2;
        pos2++;
      }
    }
    else
    {
      min1i = pos2;
      pos2++;
    }
    if (pos1 >= 0)
    {
      if (count[pos1] < count[pos2])
      {
        min2i = pos1;
        pos1--;
      }
      else
      {
        min2i = pos2;
        pos2++;
      }
    }
    else
    {
      min2i = pos2;
      pos2++;
    }
    count[vocab_size + a] = count[min1i] + count[min2i];
    parent_node[min1i] = vocab_size + a;
    parent_node[min2i] = vocab_size + a;
    binary[min2i] = 1;
  }
}

// word2vec: CBT is simpled CreateBinaryTree, in case CreateBinaryTree is too hard to map
void CBT()
{
  long long *count = new long long[vocab_size * 2 + 1];
  long long *binary = new long long[vocab_size * 2 + 1];
  long long *parent_node = new long long[vocab_size * 2 + 1];
  long long a, min1i, min2i, pos1, pos2;
  pos1 = vocab_size - 1;
  pos2 = vocab_size;
  for (a = 0; a < vocab_size - 1; a++)
  {
    if ((pos1 >= 0) && (count[pos1] < count[pos2]))
    {
      min1i = pos1;
      pos1--;
    }
    else
    {
      min1i = pos2;
      pos2++;
    }
    if ((pos1 >= 0) && (count[pos1] < count[pos2]))
    {
      min2i = pos1;
      pos1--;
    }
    else
    {
      min2i = pos2;
      pos2++;
    }
    count[vocab_size + a] = count[min1i] + count[min2i];
    parent_node[min1i] = vocab_size + a;
    parent_node[min2i] = vocab_size + a;
    binary[min2i] = 1;
  }
}

// dwt::wtree_per
void wtree_per(wt_object wt, double *inp, int N, double *cA, int len_cA, double *cD)
{
  int l, l2, isodd, i, t, len_avg;

  len_avg = 128;
  l2 = len_avg / 2;
  isodd = N % 2;

  for (i = 0; i < len_cA; ++i)
  {
    t = 2 * i + l2;
    cA[i] = 0.0;
    cD[i] = 0.0;
    for (l = 0; l < len_avg; ++l)
    {
      if ((t - l) >= l2 && (t - l) < N)
      {
        cA[i] += wt->wave * inp[t - l];
        cD[i] += wt->wave * inp[t - l];
      }
      else if ((t - l) < l2 && (t - l) >= 0)
      {
        cA[i] += wt->wave * inp[t - l];
        cD[i] += wt->wave * inp[t - l];
      }
      else if ((t - l) < 0 && isodd == 0)
      {
        cA[i] += wt->wave * inp[t - l + N];
        cD[i] += wt->wave * inp[t - l + N];
      }
      else if ((t - l) < 0 && isodd == 1)
      {
        if ((t - l) != -1)
        {
          cA[i] += wt->wave * inp[t - l + N + 1];
          cD[i] += wt->wave * inp[t - l + N + 1];
        }
        else
        {
          cA[i] += wt->wave * inp[N - 1];
          cD[i] += wt->wave * inp[N - 1];
        }
      }
      else if ((t - l) >= N && isodd == 0)
      {
        cA[i] += wt->wave * inp[t - l - N];
        cD[i] += wt->wave * inp[t - l - N];
      }
      else if ((t - l) >= N && isodd == 1)
      {
        if (t - l != N)
        {
          cA[i] += wt->wave * inp[t - l - (N + 1)];
          cD[i] += wt->wave * inp[t - l - (N + 1)];
        }
        else
        {
          cA[i] += wt->wave * inp[N - 1];
          cD[i] += wt->wave * inp[N - 1];
        }
      }
    }
  }
}
