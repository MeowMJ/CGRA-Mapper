typedef struct wt_set *wt_object;

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

wt_object wt_init(int wave, const char *method, int siglength, int J);

// wtree_per_kernel
void wtree_per_kernel(wt_object wt, double *inp, int N, double *cA, int len_cA, double *cD)
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
