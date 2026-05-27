#define ALEN 128
#define BLEN 128
#define ALIGN '\\'
#define SKIPA '^'
#define SKIPB '<'

// nwkernel
void nwkernel_kernel(char SEQA[ALEN], char SEQB[BLEN],
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
