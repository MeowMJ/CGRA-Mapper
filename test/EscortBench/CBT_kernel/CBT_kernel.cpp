long long vocab_size = 0;

// CBT_kernel
void CBT_kernel()
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
