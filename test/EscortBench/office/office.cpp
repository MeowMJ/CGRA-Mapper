#define PI 3.14159265358979323846
#include <math.h>
#include <stdlib.h>
// office/stringsearch/pbmsrch_large.c
#define UCHAR_MAX	255
typedef long unsigned int size_t;
#include <string.h>
#include <ctype.h>
typedef unsigned char uchar;


#define LARGE 32767

static int patlen;
static int skip[UCHAR_MAX+1];                               /* rdg 10/93 */
static int skip2;
static uchar *pat;

void bmh_init(const char *pattern)
{
          int i, lastpatchar;

          pat = (uchar *)pattern;
          patlen = strlen(pattern);
          for (i = 0; i <= UCHAR_MAX; ++i)                  /* rdg 10/93 */
                skip[i] = patlen;
          for (i = 0; i < patlen; ++i)
                skip[pat[i]] = patlen - i - 1;
          lastpatchar = pat[patlen - 1];
          skip[lastpatchar] = LARGE;
          skip2 = patlen;                 /* Horspool's fixed second shift */
          for (i = 0; i < patlen - 1; ++i)
          {
                if (pat[i] == lastpatchar)
                      skip2 = patlen - i - 1;
          }
}

char *bmh_search_kernel(const char *string, const int stringlen)
{
      int i, j;
      char *s;

      i = patlen - 1 - stringlen;
      if (i >= 0)
            return NULL;
      string += stringlen;
      for ( ;; )
      {
            while ( (i += skip[((uchar *)string)[i]]) < 0 )
                  ;                           /* mighty fast inner loop */
            if (i < (LARGE - stringlen))
                  return NULL;
            i -= LARGE;
            j = patlen - 1;
            s = (char *)string + (i - j);
            while (--j >= 0 && s[j] == pat[j])
                  ;
            if ( j < 0 )                                    /* rdg 10/93 */
                  return s;                                 /* rdg 10/93 */
            if ( (i += skip2) >= 0 )                        /* rdg 10/93 */
                  return NULL;                              /* rdg 10/93 */
      }
}
