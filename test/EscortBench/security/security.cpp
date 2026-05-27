#define PI 3.14159265358979323846
#include <math.h>
#include <stdlib.h>
// rijndael, aes.c
#define aes_good    1
typedef short           aes_ret;
#define cf_dec  aes_ret
#define word_out(x,v)   *(word*)(x) = (v)
#define so(y,x,c)   word_out(y + 4 * c, s(x,c))
#define state_out(y,x)  so(y,x,0); so(y,x,1); so(y,x,2); so(y,x,3)
typedef unsigned char   byte;           /* must be an 8-bit storage unit */
typedef unsigned long   word;           /* must be a 32-bit storage unit */
#define BLOCK_SIZE  16
#define RC_LENGTH   5 * BLOCK_SIZE / 4 - (BLOCK_SIZE == 16 ? 10 : 11)
#define KS_LENGTH   4 * BLOCK_SIZE
#define locals(y,x)     x[4],y[4]
#define word_in(x)      *(word*)(x)
#define s(x,c) x[c]
#define si(y,x,k,c) s(y,c) = word_in(x + 4 * c) ^ k[c]
#define astate_in(y,x,k) si(y,x,k,0); si(y,x,k,1); si(y,x,k,2); si(y,x,k,3)
#define Ncol   BLOCK_SIZE / 4
#define nc   Ncol
#define aes_bad     0
#define bval(x,n)       ((byte*)((x) >> 8 * (n)))
#define four_tables(x,tab,vf,rf,c)  (  tab[0][bval(vf(x,0,c),rf(0,c))]   ^ tab[1][bval(vf(x,1,c),rf(1,c))]   ^ tab[2][bval(vf(x,2,c),rf(2,c))]   ^ tab[3][bval(vf(x,3,c),rf(3,c))])
#define fwd_lrnd(y,x,k,c)   s(y,c)= (k)[c] ^ four_tables(x,fl_tab,fwd_var,rf1,c)
static word  ft_tab[256];
static word  it_tab[256];
#define fwd_var(x,r,c)  ( r==0 ?               ( c==0 ? s(x,0)     : c==1 ? s(x,1)     : c==2 ? s(x,2)     : c==3 ? s(x,3)     : c==4 ? s(x,4)     : c==5 ? s(x,5)     : c==6 ? s(x,6)     : s(x,7))        : r==1 ?               ( c==0 ? s(x,1)     : c==1 ? s(x,2)     : c==2 ? s(x,3)     : c==3 ? nc==4 ? s(x,0) : s(x,4)     : c==4 ? s(x,5)     : c==5 ? nc==8 ? s(x,6) : s(x,0)     : c==6 ? s(x,7)     : s(x,0))        : r==2 ?               ( c==0 ? nc==8 ? s(x,3) : s(x,2)     : c==1 ? nc==8 ? s(x,4) : s(x,3)     : c==2 ? nc==4 ? s(x,0) : nc==8 ? s(x,5) : s(x,4)     : c==3 ? nc==4 ? s(x,1) : nc==8 ? s(x,6) : s(x,5)     : c==4 ? nc==8 ? s(x,7) : s(x,0)     : c==5 ? nc==8 ? s(x,0) : s(x,1)     : c==6 ? s(x,1)     : s(x,2))        :                      ( c==0 ? nc==8 ? s(x,4) : s(x,3)     : c==1 ? nc==4 ? s(x,0) : nc==8 ? s(x,5) : s(x,4)     : c==2 ? nc==4 ? s(x,1) : nc==8 ? s(x,6) : s(x,5)     : c==3 ? nc==4 ? s(x,2) : nc==8 ? s(x,7) : s(x,0)     : c==4 ? nc==8 ? s(x,0) : s(x,1)     : c==5 ? nc==8 ? s(x,1) : s(x,2)     : c==6 ? s(x,2)     : s(x,3)))
#define rf1(r,c)    (r)
#define fwd_rnd(y,x,k,c)    s(y,c)= (k)[c] ^ four_tables(x,ft_tab,fwd_var,rf1,c)
typedef struct
{
    word    Nkey;               /* the number of words in the key input block */
    word    Nrnd;               /* the number of cipher rounds                */
    word    e_key[KS_LENGTH];   /* the encryption key schedule                */
    word    d_key[KS_LENGTH];   /* the decryption key schedule                */
#if !defined(BLOCK_SIZE)
    word    Ncol;               /* the number of columns in the cipher state  */
#endif
    byte    mode;               /* encrypt, decrypt or both                   */
}aes;
#define	MAX_CMAP_SIZE	256


word  fl_tab[256];
#define c_name(x)   x
//#define round(rm,y,x,k) rm(y,x,k,0); rm(y,x,k,1); rm(y,x,k,2); rm(y,x,k,3)
void __attribute__((inline)) rount(){

}

cf_dec c_name(encrypt_kernel)(const byte in_blk[], byte out_blk[], const aes * cx)
{
  word locals(b0, b1);
  const word *kp = cx->e_key;

  if (!(cx->mode & 0x01))
    return aes_bad;
  astate_in(b0, in_blk, kp);
  kp += nc;

  switch (cx->Nrnd)
  {
  case 14:
 rount();
 rount();
    kp += 2 * nc;
  case 12:

 rount();
 rount();
    kp += 2 * nc;
  case 10:

 rount();
 rount();
 rount();
 rount();
 rount();
 rount();
 rount();
 rount();
 rount();
 rount();
  }

  state_out(out_blk, b0);
  return aes_good;
}

