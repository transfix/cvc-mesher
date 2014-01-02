#ifndef __ENDIANS_H__
#define __ENDIANS_H__

#define SWAP_64(a) \
  { \
    unsigned char tmp[8]; \
    unsigned char *ch; \
    ch = (unsigned char *)(a); \
    tmp[0] = ch[0]; tmp[1] = ch[1]; tmp[2] = ch[2]; tmp[3] = ch[3]; \
    tmp[4] = ch[4]; tmp[5] = ch[5]; tmp[6] = ch[6]; tmp[7] = ch[7]; \
    ch[0] = tmp[7]; ch[1] = tmp[6]; ch[2] = tmp[5]; ch[3] = tmp[4]; \
    ch[4] = tmp[3]; ch[5] = tmp[2]; ch[6] = tmp[1]; ch[7] = tmp[0]; \
  }
#define SWAP_32(a) \
  { \
    unsigned char tmp[4]; \
    unsigned char *ch; \
    ch = (unsigned char *)(a); \
    tmp[0] = ch[0]; tmp[1] = ch[1]; tmp[2] = ch[2]; tmp[3] = ch[3]; \
    ch[0] = tmp[3]; ch[1] = tmp[2]; ch[2] = tmp[1]; ch[3] = tmp[0]; \
  }
#define SWAP_16(a) \
  { \
    unsigned char d; \
    unsigned char *ch; \
    ch = (unsigned char *)(a); \
    d = ch[0]; \
    ch[0] = ch[1]; \
    ch[1] = d; \
  }

static inline int big_endian()
{
  long one=1;
  return !(*((char *)(&one)));
}

#endif
