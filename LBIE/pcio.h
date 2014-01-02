#ifndef _PC_IO_H
#define _PC_IO_H

#include <stdio.h>

namespace LBIE
{

size_t getFloat(float *, size_t, FILE *);
size_t getInt(int *, size_t, FILE *);
size_t getShort(short *, size_t, FILE *);
size_t getUnChar(unsigned char *, size_t, FILE *);
size_t putFloat(float*,size_t,FILE*);

}


#endif
