#ifndef __E_FACE_H__
#define __E_FACE_H__

#include <boost/array.hpp>

namespace LBIE
{

typedef struct _edge_face {
  int faceidx;
  //int orientation[4];
  boost::array<int,4> orientation;
  int facebit;
} edge_face;

}

#endif

