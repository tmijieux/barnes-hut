#ifndef TDP_UTIL_H
#define TDP_UTIL_H

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#define ASSERT_MSG(cond, msg) assert(  ((void)(msg), (cond)) )

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifdef PASTE
#undef PASTE
#endif
#define PASTE(X, Y) PASTE2(X, Y)
#define PASTE2(X, Y) X ## Y

#ifdef DEQUAL
#undef DEQUAL
#endif
#define DEQUAL(X_, Y_, EPS_) (fabs((X_) - (Y_)) < (EPS_))

#endif // TDP_UTIL_H
