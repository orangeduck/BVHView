#ifndef COMMON_MACROS
#define COMMON_MACROS

#include <limits.h>

#ifdef _WIN32
    #define realpath(N,R) _fullpath((R),(N),_MAX_PATH)

    #define popen _popen
    #define pclose _pclose

    #undef PATH_MAX
    #define PATH_MAX _MAX_PATH
#endif

#endif
