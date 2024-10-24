#ifndef COMMON_MACROS
#define COMMON_MACROS

#include <limits.h>

#ifdef _WIN32
    #define realpath(N,R) _fullpath((R),(N),_MAX_PATH)

    #undef PATH_MAX
    #define PATH_MAX _MAX_PATH
#endif

#ifdef _WIN32
    #define mkdir_p(dir) mkdir(dir)
    #define PATH_SEPARATOR '\\'
#else
    #define mkdir_p(dir) mkdir(dir, 0755)
    #define PATH_SEPARATOR '/'
#endif

#define WAIT_FOR_INPUT_ON_EXIT false // Useful to inspect console window output before it closes

#endif
