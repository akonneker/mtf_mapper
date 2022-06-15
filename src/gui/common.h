#ifndef COMMON_H
#define COMMON_H

#ifndef _WIN32
    // assume non-win32 means gcc. not true in general
    #ifndef ATTRIBUTE_UNUSED
        #define ATTRIBUTE_UNUSED __attribute__ ((__unused__))
    #endif /* ATTRIBUTE_UNUSED */
#else
    #define ATTRIBUTE_UNUSED
#endif

#endif

