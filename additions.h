#include "raylib.h"

// Function specifiers in case library is build/used as a shared library (Windows)
// NOTE: Microsoft specifiers to tell compiler that symbols are imported/exported from a .dll
#if defined(_WIN32)
    #if defined(BUILD_LIBTYPE_SHARED)
        #if defined(__TINYC__)
            #define __declspec(x) __attribute__((x))
        #endif
        #define BVHA __declspec(dllexport)     // We are building the library as a Win32 shared library (.dll)
    #elif defined(USE_LIBTYPE_SHARED)
        #define BVHA __declspec(dllimport)     // We are using the library as a Win32 shared library (.dll)
    #endif
#endif

#ifndef BVHA
    #define BVHA       // Functions defined as 'extern' by default (implicit specifiers)
#endif

typedef struct {
    Model model;
    bool isLoaded;
} CharacterModel;

BVHA bool BVHALoadCharacterModelFromFile(CharacterModel* characterModel, const char *fileName);
BVHA void BVHAUnloadCharacterModel(CharacterModel* characterModel);
BVHA bool boneExistsInModel(const char* boneName, const Model* model);
