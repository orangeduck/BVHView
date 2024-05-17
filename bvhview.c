/*******************************************************************************************
*
*    BVHView - A simple BVH animation viewer written using raylib
*
*  This is a simple viewer for the .bvh animation file format made using raylib. For more
*  info on the motivation behind it and information on features and documentation please 
*  see: https://theorangeduck.com/page/bvhview
*
*  The program itself essentially consists of the following components:
*
*     - A parser for the BVH file format
*     - A set of functions for sampling data from particular frames of the BVH file.
*     - A set of functions for creating capsules from the skeleton structure of the BVH data 
*       and animation transforms.
*     - A (relatively) efficient and high quality shader for rendering capsules that includes 
*       nice lighting, soft shadows, and some CPU based culling to limit the amount of work
*       required by the GPU.
*
*  Coding style is roughly meant to follow the rest of raylib and community contributions
*  are very welcome.
*
*******************************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <float.h>
#include <errno.h>

#include "raylib.h"
#include "rcamera.h"
#include "raymath.h"
#include "rlgl.h"
#define RAYGUI_WINDOWBOX_STATUSBAR_HEIGHT 24
#define GUI_WINDOW_FILE_DIALOG_IMPLEMENTATION
#include "../examples/custom_file_dialog/gui_window_file_dialog.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#include "additions.h"
#include "rmodels.c"

#if defined(PLATFORM_WEB)
#include <emscripten/emscripten.h>
#endif

//----------------------------------------------------------------------------------
// Profiling
//----------------------------------------------------------------------------------

// Un-comment to enable profiling
//#define ENABLE_PROFILE

// Profiling only available on Windows
#if defined(ENABLE_PROFILE) && defined(_WIN32)

#include <profileapi.h>

enum
{
    // Max number of profile records (profiled code locations) 
    PROFILE_RECORD_MAX = 512,
    
    // Maximum number of timer samples per record
    PROFILE_RECORD_SAMPLE_MAX = 128,
};

// A single record for a profiled code location with a cyclic buffer of start and end times.
typedef struct 
{
    const char* name;
    uint32_t idx;
    uint32_t num;
    
    struct {
        LARGE_INTEGER start;
        LARGE_INTEGER end;
    } samples[PROFILE_RECORD_SAMPLE_MAX];
    
} ProfileRecord;

// Structure containing space for all profiled code locations
typedef struct
{
    uint32_t num;
    LARGE_INTEGER freq;
    ProfileRecord* records[PROFILE_RECORD_MAX];
    
} ProfileRecordData;

// Global variable storing all the profile record data
static ProfileRecordData globalProfileRecords;

// Init the Profile Record Data. Must be called at program start
static void ProfileRecordDataInit()
{
    globalProfileRecords.num = 0;
    QueryPerformanceFrequency(&globalProfileRecords.freq);
    memset(globalProfileRecords.records, 0, sizeof(ProfileRecord*) * PROFILE_RECORD_MAX);
}

// If uninitialized, then initialize the profile record, then store the start time
static inline void ProfileRecordBegin(ProfileRecord* record, const char* name)
{
    if (!record->name && globalProfileRecords.num < PROFILE_RECORD_MAX)
    {
        record->name = name;
        record->idx = 0;
        record->num = 0;
        globalProfileRecords.records[globalProfileRecords.num] = record;
        globalProfileRecords.num++;
    }
  
    QueryPerformanceCounter(&record->samples[record->idx].start);
}

// Store the end time and increment the record sample num
static inline void ProfileRecordEnd(ProfileRecord* record)
{
    QueryPerformanceCounter(&record->samples[record->idx].end);
    record->idx = (record->idx + 1) % PROFILE_RECORD_SAMPLE_MAX;
    record->num++;
}

// Tickers record a rolling average of Profile Record durations in microseconds
typedef struct
{
    uint64_t unitScale; 
    double alpha;
    uint32_t samples[PROFILE_RECORD_MAX];
    uint64_t iterations[PROFILE_RECORD_MAX];
    double averages[PROFILE_RECORD_MAX];
    double times[PROFILE_RECORD_MAX];
    
} ProfileTickers;

// Global profile tickers data
static ProfileTickers globalProfileTickers;

// Initialize ticker data
static inline void ProfileTickersInit()
{
    globalProfileTickers.unitScale = 1000000; // Microseconds
    globalProfileTickers.alpha = 0.9f;
    memset(globalProfileTickers.samples, 0, sizeof(uint32_t) * PROFILE_RECORD_MAX);
    memset(globalProfileTickers.iterations, 0, sizeof(uint64_t) * PROFILE_RECORD_MAX);
    memset(globalProfileTickers.averages, 0, sizeof(double) * PROFILE_RECORD_MAX);
    memset(globalProfileTickers.times, 0, sizeof(double) * PROFILE_RECORD_MAX);
}

// Update tickers and compute the rolling average of the duration
static inline void ProfileTickersUpdate()
{
    for (int i = 0; i < globalProfileRecords.num; i++)
    {
        ProfileRecord* record = globalProfileRecords.records[i];
        
        if (record && record->name)
        {
            globalProfileTickers.samples[i] = record->num;
            
            int bufferedSampleNum = record->num < PROFILE_RECORD_SAMPLE_MAX ? record->num : PROFILE_RECORD_SAMPLE_MAX;
            
            for (int j = 0; j < bufferedSampleNum; j++)
            {
                double time = (double)((
                    record->samples[j].end.QuadPart - 
                    record->samples[j].start.QuadPart) * globalProfileTickers.unitScale) / 
                        (double)globalProfileRecords.freq.QuadPart;
                    
                globalProfileTickers.iterations[i]++;
                globalProfileTickers.averages[i] = globalProfileTickers.alpha * globalProfileTickers.averages[i] + (1.0 - globalProfileTickers.alpha) * time;
                globalProfileTickers.times[i] = globalProfileTickers.averages[i] / (1.0 - pow(globalProfileTickers.alpha, globalProfileTickers.iterations[i]));
            }
            
            // Flush Samples
            record->idx = 0;
            record->num = 0;
        }
    }
}

#define PROFILE_INIT() ProfileRecordDataInit(); 
#define PROFILE_BEGIN(NAME) static ProfileRecord __PROFILE_RECORD_##NAME; ProfileRecordBegin(&__PROFILE_RECORD_##NAME, #NAME);
#define PROFILE_END(NAME) ProfileRecordEnd(&__PROFILE_RECORD_##NAME);

#define PROFILE_TICKERS_INIT() ProfileTickersInit();
#define PROFILE_TICKERS_UPDATE() ProfileTickersUpdate()

#else
#define PROFILE_INIT() 
#define PROFILE_BEGIN(NAME) 
#define PROFILE_END(NAME) 

#define PROFILE_TICKERS_INIT()
#define PROFILE_TICKERS_UPDATE()
#endif

//----------------------------------------------------------------------------------
// Additional Raylib Functions
//----------------------------------------------------------------------------------

static inline float Max(float x, float y)
{
    return x > y ? x : y;
}

static inline float Min(float x, float y)
{
    return x < y ? x : y;
}

static inline float Saturate(float x)
{
    return Clamp(x, 0.0f, 1.0f);
}

static inline float Square(float x)
{
    return x * x;
}

static inline int ClampInt(int x, int min, int max)
{
    return x < min ? min : x > max ? max : x;
}

static inline int MaxInt(int x, int y)
{
    return x > y ? x : y;
}

static inline int MinInt(int x, int y)
{
    return x < y ? x : y;
}

// This is a safe version of QuaternionBetween which returns a 180 deg rotation
// at the singularity where vectors are facing exactly in opposite directions
static inline Quaternion QuaternionBetween(Vector3 p, Vector3 q)
{
    Vector3 c = Vector3CrossProduct(p, q);

    Quaternion o = {
        c.x,
        c.y,
        c.z,
        sqrtf(Vector3DotProduct(p, p) * Vector3DotProduct(q, q)) + Vector3DotProduct(p, q),
    };
    
    return QuaternionLength(o) < 1e-8f ?
        QuaternionFromAxisAngle((Vector3){ 1.0f, 0.0f, 0.0f }, PI) :
        QuaternionNormalize(o);
}

// Puts the quaternion in the hemisphere closest to the identity
static inline Quaternion QuaternionAbsolute(Quaternion q)
{
    if (q.w < 0.0f)
    {
        q.x = -q.x;
        q.y = -q.y;
        q.z = -q.z;
        q.w = -q.w;
    }

    return q;
}

// Quaternion exponent, log, and angle axis functions (see: https://theorangeduck.com/page/exponential-map-angle-axis-angular-velocity)

static inline Quaternion QuaternionExp(Vector3 v)
{
    float halfangle = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);

    if (halfangle < 1e-4f)
    {
        return QuaternionNormalize((Quaternion){ v.x, v.y, v.z, 1.0f });
    }
    else
    {
        float c = cosf(halfangle);
        float s = sinf(halfangle) / halfangle;
        return (Quaternion){ s * v.x, s * v.y, s * v.z, c };
    }
}

static inline Vector3 QuaternionLog(Quaternion q)
{
    float length = sqrtf(q.x*q.x + q.y*q.y + q.z*q.z);

    if (length < 1e-4f)
    {
        return (Vector3){ q.x, q.y, q.z };
    }
    else
    {
        float halfangle = acosf(Clamp(q.w, -1.0f, 1.0f));
        return Vector3Scale((Vector3){ q.x, q.y, q.z }, halfangle / length);
    }
}

static inline Vector3 QuaternionToScaledAngleAxis(Quaternion q)
{
    return Vector3Scale(QuaternionLog(q), 2.0f);
}

static inline Quaternion QuaternionFromScaledAngleAxis(Vector3 v)
{
    return QuaternionExp(Vector3Scale(v, 0.5f));
}

// Cubic Interpolation (see: https://theorangeduck.com/page/cubic-interpolation-quaternions)

static inline Vector3 Vector3Hermite(Vector3 p0, Vector3 p1, Vector3 v0, Vector3 v1, float alpha)
{
    float x = alpha;
    float w0 = 2*x*x*x - 3*x*x + 1;
    float w1 = 3*x*x - 2*x*x*x;
    float w2 = x*x*x - 2*x*x + x;
    float w3 = x*x*x - x*x;

    return Vector3Add(
        Vector3Add(Vector3Scale(p0, w0), Vector3Scale(p1, w1)),
        Vector3Add(Vector3Scale(v0, w2), Vector3Scale(v1, w3)));
}

static inline Vector3 Vector3InterpolateCubic(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, float alpha)
{
    Vector3 v1 = Vector3Scale(Vector3Add(Vector3Subtract(p1, p0), Vector3Subtract(p2, p1)), 0.5f);
    Vector3 v2 = Vector3Scale(Vector3Add(Vector3Subtract(p2, p1), Vector3Subtract(p3, p2)), 0.5f);
    return Vector3Hermite(p1, p2, v1, v2, alpha);
}

static inline Quaternion QuaternionHermite(Quaternion r0, Quaternion r1, Vector3 v0, Vector3 v1, float alpha)
{
    float x = alpha;
    float w1 = 3*x*x - 2*x*x*x;
    float w2 = x*x*x - 2*x*x + x;
    float w3 = x*x*x - x*x;

    Vector3 r1r0 = QuaternionToScaledAngleAxis(QuaternionAbsolute(QuaternionMultiply(r1, QuaternionInvert(r0))));

    return QuaternionMultiply(QuaternionFromScaledAngleAxis(
        Vector3Add(Vector3Add(Vector3Scale(r1r0, w1), Vector3Scale(v0, w2)), Vector3Scale(v1, w3))), r0);
}

static inline Quaternion QuaternionInterpolateCubic(Quaternion r0, Quaternion r1, Quaternion r2, Quaternion r3, float alpha)
{
    Vector3 r1r0 = QuaternionToScaledAngleAxis(QuaternionAbsolute(QuaternionMultiply(r1, QuaternionInvert(r0))));
    Vector3 r2r1 = QuaternionToScaledAngleAxis(QuaternionAbsolute(QuaternionMultiply(r2, QuaternionInvert(r1))));
    Vector3 r3r2 = QuaternionToScaledAngleAxis(QuaternionAbsolute(QuaternionMultiply(r3, QuaternionInvert(r2))));

    Vector3 v1 = Vector3Scale(Vector3Add(r1r0, r2r1), 0.5f);
    Vector3 v2 = Vector3Scale(Vector3Add(r2r1, r3r2), 0.5f);
    return QuaternionHermite(r1, r2, v1, v2, alpha);
}

// Frustum culling (based off https://github.com/JeffM2501/raylibExtras)

typedef struct
{
    Vector4 back;
    Vector4 front;
    Vector4 bottom;
    Vector4 top;
    Vector4 right;
    Vector4 left;
    
} Frustum;

static inline Vector4 FrustumPlaneNormalize(Vector4 plane)
{
    float magnitude = sqrtf(Square(plane.x) + Square(plane.y) + Square(plane.z));
    plane.x /= magnitude;
    plane.y /= magnitude;
    plane.z /= magnitude;
    plane.w /= magnitude;
    return plane;
}

static inline Frustum FrustumFromCameraMatrices(Matrix projection, Matrix modelview)
{
    Matrix planes = { 0 };
    planes.m0 = modelview.m0 * projection.m0 + modelview.m1 * projection.m4 + modelview.m2 * projection.m8 + modelview.m3 * projection.m12;
    planes.m1 = modelview.m0 * projection.m1 + modelview.m1 * projection.m5 + modelview.m2 * projection.m9 + modelview.m3 * projection.m13;
    planes.m2 = modelview.m0 * projection.m2 + modelview.m1 * projection.m6 + modelview.m2 * projection.m10 + modelview.m3 * projection.m14;
    planes.m3 = modelview.m0 * projection.m3 + modelview.m1 * projection.m7 + modelview.m2 * projection.m11 + modelview.m3 * projection.m15;
    planes.m4 = modelview.m4 * projection.m0 + modelview.m5 * projection.m4 + modelview.m6 * projection.m8 + modelview.m7 * projection.m12;
    planes.m5 = modelview.m4 * projection.m1 + modelview.m5 * projection.m5 + modelview.m6 * projection.m9 + modelview.m7 * projection.m13;
    planes.m6 = modelview.m4 * projection.m2 + modelview.m5 * projection.m6 + modelview.m6 * projection.m10 + modelview.m7 * projection.m14;
    planes.m7 = modelview.m4 * projection.m3 + modelview.m5 * projection.m7 + modelview.m6 * projection.m11 + modelview.m7 * projection.m15;
    planes.m8 = modelview.m8 * projection.m0 + modelview.m9 * projection.m4 + modelview.m10 * projection.m8 + modelview.m11 * projection.m12;
    planes.m9 = modelview.m8 * projection.m1 + modelview.m9 * projection.m5 + modelview.m10 * projection.m9 + modelview.m11 * projection.m13;
    planes.m10 = modelview.m8 * projection.m2 + modelview.m9 * projection.m6 + modelview.m10 * projection.m10 + modelview.m11 * projection.m14;
    planes.m11 = modelview.m8 * projection.m3 + modelview.m9 * projection.m7 + modelview.m10 * projection.m11 + modelview.m11 * projection.m15;
    planes.m12 = modelview.m12 * projection.m0 + modelview.m13 * projection.m4 + modelview.m14 * projection.m8 + modelview.m15 * projection.m12;
    planes.m13 = modelview.m12 * projection.m1 + modelview.m13 * projection.m5 + modelview.m14 * projection.m9 + modelview.m15 * projection.m13;
    planes.m14 = modelview.m12 * projection.m2 + modelview.m13 * projection.m6 + modelview.m14 * projection.m10 + modelview.m15 * projection.m14;
    planes.m15 = modelview.m12 * projection.m3 + modelview.m13 * projection.m7 + modelview.m14 * projection.m11 + modelview.m15 * projection.m15;

    Frustum frustum;
    frustum.back = FrustumPlaneNormalize((Vector4){ planes.m3 - planes.m2, planes.m7 - planes.m6, planes.m11 - planes.m10, planes.m15 - planes.m14 });
    frustum.front = FrustumPlaneNormalize((Vector4){ planes.m3 + planes.m2, planes.m7 + planes.m6, planes.m11 + planes.m10, planes.m15 + planes.m14 });
    frustum.bottom = FrustumPlaneNormalize((Vector4){ planes.m3 + planes.m1, planes.m7 + planes.m5, planes.m11 + planes.m9, planes.m15 + planes.m13 });
    frustum.top = FrustumPlaneNormalize((Vector4){ planes.m3 - planes.m1, planes.m7 - planes.m5, planes.m11 - planes.m9, planes.m15 - planes.m13 });
    frustum.left = FrustumPlaneNormalize((Vector4){ planes.m3 + planes.m0, planes.m7 + planes.m4, planes.m11 + planes.m8, planes.m15 + planes.m12 });
    frustum.right = FrustumPlaneNormalize((Vector4){ planes.m3 - planes.m0, planes.m7 - planes.m4, planes.m11 - planes.m8, planes.m15 - planes.m12 });
    return frustum;
}

static inline float FrustumPlaneDistanceToPoint(Vector4 plane, Vector3 position)
{
    return (plane.x * position.x + plane.y * position.y + plane.z * position.z + plane.w);
}

static inline bool FrustumContainsSphere(Frustum frustum, Vector3 position, float radius)
{
    if (FrustumPlaneDistanceToPoint(frustum.back, position) < -radius) { return false; }
    if (FrustumPlaneDistanceToPoint(frustum.front, position) < -radius) { return false; }
    if (FrustumPlaneDistanceToPoint(frustum.bottom, position) < -radius) { return false; }
    if (FrustumPlaneDistanceToPoint(frustum.top, position) < -radius) { return false; }
    if (FrustumPlaneDistanceToPoint(frustum.left, position) < -radius) { return false; }
    if (FrustumPlaneDistanceToPoint(frustum.right, position) < -radius) { return false; }
    return true;
}

//----------------------------------------------------------------------------------
// Command Line Args
//----------------------------------------------------------------------------------

// Finds an argument on the command line with the given name (in the format "--argName=argValue") and returns the argValue as a string
static inline char* ArgFind(int argc, char** argv, const char* name)
{
    for (int i = 1; i < argc; i++)
    {
        if (strlen(argv[i]) > 4 &&
          argv[i][0] == '-' &&
          argv[i][1] == '-' &&
          strstr(argv[i] + 2, name) == argv[i] + 2)
        {
            char* argStart = strchr(argv[i], '=');
            return argStart ? argStart + 1 : NULL;
        }
    }

    return NULL;
}

// Parse a float argument from the command line
static inline float ArgFloat(int argc, char** argv, const char* name, float defaultValue)
{
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }

    errno = 0;
    float output = strtof(value, NULL);
    if (errno == 0) { printf("INFO: Parsed option '%s' as '%s'\n", name, value); return output; }

    printf("ERROR: Could not parse value '%s' given for option '%s' as float\n", value, name);
    return defaultValue;
}

// Parse an integer argument from the command line
static inline int ArgInt(int argc, char** argv, const char* name, int defaultValue)
{
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }

    errno = 0;
    int output = (int)strtol(value, NULL, 10);
    if (errno == 0) { printf("INFO: Parsed option '%s' as '%s'\n", name, value); return output; }

    printf("ERROR: Could not parse value '%s' given for option '%s' as int\n", value, name);
    return defaultValue;
}

// Parse a boolean argument from the command line
static inline int ArgBool(int argc, char** argv, const char* name, bool defaultValue)
{
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }
    if (strcmp(value, "true") == 0) { printf("INFO: Parsed option '%s' as '%s'\n", name, value); return true; }
    if (strcmp(value, "false") == 0) { printf("INFO: Parsed option '%s' as '%s'\n", name, value); return false; }

    printf("ERROR: Could not parse value '%s' given for option '%s' as bool\n", value, name);
    return defaultValue;
}

// Parse an enum argument from the command line
static inline int ArgEnum(int argc, char** argv, const char* name, int optionCount, const char* options[], int defaultValue)
{
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }

    for (int i = 0; i < optionCount; i++)
    {
        if (strcmp(value, options[i]) == 0)
        {
            printf("INFO: Parsed option '%s' as '%s'\n", name, value);
            return i;
        }
    }

    printf("ERROR: Could not parse value '%s' given for option '%s' as enum\n", value, name);
    return defaultValue;
}

// Parse a string argument from the command line
static inline const char* ArgStr(int argc, char** argv, const char* name, const char* defaultValue)
{
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }

    printf("INFO: Parsed option '%s' as '%s'\n", name, value);
    return value;
}

// Parse a color argument from the command line
static inline Color ArgColor(int argc, char** argv, const char* name, Color defaultValue)
{
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }

    int cx, cy, cz;
    if (sscanf(value, "%i,%i,%i", &cx, &cy, &cz) == 3)
    {
        printf("INFO: Parsed option '%s' as '%s'\n", name, value);
        return (Color){ ClampInt(cx, 0, 255), ClampInt(cy, 0, 255), ClampInt(cz, 0, 255) };
    }

    printf("ERROR: Could not parse value '%s' given for option '%s' as color\n", value, name);
    return defaultValue;
}

// Parse a vector3 argument from the command line
static inline Vector3 ArgVector3(int argc, char** argv, const char* name, Vector3 defaultValue)
{
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }

    float cx, cy, cz;
    if (sscanf(value, "%f,%f,%f", &cx, &cy, &cz) == 3)
    {
        printf("INFO: Parsed option '%s' as '%s'\n", name, value);
        return (Vector3){ cx, cy, cz };
    }

    printf("ERROR: Could not parse value '%s' given for option '%s' as color\n", value, name);
    return defaultValue;
}

//----------------------------------------------------------------------------------
// Camera
//----------------------------------------------------------------------------------

// Basic Orbit Camera with simple controls
typedef struct {

    Camera3D cam3d;
    float azimuth;
    float altitude;
    float distance;
    Vector3 offset;
    bool track;
    int trackBone;

} OrbitCamera;

static inline void OrbitCameraInit(OrbitCamera* camera, int argc, char** argv)
{
    memset(&camera->cam3d, 0, sizeof(Camera3D));
    camera->cam3d.position = (Vector3){ 2.0f, 3.0f, 5.0f };
    camera->cam3d.target = (Vector3){ -0.5f, 1.0f, 0.0f };
    camera->cam3d.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera->cam3d.fovy = ArgFloat(argc, argv, "cameraFOV", 45.0f);
    camera->cam3d.projection = CAMERA_PERSPECTIVE;

    camera->azimuth = ArgFloat(argc, argv, "cameraAzimuth", 0.0f);
    camera->altitude = ArgFloat(argc, argv, "cameraAltitude", 0.4f);
    camera->distance = ArgFloat(argc, argv, "cameraDistance", 4.0f);
    camera->offset = ArgVector3(argc, argv, "cameraOffset", Vector3Zero());
    camera->track = ArgBool(argc, argv, "cameraTrack", false);
    camera->trackBone = ArgInt(argc, argv, "cameraTrackBone", 0);
}

static inline void OrbitCameraUpdate(
    OrbitCamera* camera,
    Vector3 target,
    float azimuthDelta,
    float altitudeDelta,
    float offsetDeltaX,
    float offsetDeltaY,
    float mouseWheel,
    float dt)
{
    camera->azimuth = camera->azimuth + 1.0f * dt * -azimuthDelta;
    camera->altitude = Clamp(camera->altitude + 1.0f * dt * altitudeDelta, 0.0, 0.4f * PI);
    camera->distance = Clamp(camera->distance +  20.0f * dt * -mouseWheel, 0.1f, 100.0f);
    
    Quaternion rotationAzimuth = QuaternionFromAxisAngle((Vector3){0, 1, 0}, camera->azimuth);
    Vector3 position = Vector3RotateByQuaternion((Vector3){0, 0, camera->distance}, rotationAzimuth);
    Vector3 axis = Vector3Normalize(Vector3CrossProduct(position, (Vector3){0, 1, 0}));

    Quaternion rotationAltitude = QuaternionFromAxisAngle(axis, camera->altitude);

    Vector3 localOffset = (Vector3){ dt * offsetDeltaX, dt * -offsetDeltaY, 0.0f };
    localOffset = Vector3RotateByQuaternion(localOffset, rotationAzimuth);

    camera->offset = Vector3Add(camera->offset, Vector3RotateByQuaternion(localOffset, rotationAltitude));

    Vector3 cameraTarget = Vector3Add(camera->offset, target);
    Vector3 eye = Vector3Add(cameraTarget, Vector3RotateByQuaternion(position, rotationAltitude));

    camera->cam3d.target = cameraTarget;
    camera->cam3d.position = eye;
}

//----------------------------------------------------------------------------------
// Parser
//----------------------------------------------------------------------------------

enum
{
    PARSER_ERR_MAX = 512,
};

// Simple parser that keeps track of rows, and cols in a string and so can provide slightly 
// nicer error messages. Has ability to peek at next character and advance the input
typedef struct {

    const char* filename;
    int offset;
    const char* data;
    int row;
    int col;
    char err[PARSER_ERR_MAX];

} Parser;

// Initialize the Parser
static inline void ParserInit(Parser* par, const char* filename, const char* data)
{
    par->filename = filename;
    par->offset = 0;
    par->data = data;
    par->row = 0;
    par->col = 0;
    par->err[0] = '\0';
}

// Peek at the next character in the stream
static inline char ParserPeek(const Parser* par)
{
    return par->data[par->offset];
}

// Peek forward N steps in the stream. Does not check the stream is long enough.
static inline char ParserPeekForward(const Parser* par, int steps)
{
    return par->data[par->offset + steps];
}

// Checks the current character matches the given input
static inline bool ParserMatch(const Parser* par, char match)
{
    return match == par->data[par->offset];
}

// Checks the current character matches one of the given characters
static inline bool ParserOneOf(const Parser* par, const char* matches)
{
    return strchr(matches, par->data[par->offset]);
}

// Checks the following characters in the stream match the prefix (in a caseless way)
static inline bool ParserStartsWithCaseless(const Parser* par, const char* prefix)
{
    const char* start = par->data + par->offset;
    while (*prefix)
    {
        if (tolower(*prefix) != tolower(*start)) { return false; }
        prefix++;
        start++;
    }

    return true;
}

// Advances the stream forward one
static inline void ParserInc(Parser* par)
{
    if (par->data[par->offset] == '\n')
    {
        par->row++;
        par->col = 0;
    }
    else
    {
        par->col++;
    }

    par->offset++;
}

// Advances the stream forward "num" characters
static inline void ParserAdvance(Parser* par, int num)
{
    for (int i = 0; i < num; i++) { ParserInc(par); }
}

// Gets the human readable name of a particular character
static inline char* ParserCharName(char c)
{
    static char parserCharName[2];

    switch (c)
    {
        case '\0': return "end of file";
        case '\r': return "new line";
        case '\n': return "new line";
        case '\t': return "tab";
        case '\v': return "vertical tab";
        case '\b': return "backspace";
        case '\f': return "form feed";
        default:
            parserCharName[0] = c;
            parserCharName[1] = '\0';
            return parserCharName;
    }
}

// Prints a formatted error to the parser error buffer
#define ParserError(par, fmt, ...) \
    snprintf(par->err, PARSER_ERR_MAX, "%s:%i:%i: error: " fmt, par->filename, par->row, par->col, ##__VA_ARGS__)

//----------------------------------------------------------------------------------
// BVH File Data
//----------------------------------------------------------------------------------

// Types of "channels" that are possible in the BVH format
enum
{
    CHANNEL_X_POSITION = 0,
    CHANNEL_Y_POSITION = 1,
    CHANNEL_Z_POSITION = 2,
    CHANNEL_X_ROTATION = 3,
    CHANNEL_Y_ROTATION = 4,
    CHANNEL_Z_ROTATION = 5,
    CHANNELS_MAX = 6,
};

// Data associated with a single "joint" in the BVH format
typedef struct
{
    int parent;
    char* name;
    Vector3 offset;
    int channelCount;
    char channels[CHANNELS_MAX];
    bool endSite;

} BVHJointData;

static inline void BVHJointDataInit(BVHJointData* data)
{
    data->parent = -1;
    data->name = NULL;
    data->offset = (Vector3){ 0.0f, 0.0f, 0.0f };
    data->channelCount = 0;
    data->endSite = false;
}

static inline void BVHJointDataRename(BVHJointData* data, const char* name)
{
    data->name = realloc(data->name, strlen(name) + 1);
    strcpy(data->name, name);
}

static inline void BVHJointDataFree(BVHJointData* data)
{
    free(data->name);
}

// Data structure matching what is present in the BVH file format
typedef struct
{
    // Hierarchy Data

    int jointCount;
    BVHJointData* joints;

    // Motion Data

    int frameCount;
    int channelCount;
    float frameTime;
    float* motionData;

} BVHData;

static inline void BVHDataInit(BVHData* bvh)
{
    bvh->jointCount = 0;
    bvh->joints = NULL;
    bvh->frameCount = 0;
    bvh->channelCount = 0;
    bvh->frameTime = 0.0f;
    bvh->motionData = NULL;
}

static inline void BVHDataFree(BVHData* bvh)
{
    for (int i = 0; i < bvh->jointCount; i++)
    {
        BVHJointDataFree(&bvh->joints[i]);
    }
    free(bvh->joints);

    free(bvh->motionData);
}

static inline int BVHDataAddJoint(BVHData* bvh)
{
    bvh->joints = (BVHJointData*)realloc(bvh->joints, (bvh->jointCount + 1) * sizeof(BVHJointData));
    bvh->jointCount++;
    BVHJointDataInit(&bvh->joints[bvh->jointCount - 1]);
    return bvh->jointCount - 1;
}

//----------------------------------------------------------------------------------
// BVHViewer Additions by Teodor Nikolov
//----------------------------------------------------------------------------------
static inline void CharacterModelInit(CharacterModel* model)
{
    const char* modelPath = "./assets/GENEA_Model.gltf";
    BVHALoadCharacterModelFromFile(model, modelPath);
}

static inline void CharacterModelFree(CharacterModel* model)
{
    BVHAUnloadCharacterModel(model);
}

static BoneInfo convertBVHJointToBoneInfo(BVHJointData* joint)
{
    BoneInfo bone;
    strncpy(bone.name, joint->name, sizeof(bone.name));
    bone.parent = joint->parent;
    return bone;
}

static char** getCommonBoneNames(const BVHData* bvhData, const Model* model, int* commonBoneCount)
{
    char** commonBoneNames = malloc(sizeof(char*) * bvhData->jointCount);
    *commonBoneCount = 0;

    for (int i = 0; i < bvhData->jointCount; ++i) {
        if (boneExistsInModel(bvhData->joints[i].name, model)) {
            commonBoneNames[*commonBoneCount] = malloc(sizeof(char) * (strlen(bvhData->joints[i].name) + 1));
            strcpy(commonBoneNames[*commonBoneCount], bvhData->joints[i].name);
            (*commonBoneCount)++;
        }
    }

    return commonBoneNames;
}

static ModelAnimation BVHToModelAnimation(const BVHData* bvhData, const Model* model)
{
    // int commonBoneCount = 0;
    // char** commonBoneNames = getCommonBoneNames(bvhData, model, &commonBoneCount);

    ModelAnimation animation;
    animation.boneCount = model->boneCount;
    animation.frameCount = bvhData->frameCount;
    animation.bones = model->bones;

    // Convert BVHData animation frames to ModelAnimation animation frames
    animation.framePoses = malloc(sizeof(Transform*) * animation.frameCount);
    for (int i = 0; i < animation.frameCount; ++i)
    {
        int boneCount = animation.boneCount;
        Matrix* boneTransformMats = malloc(sizeof(Matrix) * boneCount);
        animation.framePoses[i] = malloc(sizeof(Transform) * boneCount);

        int channelOffset = 0;
        for (int j = 0; j < boneCount; ++j)
        {
            int dataIndex = i * bvhData->channelCount + channelOffset;
            int boneId = j;
            int boneIdParent = animation.bones[j].parent;
            assert(boneIdParent < boneId);

            Vector3 boneOffset;
            boneOffset.x = bvhData->motionData[dataIndex+0] / 100.0f;
            boneOffset.y = bvhData->motionData[dataIndex+1] / 100.0f;
            boneOffset.z = bvhData->motionData[dataIndex+2] / 100.0f;

            Vector3 rotation; // Euler angles
            rotation.x = bvhData->motionData[dataIndex+4] * DEG2RAD;
            rotation.y = bvhData->motionData[dataIndex+5] * DEG2RAD;
            rotation.z = bvhData->motionData[dataIndex+3] * DEG2RAD;
            channelOffset += bvhData->joints[j].channelCount;

            // Calculate local transform matrix
            Matrix matTranslation = MatrixTranslate(boneOffset.x, boneOffset.y, boneOffset.z); // Bones do not slide
            Matrix matRotationX = MatrixRotateX(rotation.x);
            Matrix matRotationY = MatrixRotateY(rotation.y);
            Matrix matRotationZ = MatrixRotateZ(rotation.z);
            Matrix matRotation = MatrixMultiply(MatrixMultiply(matRotationY, matRotationX), matRotationZ);
            Matrix matScale = MatrixScale(1.0f, 1.0f, 1.0f); // Bones do not stretch
            Matrix matTransform = MatrixMultiply(MatrixMultiply(matScale, matRotation), matTranslation);

            // Compute model-level bone transform matrix from parent bone
            Matrix matTransformParent;
            if (j == 0) {
                matTransformParent = MatrixIdentity();
            }
            else {
                matTransformParent = boneTransformMats[boneIdParent];
            }
            Matrix matTransformBone = MatrixMultiply(matTransform, matTransformParent);
            boneTransformMats[j] = matTransformBone;

            // Set the final transformation for this bone in this frame
            Transform transform;
            transform.translation = Vector3Transform((Vector3){0.0f, 0.0f, 0.0f}, matTransformBone);
            transform.rotation = QuaternionFromMatrix(matTransformBone);
            transform.scale = (Vector3){1.0f, 1.0f, 1.0f};
            animation.framePoses[i][j] = transform;
        }

        free(boneTransformMats);
    }

    strcpy(animation.name, "placeholder_text");

    return animation;
}

//----------------------------------------------------------------------------------
// BVH Parser
//----------------------------------------------------------------------------------

// Parse any whitespace
static void BVHParseWhitespace(Parser* par)
{
    while (ParserOneOf(par, " \r\t\v")) { ParserInc(par); }
}

// Parse the given string (in a non-case sensitive way). I've found that in practice
// many BVH files don't respect case sensitivity so parsing any keywords in a non-case
// sensitive way seems safer.
static bool BVHParseString(Parser* par, const char* string)
{
    if (ParserStartsWithCaseless(par, string))
    {
        ParserAdvance(par, strlen(string));
        return true;
    }
    else
    {
        ParserError(par, "expected '%s' at '%s'", string, ParserCharName(ParserPeek(par)));
        return false;
    }
}

// Parse any whitespace followed by a newline
static bool BVHParseNewline(Parser* par)
{
    BVHParseWhitespace(par);

    if (ParserMatch(par, '\n'))
    {
        ParserInc(par);
        BVHParseWhitespace(par);
        return true;
    }
    else
    {
        ParserError(par, "expected newline at '%s'", ParserCharName(ParserPeek(par)));
        return false;
    }
}

// Parse any whitespace and then an identifier for the name of a joint
static bool BVHParseJointName(BVHJointData* jnt, Parser* par)
{
    BVHParseWhitespace(par);

    char buffer[256];
    int chrnum = 0;
    while (chrnum < 255 && ParserOneOf(par,
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_:-"))
    {
        buffer[chrnum] = ParserPeek(par);
        chrnum++;
        ParserInc(par);
    }
    buffer[chrnum] = '\0';

    if (chrnum > 0)
    {
        BVHJointDataRename(jnt, buffer);
        BVHParseWhitespace(par);
        return true;
    }
    else
    {
        ParserError(par, "expected joint name at '%s'", ParserCharName(ParserPeek(par)));
        return false;
    }
}

// Parse a float value
static bool BVHParseFloat(float* out, Parser* par)
{
    BVHParseWhitespace(par);

    char* end;
    errno = 0;
    (*out) = strtof(par->data + par->offset, &end);

    if (errno == 0)
    {
        ParserAdvance(par, end - (par->data + par->offset));
        return true;
    }
    else
    {
        ParserError(par, "expected float at '%s'", ParserCharName(ParserPeek(par)));
        return false;
    }
}

// Parse an integer value
static bool BVHParseInt(int* out, Parser* par)
{
    BVHParseWhitespace(par);

    char* end;
    errno = 0;
    (*out) = (int)strtol(par->data + par->offset, &end, 10);

    if (errno == 0)
    {
        ParserAdvance(par, end - (par->data + par->offset));
        return true;
    }
    else
    {
        ParserError(par, "expected integer at '%s'", ParserCharName(ParserPeek(par)));
        return false;
    }
}

// Parse the "joint offset" part of the BVH File
static bool BVHParseJointOffset(BVHJointData* jnt, Parser* par)
{
    if (!BVHParseString(par, "OFFSET")) { return false; }
    if (!BVHParseFloat(&jnt->offset.x, par)) { return false; }
    if (!BVHParseFloat(&jnt->offset.y, par)) { return false; }
    if (!BVHParseFloat(&jnt->offset.z, par)) { return false; }
    if (!BVHParseNewline(par)) { return false; }
    return true;
}

// Parse a channel type and return it in "channel"
static bool BVHParseChannelEnum(
    char* channel,
    Parser* par,
    const char* channelName,
    char channelValue)
{
    BVHParseWhitespace(par);
    if (!BVHParseString(par, channelName)) { return false; }
    BVHParseWhitespace(par);
    *channel = channelValue;
    return true;
}

// Parse a channel type and return it in "channel"
static bool BVHParseChannel(char* channel, Parser* par)
{
    BVHParseWhitespace(par);

    if (ParserPeek(par) == '\0')
    {
        ParserError(par, "expected channel at end of file");
        return false;
    }

    // Here we are safe to peek forward an extra character since we've already
    // checked the current character is not the null terminator.

    if (ParserPeek(par) == 'X' && ParserPeekForward(par, 1) == 'p')
    {
        return BVHParseChannelEnum(channel, par, "Xposition", CHANNEL_X_POSITION);
    }

    if (ParserPeek(par) == 'Y' && ParserPeekForward(par, 1) == 'p')
    {
        return BVHParseChannelEnum(channel, par, "Yposition", CHANNEL_Y_POSITION);
    }

    if (ParserPeek(par) == 'Z' && ParserPeekForward(par, 1) == 'p')
    {
        return BVHParseChannelEnum(channel, par, "Zposition", CHANNEL_Z_POSITION);
    }

    if (ParserPeek(par) == 'X' && ParserPeekForward(par, 1) == 'r')
    {
        return BVHParseChannelEnum(channel, par, "Xrotation", CHANNEL_X_ROTATION);
    }

    if (ParserPeek(par) == 'Y' && ParserPeekForward(par, 1) == 'r')
    {
        return BVHParseChannelEnum(channel, par, "Yrotation", CHANNEL_Y_ROTATION);
    }

    if (ParserPeek(par) == 'Z' && ParserPeekForward(par, 1) == 'r')
    {
        return BVHParseChannelEnum(channel, par, "Zrotation", CHANNEL_Z_ROTATION);
    }

    ParserError(par, "expected channel type");
    return false;
}

// Parse the "channels" part of the BVH file format
static bool BVHParseJointChannels(BVHJointData* jnt, Parser* par)
{
    if (!BVHParseString(par, "CHANNELS")) { return false; }
    if (!BVHParseInt(&jnt->channelCount, par)) { return false; }

    for (int i = 0; i < jnt->channelCount; i++)
    {
        if (!BVHParseChannel(&jnt->channels[i], par)) { return false; }
    }

    if (!BVHParseNewline(par)) { return false; }

    return true;
}

// Parse a joint in the BVH file format
static bool BVHParseJoints(BVHData* bvh, int parent, Parser* par)
{
    while (ParserOneOf(par, "JEje")) // Either "JOINT" or "End Site"
    {
        int j = BVHDataAddJoint(bvh);
        bvh->joints[j].parent = parent;

        if (ParserMatch(par, 'J'))
        {
            if (!BVHParseString(par, "JOINT")) { return false; }
            if (!BVHParseJointName(&bvh->joints[j], par)) { return false; }
            if (!BVHParseNewline(par)) { return false; }
            if (!BVHParseString(par, "{")) { return false; }
            if (!BVHParseNewline(par)) { return false; }
            if (!BVHParseJointOffset(&bvh->joints[j], par)) { return false; }
            if (!BVHParseJointChannels(&bvh->joints[j], par)) { return false; }
            if (!BVHParseJoints(bvh, j, par)) { return false; }
            if (!BVHParseString(par, "}")) { return false; }
            if (!BVHParseNewline(par)) { return false; }
        }
        else if (ParserMatch(par, 'E'))
        {
            bvh->joints[j].endSite = true;

            if (!BVHParseString(par, "End Site")) { return false; }
            BVHJointDataRename(&bvh->joints[j], "End Site");
            if (!BVHParseNewline(par)) { return false; }
            if (!BVHParseString(par, "{")) { return false; }
            if (!BVHParseNewline(par)) { return false; }
            if (!BVHParseJointOffset(&bvh->joints[j], par)) { return false; }
            if (!BVHParseString(par, "}")) { return false; }
            if (!BVHParseNewline(par)) { return false; }
        }
    }

    return true;
}

// Parse the frame count
static bool BVHParseFrames(BVHData* bvh, Parser* par)
{
    if (!BVHParseString(par, "Frames:")) { return false; }
    if (!BVHParseInt(&bvh->frameCount, par)) { return false; }
    if (!BVHParseNewline(par)) { return false; }
    return true;
}

// Parse the frame time
static bool BVHParseFrameTime(BVHData* bvh, Parser* par)
{
    if (!BVHParseString(par, "Frame Time:")) { return false; }
    if (!BVHParseFloat(&bvh->frameTime, par)) { return false; }
    if (!BVHParseNewline(par)) { return false; }
    if (bvh->frameTime == 0.0f) { bvh->frameTime = 1.0f / 60.0f; }
    return true;
}

// Parse the motion data part of the BVH file format
static bool BVHParseMotionData(BVHData* bvh, Parser* par)
{
    int channelCount = 0;
    for (int i = 0; i < bvh->jointCount; i++)
    {
        channelCount += bvh->joints[i].channelCount;
    }

    bvh->channelCount = channelCount;
    bvh->motionData = malloc(bvh->frameCount * channelCount * sizeof(float));

    for (int i = 0; i < bvh->frameCount; i++)
    {
        for (int j = 0; j < channelCount; j++)
        {
            if (!BVHParseFloat(&bvh->motionData[i * channelCount + j], par)) { return false; }
        }

        if (!BVHParseNewline(par)) { return false; }
    }

    return true;
}

// Parse the entire BVH file format
static bool BVHParse(BVHData* bvh, Parser* par)
{
    // Hierarchy Data

    if (!BVHParseString(par, "HIERARCHY")) { return false; }
    if (!BVHParseNewline(par)) { return false; }

    int j = BVHDataAddJoint(bvh);

    if (!BVHParseString(par, "ROOT")) { return false; }
    if (!BVHParseJointName(&bvh->joints[j], par)) { return false; }
    if (!BVHParseNewline(par)) { return false; }
    if (!BVHParseString(par, "{")) { return false; }
    if (!BVHParseNewline(par)) { return false; }
    if (!BVHParseJointOffset(&bvh->joints[j], par)) { return false; }
    if (!BVHParseJointChannels(&bvh->joints[j], par)) { return false; }
    if (!BVHParseJoints(bvh, j, par)) { return false; }
    if (!BVHParseString(par, "}")) { return false; }
    if (!BVHParseNewline(par)) { return false; }

    // Motion Data

    if (!BVHParseString(par, "MOTION")) { return false; }
    if (!BVHParseNewline(par)) { return false; }
    if (!BVHParseFrames(bvh, par)) { return false; }
    if (!BVHParseFrameTime(bvh, par)) { return false; }
    if (!BVHParseMotionData(bvh, par)) { return false; }

    return true;
}

// Load the given file and parse the contents as a BVH file.
static bool BVHDataLoad(BVHData* bvh, const char* filename, char* errMsg, int errMsgSize)
{
    // Read file Contents

    FILE* f = fopen(filename, "rb");

    if (f == NULL)
    {
        snprintf(errMsg, errMsgSize, "Error: Could not find file '%s'\n", filename);
        return false;
    }

    fseek(f, 0, SEEK_END);
    long int length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buffer = malloc(length + 2);
    fread(buffer, 1, length, f);
    buffer[length] = '\n';
    buffer[length+1] = '\n'; // Prevents segmentation fault during parsing of .bvh file that does not end with an empty line
    fclose(f);
    
    // Free and re-init in case we are re-using an old buffer 
    BVHDataFree(bvh); 
    BVHDataInit(bvh);

    // Parse BVH
    Parser par;
    ParserInit(&par, filename, buffer);
    bool result = BVHParse(bvh, &par);

    // Free contents and return result
    free(buffer);

    if (!result)
    {
        snprintf(errMsg, errMsgSize, "Error: Could not parse BVH file:\n    %s", par.err);
    }
    else
    {
        errMsg[0] = '\0';
        printf("INFO: parsed '%s' successfully\n", filename);
    }

    return result;
}

//----------------------------------------------------------------------------------
// Transform Data
//----------------------------------------------------------------------------------

// Structure for containing a sampled pose as joint transforms
typedef struct
{
    int jointCount;
    int* parents;
    bool* endSite;
    Vector3* localPositions;
    Quaternion* localRotations;
    Vector3* globalPositions;
    Quaternion* globalRotations;

} TransformData;

static inline void TransformDataInit(TransformData* data)
{
    data->jointCount = 0;
    data->parents = NULL;
    data->endSite = NULL;
    data->localPositions = NULL;
    data->localRotations = NULL;
    data->globalPositions = NULL;
    data->globalRotations = NULL;
}

// Resize the transform buffer according to the given BVH data and record the joint
// parents and end-sites.
static inline void TransformDataResize(TransformData* data, BVHData* bvh)
{
    data->jointCount = bvh->jointCount;
    data->parents = realloc(data->parents, data->jointCount * sizeof(int));
    data->endSite = realloc(data->endSite, data->jointCount * sizeof(bool));
    data->localPositions = realloc(data->localPositions, data->jointCount * sizeof(Vector3));
    data->localRotations = realloc(data->localRotations, data->jointCount * sizeof(Quaternion));
    data->globalPositions = realloc(data->globalPositions, data->jointCount * sizeof(Vector3));
    data->globalRotations = realloc(data->globalRotations, data->jointCount * sizeof(Quaternion));

    for (int i = 0; i < data->jointCount; i++)
    {
        data->endSite[i] = bvh->joints[i].endSite;
        data->parents[i] = bvh->joints[i].parent;
    }
}

static inline void TransformDataFree(TransformData* data)
{
    free(data->parents);
    free(data->endSite);
    free(data->localPositions);
    free(data->localRotations);
    free(data->globalPositions);
    free(data->globalRotations);
}

// Sample joint transforms from a given frame of the BVH file and with a given scale
static void TransformDataSampleFrame(TransformData* data, BVHData* bvh, int frame, float scale)
{
    // Clamp the frame index in range.
    frame = frame < 0 ? 0 : frame >= bvh->frameCount ? bvh->frameCount - 1 : frame;

    int offset = 0;
    for (int i = 0; i < bvh->jointCount; i++)
    {
        Vector3 position = Vector3Scale(bvh->joints[i].offset, scale);
        Quaternion rotation = QuaternionIdentity();

        for (int c = 0; c < bvh->joints[i].channelCount; c++)
        {
            switch (bvh->joints[i].channels[c])
            {
                case CHANNEL_X_POSITION:
                    position.x = scale * bvh->motionData[frame * bvh->channelCount + offset];
                    offset++;
                    break;

                case CHANNEL_Y_POSITION:
                    position.y = scale * bvh->motionData[frame * bvh->channelCount + offset];
                    offset++;
                    break;

                case CHANNEL_Z_POSITION:
                    position.z = scale * bvh->motionData[frame * bvh->channelCount + offset];
                    offset++;
                    break;

                case CHANNEL_X_ROTATION:
                    rotation = QuaternionMultiply(rotation, QuaternionFromAxisAngle(
                        (Vector3){1, 0, 0}, DEG2RAD * bvh->motionData[frame * bvh->channelCount + offset]));
                    offset++;
                    break;

                case CHANNEL_Y_ROTATION:
                    rotation = QuaternionMultiply(rotation, QuaternionFromAxisAngle(
                        (Vector3){0, 1, 0}, DEG2RAD * bvh->motionData[frame * bvh->channelCount + offset]));
                    offset++;
                    break;

                case CHANNEL_Z_ROTATION:
                    rotation = QuaternionMultiply(rotation, QuaternionFromAxisAngle(
                        (Vector3){0, 0, 1}, DEG2RAD * bvh->motionData[frame * bvh->channelCount + offset]));
                    offset++;
                    break;
            }
        }

        data->localPositions[i] = position;
        data->localRotations[i] = rotation;
    }

    assert(offset == bvh->channelCount);
}

// Sample the nearest frame to the given time
static void TransformDataSampleFrameNearest(TransformData* data, BVHData* bvh, float time, float scale)
{
    int frame = ClampInt((int)(time / bvh->frameTime + 0.5f), 0, bvh->frameCount - 1);
    TransformDataSampleFrame(data, bvh, frame, scale);
}

// Perform a basic linear interpolation of the frame data in the BVH file
static void TransformDataSampleFrameLinear(
    TransformData* data,
    TransformData* tmp0,
    TransformData* tmp1,
    BVHData* bvh,
    float time,
    float scale)
{
    const float alpha = fmod(time / bvh->frameTime, 1.0f);
    int frame0 = ClampInt((int)(time / bvh->frameTime) + 0, 0, bvh->frameCount - 1);
    int frame1 = ClampInt((int)(time / bvh->frameTime) + 1, 0, bvh->frameCount - 1);

    TransformDataSampleFrame(tmp0, bvh, frame0, scale);
    TransformDataSampleFrame(tmp1, bvh, frame1, scale);

    for (int i = 0; i < data->jointCount; i++)
    {
        data->localPositions[i] = Vector3Lerp(tmp0->localPositions[i], tmp1->localPositions[i], alpha);
        data->localRotations[i] = QuaternionSlerp(tmp0->localRotations[i], tmp1->localRotations[i], alpha);
    }
}

// Perform a cubic interpolation of the frame data in the BVH file
static void TransformDataSampleFrameCubic(
    TransformData* data,
    TransformData* tmp0,
    TransformData* tmp1,
    TransformData* tmp2,
    TransformData* tmp3,
    BVHData* bvh,
    float time,
    float scale)
{
    const float alpha = fmod(time / bvh->frameTime, 1.0f);
    int frame0 = ClampInt((int)(time / bvh->frameTime) - 1, 0, bvh->frameCount - 1);
    int frame1 = ClampInt((int)(time / bvh->frameTime) + 0, 0, bvh->frameCount - 1);
    int frame2 = ClampInt((int)(time / bvh->frameTime) + 1, 0, bvh->frameCount - 1);
    int frame3 = ClampInt((int)(time / bvh->frameTime) + 2, 0, bvh->frameCount - 1);

    TransformDataSampleFrame(tmp0, bvh, frame0, scale);
    TransformDataSampleFrame(tmp1, bvh, frame1, scale);
    TransformDataSampleFrame(tmp2, bvh, frame2, scale);
    TransformDataSampleFrame(tmp3, bvh, frame3, scale);

    for (int i = 0; i < data->jointCount; i++)
    {
        data->localPositions[i] = Vector3InterpolateCubic(
            tmp0->localPositions[i], tmp1->localPositions[i],
            tmp2->localPositions[i], tmp3->localPositions[i], alpha);

        data->localRotations[i] = QuaternionInterpolateCubic(
            tmp0->localRotations[i], tmp1->localRotations[i],
            tmp2->localRotations[i], tmp3->localRotations[i], alpha);
    }
}

// Compute format kinematics on the transform buffer
static void TransformDataForwardKinematics(TransformData* data)
{
    for (int i = 0; i < data->jointCount; i++)
    {
        int p = data->parents[i];
        assert(p <= i);

        if (p == -1)
        {
            data->globalPositions[i] = data->localPositions[i];
            data->globalRotations[i] = data->localRotations[i];
        }
        else
        {
            data->globalPositions[i] = Vector3Add(Vector3RotateByQuaternion(data->localPositions[i], data->globalRotations[p]), data->globalPositions[p]);
            data->globalRotations[i] = QuaternionMultiply(data->globalRotations[p], data->localRotations[i]);
        }
    }
}

//----------------------------------------------------------------------------------
// Character Data
//----------------------------------------------------------------------------------

// Maximum number of characters to allow in the scene
enum
{
    CHARACTERS_MAX = 6,
};

// All the data required for all of the characters we want to have in the scene
typedef struct {

    // Total number of characters
    int count;
    
    // Character which is "active" or selected
    int active;

    // Character BVH Data
    BVHData bvhData[CHARACTERS_MAX];
    ModelAnimation animData[CHARACTERS_MAX];
    
    // Scales of each character
    float scales[CHARACTERS_MAX];
    
    // Names of each character
    char names[CHARACTERS_MAX][128];
    
    // Automatic scaling for each character
    float autoScales[CHARACTERS_MAX];
    
    // Color of each character
    Color colors[CHARACTERS_MAX];
    
    // Opacity of each character
    float opacities[CHARACTERS_MAX];
    
    // Maximum capsule radius of each character
    float radii[CHARACTERS_MAX];
    
    // Original file path for each character
    char filePaths[CHARACTERS_MAX][512];
    
    // Transform buffers for each character
    TransformData xformData[CHARACTERS_MAX];
    TransformData xformTmp0[CHARACTERS_MAX];
    TransformData xformTmp1[CHARACTERS_MAX];
    TransformData xformTmp2[CHARACTERS_MAX];
    TransformData xformTmp3[CHARACTERS_MAX];
    
    // Joint combo string for each character
    char* jointNamesCombo[CHARACTERS_MAX];

    // If the color picker is active
    bool colorPickerActive;

} CharacterData;

// Initializes all the CharacterData to a safe state
static inline void CharacterDataInit(CharacterData* data, int argc, char** argv)
{  
    data->count = 0;
    data->active = 0;

    data->colors[0] = ArgColor(argc, argv, "-color0", ORANGE);
    data->colors[1] = ArgColor(argc, argv, "-color1", (Color){ 38, 134, 157, 255 });
    data->colors[2] = ArgColor(argc, argv, "-color2", PINK);
    data->colors[3] = ArgColor(argc, argv, "-color3", LIME);
    data->colors[4] = ArgColor(argc, argv, "-color4", VIOLET);
    data->colors[5] = ArgColor(argc, argv, "-color5", MAROON);

    srand(1234);
    for (int i = 6; i < CHARACTERS_MAX; i++)
    {
        data->colors[i] = (Color){ rand() % 255, rand() % 255, rand() % 255  };
    }

    for (int i = 0; i < CHARACTERS_MAX; i++)
    {
        BVHDataInit(&data->bvhData[i]);
        data->scales[i] = 1.0f;
        data->names[i][0] = '\0';
        data->autoScales[i] = 1.0f;
        data->opacities[i] = ArgFloat(argc, argv, "capsuleOpacity", 1.0f);
        data->radii[i] = ArgFloat(argc, argv, "maxCapsuleRadius", 0.04f);
        data->filePaths[i][0] = '\0';
        TransformDataInit(&data->xformData[i]);
        TransformDataInit(&data->xformTmp0[i]);
        TransformDataInit(&data->xformTmp1[i]);
        TransformDataInit(&data->xformTmp2[i]);
        TransformDataInit(&data->xformTmp3[i]);
        data->jointNamesCombo[i] = NULL;
    }

    data->colorPickerActive = ArgBool(argc, argv, "colorPickerActive", false);
}

static inline void CharacterDataFree(CharacterData* data)
{
    for (int i = 0; i < data->count; i++)
    {
        TransformDataFree(&data->xformData[i]);
        TransformDataFree(&data->xformTmp0[i]);
        TransformDataFree(&data->xformTmp1[i]);
        TransformDataFree(&data->xformTmp2[i]);
        TransformDataFree(&data->xformTmp3[i]);
        BVHDataFree(&data->bvhData[i]);
        free(data->jointNamesCombo[i]);
    }
}

// Attempt to load a new character from the given file path
static bool CharacterDataLoadFromFile(
    CharacterData* data,
    const char* path,
    char* errMsg,
    int errMsgSize)
{
    printf("INFO: Loading '%s'\n", path);

    if (data->count == CHARACTERS_MAX)
    {
        snprintf(errMsg, 512, "Error: Maximum number of BVH files loaded (%i)", CHARACTERS_MAX);
        return false;
    }

    if (BVHDataLoad(&data->bvhData[data->count], path, errMsg, errMsgSize))
    {
        TransformDataResize(&data->xformData[data->count], &data->bvhData[data->count]);
        TransformDataResize(&data->xformTmp0[data->count], &data->bvhData[data->count]);
        TransformDataResize(&data->xformTmp1[data->count], &data->bvhData[data->count]);
        TransformDataResize(&data->xformTmp2[data->count], &data->bvhData[data->count]);
        TransformDataResize(&data->xformTmp3[data->count], &data->bvhData[data->count]);

        snprintf(data->filePaths[data->count], 512, "%s", path);

        const char* filename = path;
        while (strchr(filename, '/')) { filename = strchr(filename, '/') + 1; }
        while (strchr(filename, '\\')) { filename = strchr(filename, '\\') + 1; }

        snprintf(data->names[data->count], 128, "%s", filename);
        data->scales[data->count] = 1.0f;

        // Auto-Scaling and unit detection

        if (data->bvhData[data->count].frameCount > 0)
        {
            TransformDataSampleFrame(&data->xformData[data->count], &data->bvhData[data->count], 0, 1.0f);
            TransformDataForwardKinematics(&data->xformData[data->count]);

            float height = 1e-8f;
            for (int j = 0; j < data->xformData[data->count].jointCount; j++)
            {
                height = Max(height, data->xformData[data->count].globalPositions[j].y);
            }

            data->scales[data->count] = height > 10.0f ? 0.01f : 1.0f;
            data->autoScales[data->count] = 1.8 / height;
        }
        else
        {
            data->autoScales[data->count] = 1.0f;
        }
        
        // Joint names combo

        int comboTotalSize = 0;
        for (int i = 0; i < data->bvhData[data->count].jointCount; i++)
        {
            comboTotalSize += (i > 0 ? 1 : 0) + strlen(data->bvhData[data->count].joints[i].name);
        }
        comboTotalSize++;

        data->jointNamesCombo[data->count] = malloc(comboTotalSize);
        data->jointNamesCombo[data->count][0] = '\0';
        for (int i = 0; i < data->bvhData[data->count].jointCount; i++)
        {
            if (i > 0)
            {
                strcat(data->jointNamesCombo[data->count], ";");
            }
            strcat(data->jointNamesCombo[data->count], data->bvhData[data->count].joints[i].name);
        }

        // Done

        data->count++;

        return true;
    }
    else
    {
        printf("INFO: Failed to Load '%s'\n", path);
        return false;
    }
}

//----------------------------------------------------------------------------------
// Geometric Functions
//----------------------------------------------------------------------------------

// Returns the time parameter along a line segment closest to another point
static inline float NearestPointOnLineSegment(
    Vector3 lineStart,
    Vector3 lineVector,
    Vector3 point)
{
    Vector3 ap = Vector3Subtract(point, lineStart);
    float lengthsq = Vector3LengthSqr(lineVector);
    return lengthsq < 1e-8f ? 0.5f : Saturate(Vector3DotProduct(lineVector, ap) / lengthsq);
}

// Returns the time parameters along two line segments at the closest point between the two
static inline void NearestPointBetweenLineSegments(
    float* nearestTime0,
    float* nearestTime1,
    Vector3 line0Start,
    Vector3 line0End,
    Vector3 line1Start,
    Vector3 line1End)
{
    Vector3 line0Vec = Vector3Subtract(line0End, line0Start);
    Vector3 line1Vec = Vector3Subtract(line1End, line1Start);
    float d0 = Vector3LengthSqr(Vector3Subtract(line1Start, line0Start));
    float d1 = Vector3LengthSqr(Vector3Subtract(line1End, line0Start));
    float d2 = Vector3LengthSqr(Vector3Subtract(line1Start, line0End));
    float d3 = Vector3LengthSqr(Vector3Subtract(line1End, line0End));

    *nearestTime0 = (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1) ? 1.0f : 0.0f;
    *nearestTime1 = NearestPointOnLineSegment(line1Start, line1Vec, Vector3Add(line0Start, Vector3Scale(line0Vec, *nearestTime0)));
    *nearestTime0 = NearestPointOnLineSegment(line0Start, line0Vec, Vector3Add(line1Start, Vector3Scale(line1Vec, *nearestTime1)));
}

// Returns the time parameter for a line segment closest to the plane
static inline float NearestPointBetweenLineSegmentAndPlane(Vector3 lineStart, Vector3 lineVector, Vector3 planePosition, Vector3 planeNormal)
{
    float denom = Vector3DotProduct(planeNormal, lineVector);
    if (fabs(denom) < 1e-8f)
    {
        return 0.5f;
    }
  
    return Saturate(Vector3DotProduct(Vector3Subtract(planePosition, lineStart), planeNormal) / denom);
}

// Returns the time parameter for a line segment closest to the ground plane
static inline float NearestPointBetweenLineSegmentAndGroundPlane(Vector3 lineStart, Vector3 lineVector)
{
    return fabs(lineVector.y) < 1e-8f ? 0.5f : Saturate((-lineStart.y) / lineVector.y);
}

// Returns the time parameter and nearest point on the ground between a line segment and ground segment 
static inline void NearestPointBetweenLineSegmentAndGroundSegment(
    float* nearestTimeOnLine,
    Vector3* nearestPointOnGround,
    Vector3 lineStart,
    Vector3 lineEnd,
    Vector3 groundMins,
    Vector3 groundMaxs)
{
    Vector3 lineVec = Vector3Subtract(lineEnd, lineStart);
  
    // Check Against Plane

    *nearestTimeOnLine = NearestPointBetweenLineSegmentAndGroundPlane(lineStart, lineVec);
    *nearestPointOnGround = (Vector3){
        lineStart.x + (*nearestTimeOnLine) * lineVec.x,
        0.0f,
        lineStart.z + (*nearestTimeOnLine) * lineVec.z,
    };

    // If point is inside plane bounds it must be the nearest

    if (nearestPointOnGround->x >= groundMins.x &&
        nearestPointOnGround->x <= groundMaxs.x &&
        nearestPointOnGround->z >= groundMins.z &&
        nearestPointOnGround->z <= groundMaxs.z)
    {
        return;
    }

    // Check against four edges

    Vector3 edgeStart0 =  (Vector3){ groundMins.x, 0.0f, groundMins.z };
    Vector3 edgeEnd0 = (Vector3){ groundMins.x, 0.0f, groundMaxs.z };
    
    Vector3 edgeStart1 = (Vector3){ groundMins.x, 0.0f, groundMaxs.z };
    Vector3 edgeEnd1 = (Vector3){ groundMaxs.x, 0.0f, groundMaxs.z };
    
    Vector3 edgeStart2 = (Vector3){ groundMaxs.x, 0.0f, groundMaxs.z };
    Vector3 edgeEnd2 = (Vector3){ groundMaxs.x, 0.0f, groundMins.z };
    
    Vector3 edgeStart3 = (Vector3){ groundMaxs.x, 0.0f, groundMins.z };
    Vector3 edgeEnd3 = (Vector3){ groundMins.x, 0.0f, groundMins.z };

    float nearestTimeOnLine0, nearestTimeOnLine1, nearestTimeOnLine2, nearestTimeOnLine3;
    float nearestTimeOnEdge0, nearestTimeOnEdge1, nearestTimeOnEdge2, nearestTimeOnEdge3;

    NearestPointBetweenLineSegments(
        &nearestTimeOnLine0,
        &nearestTimeOnEdge0,
        lineStart, lineEnd,
        edgeStart0, edgeEnd0);

    NearestPointBetweenLineSegments(
        &nearestTimeOnLine1,
        &nearestTimeOnEdge1,
        lineStart, lineEnd,
        edgeStart1, edgeEnd1);

    NearestPointBetweenLineSegments(
        &nearestTimeOnLine2,
        &nearestTimeOnEdge2,
        lineStart, lineEnd,
        edgeStart2, edgeEnd2);

    NearestPointBetweenLineSegments(
        &nearestTimeOnLine3,
        &nearestTimeOnEdge3,
        lineStart, lineEnd,
        edgeStart3, edgeEnd3);

    Vector3 nearestPointOnLine0 = Vector3Add(lineStart, Vector3Scale(lineVec, nearestTimeOnLine0));
    Vector3 nearestPointOnLine1 = Vector3Add(lineStart, Vector3Scale(lineVec, nearestTimeOnLine1));
    Vector3 nearestPointOnLine2 = Vector3Add(lineStart, Vector3Scale(lineVec, nearestTimeOnLine2));
    Vector3 nearestPointOnLine3 = Vector3Add(lineStart, Vector3Scale(lineVec, nearestTimeOnLine3));
    
    Vector3 nearestPointOnEdge0 = Vector3Add(edgeStart0, Vector3Scale(Vector3Subtract(edgeEnd0, edgeStart0), nearestTimeOnEdge0));
    Vector3 nearestPointOnEdge1 = Vector3Add(edgeStart1, Vector3Scale(Vector3Subtract(edgeEnd1, edgeStart1), nearestTimeOnEdge1));
    Vector3 nearestPointOnEdge2 = Vector3Add(edgeStart2, Vector3Scale(Vector3Subtract(edgeEnd2, edgeStart2), nearestTimeOnEdge2));
    Vector3 nearestPointOnEdge3 = Vector3Add(edgeStart3, Vector3Scale(Vector3Subtract(edgeEnd3, edgeStart3), nearestTimeOnEdge3));

    float distance0 = Vector3Distance(nearestPointOnLine0, nearestPointOnEdge0);
    float distance1 = Vector3Distance(nearestPointOnLine1, nearestPointOnEdge1);
    float distance2 = Vector3Distance(nearestPointOnLine2, nearestPointOnEdge2);
    float distance3 = Vector3Distance(nearestPointOnLine3, nearestPointOnEdge3);

    if (distance0 <= distance1 && distance0 <= distance2 && distance0 <= distance3)
    {
        *nearestTimeOnLine = nearestTimeOnLine0;
        *nearestPointOnGround = nearestPointOnEdge0;
        return;
    }

    if (distance1 <= distance0 && distance1 <= distance2 && distance1 <= distance3)
    {
        *nearestTimeOnLine = nearestTimeOnLine1;
        *nearestPointOnGround = nearestPointOnEdge1;
        return;
    }

    if (distance2 <= distance0 && distance2 <= distance1 && distance2 <= distance3)
    {
        *nearestTimeOnLine = nearestTimeOnLine2;
        *nearestPointOnGround = nearestPointOnEdge2;
        return;
    }

    if (distance3 <= distance0 && distance3 <= distance1 && distance3 <= distance2)
    {
        *nearestTimeOnLine = nearestTimeOnLine3;
        *nearestPointOnGround = nearestPointOnEdge3;
        return;
    }
    
    assert(false);
    *nearestTimeOnLine = nearestTimeOnLine0;
    *nearestPointOnGround = nearestPointOnEdge1;
    return;
}

static inline Vector3 ProjectPointOntoSweptLine(Vector3 sweptLineStart, Vector3 sweptLineVec, Vector3 sweptLineSweepVec, Vector3 position)
{
    Vector3 w = Vector3Subtract(position, sweptLineStart);
    Vector3 u = Vector3Normalize(sweptLineVec);
    Vector3 v = Vector3Normalize(sweptLineSweepVec);
    
    // x (u * u) + y (u * v) = w * u
    // x (v * u) + y (v * v) = w * v
    
    // Solved using Cramer's Rule in 2D
    float a1 = Vector3DotProduct(u, u);
    float b1 = Vector3DotProduct(u, v);
    float c1 = Vector3DotProduct(w, u);
    float a2 = Vector3DotProduct(v, u);
    float b2 = Vector3DotProduct(v, v);
    float c2 = Vector3DotProduct(w, v);
    
    float x = ((c1 * b2) - (b1 * c2)) / (a1 * b2 - b1 * a2);
    float y = (c1 - x * a1) / b1;
    
    x = Clamp(x, 0.0f, Vector3Length(sweptLineVec));
    y = Clamp(y, 0.0f, Vector3Length(sweptLineSweepVec));
    
    return Vector3Add(sweptLineStart, Vector3Add(Vector3Scale(u, x), Vector3Scale(v, y)));
}

// Returns the time parameter and nearest point on between a line segment and swept line segment
static inline void NearestPointBetweenLineSegmentAndSweptLine(
    float* nearestTimeOnLine,
    Vector3* nearestPointOnSweptLine,
    Vector3 lineStart,
    Vector3 lineEnd,
    Vector3 sweptLineStart,
    Vector3 sweptLineEnd,
    Vector3 sweptLineSweepVector)
{
    Vector3 lineVec = Vector3Subtract(lineEnd, lineStart);
    Vector3 sweptLineVec = Vector3Subtract(sweptLineEnd, sweptLineStart);
   
    Vector3 planeNormal = Vector3Length(sweptLineVec) < 1e-8f ? 
        Vector3Normalize(Vector3CrossProduct((Vector3){ 0.0f, 1.0f, 0.0f }, sweptLineSweepVector)) :
        Vector3Normalize(Vector3CrossProduct(sweptLineVec, sweptLineSweepVector));
    
    // Check Against Plane

    float nearestTimeOnLine0 = NearestPointBetweenLineSegmentAndPlane(lineStart, lineVec, sweptLineStart, planeNormal);
    Vector3 nearestPointOnLine0 = Vector3Add(lineStart, Vector3Scale(lineVec, nearestTimeOnLine0));
    
    Vector3 nearestPointOnSweptLine0;
    
    if (Vector3Length(sweptLineVec) > 1e-8f)
    {
        nearestPointOnSweptLine0 = ProjectPointOntoSweptLine(
            sweptLineStart, 
            sweptLineVec, 
            sweptLineSweepVector, 
            nearestPointOnLine0);
    }
    else
    {
        float nearestTimeOnSweptLine = NearestPointOnLineSegment(
            sweptLineStart,
            sweptLineSweepVector,
            nearestPointOnLine0);
            
        nearestPointOnSweptLine0 = Vector3Add(sweptLineStart, Vector3Scale(sweptLineSweepVector, nearestTimeOnSweptLine));
    }
    
    float distance0 = Vector3Distance(nearestPointOnLine0, nearestPointOnSweptLine0);
    
    // Check against three edges

    Vector3 edgeStart1 = sweptLineStart;
    Vector3 edgeEnd1 = Vector3Add(sweptLineStart, sweptLineSweepVector);
    
    Vector3 edgeStart2 = sweptLineEnd;
    Vector3 edgeEnd2 = Vector3Add(sweptLineEnd, sweptLineSweepVector);
    
    Vector3 edgeStart3 = sweptLineStart;
    Vector3 edgeEnd3 = sweptLineEnd;

    float nearestTimeOnLine1, nearestTimeOnLine2, nearestTimeOnLine3;
    float nearestTimeOnEdge1, nearestTimeOnEdge2, nearestTimeOnEdge3;

    NearestPointBetweenLineSegments(
        &nearestTimeOnLine1,
        &nearestTimeOnEdge1,
        lineStart, lineEnd,
        edgeStart1, edgeEnd1);

    NearestPointBetweenLineSegments(
        &nearestTimeOnLine2,
        &nearestTimeOnEdge2,
        lineStart, lineEnd,
        edgeStart2, edgeEnd2);

    NearestPointBetweenLineSegments(
        &nearestTimeOnLine3,
        &nearestTimeOnEdge3,
        lineStart, lineEnd,
        edgeStart3, edgeEnd3);

    Vector3 nearestPointOnLine1 = Vector3Add(lineStart, Vector3Scale(lineVec, nearestTimeOnLine1));
    Vector3 nearestPointOnLine2 = Vector3Add(lineStart, Vector3Scale(lineVec, nearestTimeOnLine2));
    Vector3 nearestPointOnLine3 = Vector3Add(lineStart, Vector3Scale(lineVec, nearestTimeOnLine3));
    
    Vector3 nearestPointOnSweptLine1 = Vector3Add(edgeStart1, Vector3Scale(Vector3Subtract(edgeEnd1, edgeStart1), nearestTimeOnEdge1));
    Vector3 nearestPointOnSweptLine2 = Vector3Add(edgeStart2, Vector3Scale(Vector3Subtract(edgeEnd2, edgeStart2), nearestTimeOnEdge2));
    Vector3 nearestPointOnSweptLine3 = Vector3Add(edgeStart3, Vector3Scale(Vector3Subtract(edgeEnd3, edgeStart3), nearestTimeOnEdge3));

    float distance1 = Vector3Distance(nearestPointOnLine1, nearestPointOnSweptLine1);
    float distance2 = Vector3Distance(nearestPointOnLine2, nearestPointOnSweptLine2);
    float distance3 = Vector3Distance(nearestPointOnLine3, nearestPointOnSweptLine3);

    if (distance0 <= distance1 && distance0 <= distance2 && distance0 <= distance3)
    {
        *nearestTimeOnLine = nearestTimeOnLine0;
        *nearestPointOnSweptLine = nearestPointOnSweptLine0;
        return;
    }

    if (distance1 <= distance0 && distance1 <= distance2 && distance1 <= distance3)
    {
        *nearestTimeOnLine = nearestTimeOnLine1;
        *nearestPointOnSweptLine = nearestPointOnSweptLine1;
        return;
    }

    if (distance2 <= distance0 && distance2 <= distance1 && distance2 <= distance3)
    {
        *nearestTimeOnLine = nearestTimeOnLine2;
        *nearestPointOnSweptLine = nearestPointOnSweptLine2;
        return;
    }

    if (distance3 <= distance0 && distance3 <= distance1 && distance3 <= distance2)
    {
        *nearestTimeOnLine = nearestTimeOnLine3;
        *nearestPointOnSweptLine = nearestPointOnSweptLine3;
        return;
    }
    
    // Unreachable
    assert(false);
    *nearestTimeOnLine = nearestTimeOnLine0;
    *nearestPointOnSweptLine = nearestPointOnSweptLine0;
    return;
}

// Analytical capsule and sphere occlusion functions taken from here:
// https://www.shadertoy.com/view/3stcD4


// This is the number of times the radius away from the sphere where
// the ambient occlusion drops off to zero. This is important for various
// acceleration methods to filter out capsules which are too far away and
// so not casting any ambient occlusion
#define AO_RATIO_MAX 4.0

static inline float SphereOcclusionLookup(float nlAngle, float h)
{
    float nl = cosf(nlAngle);
    float h2 = h*h;
    
    float res = Max(nl, 0.0) / h2;
    float k2 = 1.0 - h2*nl*nl;
    if (k2 > 1e-4f)
    {
        res = nl * acosf(Clamp(-nl*sqrtf((h2 - 1.0f) / Max(1.0f - nl*nl, 1e-8f)), -1.0f, 1.0f)) - sqrtf(k2*(h2 - 1.0f));
        res = (res / h2 + atanf(sqrt(k2 / (h2 - 1.0f)))) / PI;
    }

    float decay = Max(1.0f - (h - 1.0f) / ((float)AO_RATIO_MAX - 1.0f), 0.0f);
    
    return 1.0f - res * decay;
}

static inline float SphereOcclusion(Vector3 pos, Vector3 nor, Vector3 sph, float rad)
{
    Vector3 di = Vector3Subtract(sph, pos);
    float l = Vector3Length(di);
    float nlAngle = acosf(Clamp(Vector3DotProduct(nor, Vector3Scale(di, 1.0f / Max(l, 1e-8f))), -1.0f, 1.0f));
    float h  = l < rad ? 1.0 : l / rad;
    return SphereOcclusionLookup(nlAngle, h);
}

static inline float SphereIntersectionArea(float r1, float r2, float d)
{
    if (Min(r1, r2) <= Max(r1, r2) - d)
    {
        return 1.0f - Max(cosf(r1), cosf(r2));
    }
    else if (r1 + r2 <= d)
    {
        return 0.0f;
    }

    float delta = fabs(r1 - r2);
    float x = 1.0f - Saturate((d - delta) / Max(r1 + r2 - delta, 1e-8f));
    float area = Square(x) * (-2.0f * x + 3.0f);

    return area * (1.0f - Max(cosf(r1), cosf(r2)));
}

static inline float SphereDirectionalOcclusionLookup(float phi, float theta, float coneAngle)
{
    return 1.0f - SphereIntersectionArea(theta, coneAngle / 2.0f, phi) / (1.0f - cosf(coneAngle / 2.0f));
}

static inline float SphereDirectionalOcclusion(
    Vector3 pos, 
    Vector3 sphere, 
    float radius,
    Vector3 coneDir, 
    float coneAngle)
{
    Vector3 occluder = Vector3Subtract(sphere, pos);
    float occluderLen2 = Vector3DotProduct(occluder, occluder);
    Vector3 occluderDir = Vector3Scale(occluder, 1.0f / Max(sqrtf(occluderLen2), 1e-8f));

    float phi = acosf(Clamp(Vector3DotProduct(occluderDir, Vector3Negate(coneDir)), -1.0f, 1.0f));
    float theta = acosf(Clamp(sqrtf(occluderLen2 / (Square(radius) + occluderLen2)), -1.0f, 1.0f));
    
    return SphereDirectionalOcclusionLookup(phi, theta, coneAngle);
}

// Get the start point of the capsule line segment
static inline Vector3 CapsuleStart(Vector3 capsulePosition, Quaternion capsuleRotation, float capsuleHalfLength)
{
    return Vector3Add(capsulePosition,
        Vector3RotateByQuaternion((Vector3){+capsuleHalfLength, 0.0f, 0.0f}, capsuleRotation));
}

// Get the end point of the capsule line segment
static inline Vector3 CapsuleEnd(Vector3 capsulePosition, Quaternion capsuleRotation, float capsuleHalfLength)
{
    return Vector3Add(capsulePosition,
        Vector3RotateByQuaternion((Vector3){-capsuleHalfLength, 0.0f, 0.0f}, capsuleRotation));
}

// Get the vector from the start to the end of the capsule line segment
static inline Vector3 CapsuleVector(Vector3 capsulePosition, Quaternion capsuleRotation, float capsuleHalfLength)
{
    Vector3 capsuleStart = CapsuleStart(capsulePosition, capsuleRotation, capsuleHalfLength);

    return Vector3Subtract(Vector3Add(capsulePosition,
        Vector3RotateByQuaternion((Vector3){-capsuleHalfLength, 0.0f, 0.0f}, capsuleRotation)), capsuleStart);
}

static inline float CapsuleDirectionalOcclusion(
    Vector3 pos, Vector3 capStart, Vector3 capVec,
    float capRadius, Vector3 coneDir, float coneAngle)
{
    Vector3 ba = capVec;
    Vector3 pa = Vector3Subtract(capStart, pos);
    Vector3 cba = Vector3Subtract(Vector3Scale(Vector3Negate(coneDir), Vector3DotProduct(Vector3Negate(coneDir), ba)), ba);
    float t = Saturate(Vector3DotProduct(pa, cba) / Max(Vector3DotProduct(cba, cba), 1e-8f));

    return SphereDirectionalOcclusion(pos, Vector3Add(capStart, Vector3Scale(ba, t)), capRadius, coneDir, coneAngle);
}

//----------------------------------------------------------------------------------
// Capsule Data
//----------------------------------------------------------------------------------

// Basic type useful for sorting according to some value
typedef struct
{
    int index;
    float value;
} CapsuleSort;

static inline int CapsuleSortCompareGreater(const void* lhs, const void* rhs)
{
    const CapsuleSort* lhsSort = lhs;
    const CapsuleSort* rhsSort = rhs;
    return lhsSort->value > rhsSort->value ? 1 : -1;
}

static inline int CapsuleSortCompareLess(const void* lhs, const void* rhs)
{
    const CapsuleSort* lhsSort = lhs;
    const CapsuleSort* rhsSort = rhs;
    return lhsSort->value < rhsSort->value ? 1 : -1;
}

// Structure containing all of the data required for all of the capsules which are to be rendered.
typedef struct
{
    // Data for all the capsules which are in the scene
    int capsuleCount;
    Vector3* capsulePositions;
    Quaternion* capsuleRotations;
    float* capsuleRadii;
    float* capsuleHalfLengths;
    Vector3* capsuleColors;
    float* capsuleOpacities;
    CapsuleSort* capsuleSort;

    // Buffers for all the capsules casting ambient occlusion
    int aoCapsuleCount;
    Vector3* aoCapsuleStarts;
    Vector3* aoCapsuleVectors;
    float* aoCapsuleRadii;
    CapsuleSort* aoCapsuleSort;

    // Buffers for all the capsules casting shadows
    int shadowCapsuleCount;
    Vector3* shadowCapsuleStarts;
    Vector3* shadowCapsuleVectors;
    float* shadowCapsuleRadii;
    CapsuleSort* shadowCapsuleSort;
    
    // Lookup table for the capsule ambient occlusion function
    Image aoLookupImage;
    Texture2D aoLookupTable;
    Vector2 aoLookupResolution;

    // Lookup table for the capsule shadow function
    Image shadowLookupImage;
    Texture2D shadowLookupTable;
    Vector2 shadowLookupResolution;

} CapsuleData;

// Compute the capsule ambient occlusion lookup table
static void CapsuleDataUpdateAOLookupTable(CapsuleData* data)
{
    int width = (int)data->aoLookupResolution.x;
    int height = (int)data->aoLookupResolution.y;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            float nlAngle = (((float)x) / (width - 1)) * PI;
            float h = 1.0f + (AO_RATIO_MAX - 1.0f) * (((float)y) / (height - 1));
            ((unsigned char*)data->aoLookupImage.data)[y * width + x] = (unsigned char)Clamp(255.0 * SphereOcclusionLookup(nlAngle, h), 0.0, 255.0);
        }
    }

    UpdateTexture(data->aoLookupTable, data->aoLookupImage.data);
}

// Compute the capsule shadow lookup table for a given coneAngle
static void CapsuleDataUpdateShadowLookupTable(CapsuleData* data, float coneAngle)
{
    int width = (int)data->shadowLookupResolution.x;
    int height = (int)data->shadowLookupResolution.y;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            float phi = (((float)x) / (width - 1)) * PI;
            float theta = ((float)y) / (height - 1) * (PI / 2.0f);
            ((unsigned char*)data->shadowLookupImage.data)[y * width + x] = (unsigned char)Clamp(255.0 * SphereDirectionalOcclusionLookup(phi, theta, coneAngle), 0.0, 255.0);
        }
    }
    
    UpdateTexture(data->shadowLookupTable, data->shadowLookupImage.data);
}

static void CapsuleDataInit(CapsuleData* data)
{
    // Init
  
    data->capsuleCount = 0;
    data->capsulePositions = NULL;
    data->capsuleRotations = NULL;
    data->capsuleRadii = NULL;
    data->capsuleHalfLengths = NULL;
    data->capsuleColors = NULL;
    data->capsuleOpacities = NULL;
    data->capsuleSort = NULL;

    data->aoCapsuleCount = 0;
    data->aoCapsuleStarts = NULL;
    data->aoCapsuleVectors = NULL;
    data->aoCapsuleRadii = NULL;
    data->aoCapsuleSort = NULL;

    data->shadowCapsuleCount = 0;
    data->shadowCapsuleStarts = NULL;
    data->shadowCapsuleVectors = NULL;
    data->shadowCapsuleRadii = NULL;
    data->shadowCapsuleSort = NULL;
    
    // Capsule AO Lookup Table
    
    data->aoLookupImage = (Image){
        .data = calloc(32 * 32, 1),
        .width = 32,
        .height = 32,
        .format = PIXELFORMAT_UNCOMPRESSED_GRAYSCALE,
        .mipmaps = 1
    };
    
    data->aoLookupTable = LoadTextureFromImage(data->aoLookupImage);
    data->aoLookupResolution = (Vector2){ (float)data->aoLookupImage.width, (float)data->aoLookupImage.height };
    SetTextureWrap(data->aoLookupTable, TEXTURE_WRAP_CLAMP);
    SetTextureFilter(data->aoLookupTable, TEXTURE_FILTER_BILINEAR);
    
    CapsuleDataUpdateAOLookupTable(data);
    
    // Capsule Shadow Lookup Table
    
    data->shadowLookupImage = (Image){
        .data = calloc(256 * 128, 1),
        .width = 256,
        .height = 128,
        .format = PIXELFORMAT_UNCOMPRESSED_GRAYSCALE,
        .mipmaps = 1
    };
    
    data->shadowLookupTable = LoadTextureFromImage(data->shadowLookupImage);
    data->shadowLookupResolution = (Vector2){ (float)data->shadowLookupImage.width, (float)data->shadowLookupImage.height };
    SetTextureWrap(data->shadowLookupTable, TEXTURE_WRAP_CLAMP);
    SetTextureFilter(data->shadowLookupTable, TEXTURE_FILTER_BILINEAR);
    
    CapsuleDataUpdateShadowLookupTable(data, 0.2f);
}

static void CapsuleDataResize(CapsuleData* data, int maxCapsuleCount)
{
    data->capsulePositions = realloc(data->capsulePositions, maxCapsuleCount * sizeof(Vector3));
    data->capsuleRotations = realloc(data->capsuleRotations, maxCapsuleCount * sizeof(Quaternion));
    data->capsuleRadii = realloc(data->capsuleRadii, maxCapsuleCount * sizeof(float));
    data->capsuleHalfLengths = realloc(data->capsuleHalfLengths, maxCapsuleCount * sizeof(float));
    data->capsuleColors = realloc(data->capsuleColors, maxCapsuleCount * sizeof(Vector3));
    data->capsuleOpacities = realloc(data->capsuleOpacities, maxCapsuleCount * sizeof(float));
    data->capsuleSort = realloc(data->capsuleSort, maxCapsuleCount * sizeof(CapsuleSort));

    data->aoCapsuleStarts = realloc(data->aoCapsuleStarts, maxCapsuleCount * sizeof(Vector3));
    data->aoCapsuleVectors = realloc(data->aoCapsuleVectors, maxCapsuleCount * sizeof(Vector3));
    data->aoCapsuleRadii = realloc(data->aoCapsuleRadii, maxCapsuleCount * sizeof(float));
    data->aoCapsuleSort = realloc(data->aoCapsuleSort, maxCapsuleCount * sizeof(CapsuleSort));

    data->shadowCapsuleStarts = realloc(data->shadowCapsuleStarts, maxCapsuleCount * sizeof(Vector3));
    data->shadowCapsuleVectors = realloc(data->shadowCapsuleVectors, maxCapsuleCount * sizeof(Vector3));
    data->shadowCapsuleRadii = realloc(data->shadowCapsuleRadii, maxCapsuleCount * sizeof(float));
    data->shadowCapsuleSort = realloc(data->shadowCapsuleSort, maxCapsuleCount * sizeof(CapsuleSort));
}

static void CapsuleDataFree(CapsuleData* data)
{
    free(data->capsulePositions);
    free(data->capsuleRotations);
    free(data->capsuleRadii);
    free(data->capsuleHalfLengths);
    free(data->capsuleColors);
    free(data->capsuleOpacities);
    free(data->capsuleSort);

    free(data->aoCapsuleStarts);
    free(data->aoCapsuleVectors);
    free(data->aoCapsuleRadii);
    free(data->aoCapsuleSort);

    free(data->shadowCapsuleStarts);
    free(data->shadowCapsuleVectors);
    free(data->shadowCapsuleRadii);
    free(data->shadowCapsuleSort);
    
    UnloadImage(data->aoLookupImage);
    UnloadTexture(data->aoLookupTable);
    UnloadImage(data->shadowLookupImage);
    UnloadTexture(data->shadowLookupTable);
}

static inline void CapsuleDataReset(CapsuleData* data)
{
    data->capsuleCount = 0;
    data->aoCapsuleCount = 0;
    data->shadowCapsuleCount = 0;
}

// Append capsules to the capsule data based off the joint transforms
static void CapsuleDataAppendFromTransformData(CapsuleData* data, TransformData* xforms, float maxCapsuleRadius, Color color, float opacity, bool ignoreEndSite)
{
    for (int i = 0; i < xforms->jointCount; i++)
    {
        int p = xforms->parents[i];
        
        if (p == -1) { continue; }
        if (ignoreEndSite && xforms->endSite[i]) { continue; }

        float capsuleHalfLength = Vector3Length(xforms->localPositions[i]) / 2.0f;
        float capsuleRadius = Min(maxCapsuleRadius, capsuleHalfLength) + (i % 2) * 0.001f;

        if (capsuleRadius < 0.001f) { continue; }

        Vector3 capsulePosition = Vector3Scale(Vector3Add(xforms->globalPositions[i], xforms->globalPositions[p]), 0.5f);
        Quaternion capsuleRotation = QuaternionMultiply(
            xforms->globalRotations[p],
            QuaternionBetween((Vector3){ 1.0f, 0.0f, 0.0f }, Vector3Normalize(xforms->localPositions[i])));

        data->capsulePositions[data->capsuleCount] = capsulePosition;
        data->capsuleRotations[data->capsuleCount] = capsuleRotation;
        data->capsuleHalfLengths[data->capsuleCount] = capsuleHalfLength;
        data->capsuleRadii[data->capsuleCount] = capsuleRadius;
        data->capsuleColors[data->capsuleCount] = (Vector3){ color.r / 255.0f, color.g / 255.0f, color.b / 255.0f };
        data->capsuleOpacities[data->capsuleCount] = opacity;
        data->capsuleCount++;
    }
}

// Gather all of the capsules which are potentially casting ambient occlusion on a ground segment.
static void CapsuleDataUpdateAOCapsulesForGroundSegment(CapsuleData* data, Vector3 groundSegmentPosition)
{
    data->aoCapsuleCount = 0;

    for (int i = 0; i < data->capsuleCount; i++)
    {
        Vector3 capsulePosition = data->capsulePositions[i];
        float capsuleHalfLength = data->capsuleHalfLengths[i];
        float capsuleRadius = data->capsuleRadii[i];
        
        // Check if bounding spheres are more than AO_RATIO_MAX away from each other
        if (Vector3Distance(groundSegmentPosition, capsulePosition) - sqrtf(2.0f) > capsuleHalfLength + AO_RATIO_MAX * capsuleRadius)
        {
            continue;
        }
        
        Quaternion capsuleRotation = data->capsuleRotations[i];
        Vector3 capsuleStart = CapsuleStart(capsulePosition, capsuleRotation, capsuleHalfLength);
        Vector3 capsuleEnd = CapsuleEnd(capsulePosition, capsuleRotation, capsuleHalfLength);
        Vector3 capsuleVector = CapsuleVector(capsulePosition, capsuleRotation, capsuleHalfLength);

        float capsuleTime;
        Vector3 groundPoint;
        NearestPointBetweenLineSegmentAndGroundSegment(
            &capsuleTime,
            &groundPoint,
            capsuleStart,
            capsuleEnd,
            (Vector3){ groundSegmentPosition.x - 1.0f, 0.0f, groundSegmentPosition.z - 1.0f },
            (Vector3){ groundSegmentPosition.x + 1.0f, 0.0f, groundSegmentPosition.z + 1.0f });
    
        Vector3 capsulePoint = Vector3Add(capsuleStart, Vector3Scale(capsuleVector, capsuleTime));
        
        // Check if the nearest point on the ground is more than AO_RATIO_MAX away
        if (Vector3Distance(groundPoint, capsulePoint) > AO_RATIO_MAX * capsuleRadius)
        {
            continue;
        }
        
        // Compute the actual occlusion for the closest point on the ground
        float capsuleOcclusion = Vector3Distance(groundPoint, capsulePoint) < capsuleRadius ? 0.0f :
            SphereOcclusion(groundPoint, (Vector3){ 0.0f, 1.0f, 0.0f }, capsulePoint, capsuleRadius);

        if (capsuleOcclusion < 0.99f)
        {
            data->aoCapsuleSort[data->aoCapsuleCount] = (CapsuleSort){ i, capsuleOcclusion };
            data->aoCapsuleCount++;
        }
    }

    qsort(data->aoCapsuleSort, data->aoCapsuleCount, sizeof(CapsuleSort), CapsuleSortCompareGreater);

    for (int i = 0; i < data->aoCapsuleCount; i++)
    {
        int j = data->aoCapsuleSort[i].index;
        data->aoCapsuleStarts[i] = CapsuleStart(data->capsulePositions[j], data->capsuleRotations[j], data->capsuleHalfLengths[j]);
        data->aoCapsuleVectors[i] = CapsuleVector(data->capsulePositions[j], data->capsuleRotations[j], data->capsuleHalfLengths[j]);
        data->aoCapsuleRadii[i] = data->capsuleRadii[j];
    }
}

// Gather all of the capsules which are potentially casting ambient occlusion on another capsule.
static void CapsuleDataUpdateAOCapsulesForCapsule(CapsuleData* data, int capsuleIndex)
{
    Vector3 queryCapsulePosition = data->capsulePositions[capsuleIndex];
    float queryCapsuleHalfLength = data->capsuleHalfLengths[capsuleIndex];
    float queryCapsuleRadius = data->capsuleRadii[capsuleIndex];
    Quaternion queryCapsuleRotation = data->capsuleRotations[capsuleIndex];
    Vector3 queryCapsuleStart = CapsuleStart(queryCapsulePosition, queryCapsuleRotation, queryCapsuleHalfLength);
    Vector3 queryCapsuleEnd = CapsuleEnd(queryCapsulePosition, queryCapsuleRotation, queryCapsuleHalfLength);
    Vector3 queryCapsuleVector = CapsuleVector(queryCapsulePosition, queryCapsuleRotation, queryCapsuleHalfLength);

    data->aoCapsuleCount = 0;

    for (int i = 0; i < data->capsuleCount; i++)
    {
        if (i == capsuleIndex) { continue; }
        
        Vector3 capsulePosition = data->capsulePositions[i];
        float capsuleRadius = data->capsuleRadii[i];
        float capsuleHalfLength = data->capsuleHalfLengths[i];

        // Check if the bounding sphers are more than AO_RATIO_MAX away from each other
        if (Vector3Distance(queryCapsulePosition, capsulePosition) - queryCapsuleHalfLength - queryCapsuleRadius > 
            capsuleHalfLength + AO_RATIO_MAX * capsuleRadius)
        {
            continue;
        }

        Quaternion capsuleRotation = data->capsuleRotations[i];
        Vector3 capsuleStart = CapsuleStart(capsulePosition, capsuleRotation, capsuleHalfLength);
        Vector3 capsuleEnd = CapsuleEnd(capsulePosition, capsuleRotation, capsuleHalfLength);
        Vector3 capsuleVector = CapsuleVector(capsulePosition, capsuleRotation, capsuleHalfLength);
        
        float capsuleTime, queryTime;
        NearestPointBetweenLineSegments(
            &capsuleTime,
            &queryTime,
            capsuleStart,
            capsuleEnd,
            queryCapsuleStart,
            queryCapsuleEnd);

        Vector3 capsulePoint = Vector3Add(capsuleStart, Vector3Scale(capsuleVector, capsuleTime));
        Vector3 queryPoint = Vector3Add(queryCapsuleStart, Vector3Scale(queryCapsuleVector, queryTime));
        
        // Check if the nearest points on the two capsules are more than AO_RATIO_MAX away
        if (Vector3Distance(queryPoint, capsulePoint) - queryCapsuleRadius > AO_RATIO_MAX * capsuleRadius)
        {
            continue;
        }
        
        // Compute the actual occlusion at the nearest point
        Vector3 surfaceNormal = Vector3Normalize(Vector3Subtract(capsulePoint, queryPoint));
        Vector3 surfacePoint = Vector3Add(queryPoint, Vector3Scale(surfaceNormal, queryCapsuleRadius));
        float capsuleOcclusion = Vector3Distance(queryPoint, capsulePoint) <= queryCapsuleRadius + capsuleRadius ? 0.0f :
            SphereOcclusion(surfacePoint, surfaceNormal, capsulePoint, capsuleRadius);

        if (capsuleOcclusion < 0.99f)
        {
            data->aoCapsuleSort[data->aoCapsuleCount] = (CapsuleSort){ i, capsuleOcclusion };
            data->aoCapsuleCount++;
        }
    }

    qsort(data->aoCapsuleSort, data->aoCapsuleCount, sizeof(CapsuleSort), CapsuleSortCompareGreater);

    for (int i = 0; i < data->aoCapsuleCount; i++)
    {
        int j = data->aoCapsuleSort[i].index;
        data->aoCapsuleStarts[i] = CapsuleStart(data->capsulePositions[j], data->capsuleRotations[j], data->capsuleHalfLengths[j]);
        data->aoCapsuleVectors[i] = CapsuleVector(data->capsulePositions[j], data->capsuleRotations[j], data->capsuleHalfLengths[j]);
        data->aoCapsuleRadii[i] = data->capsuleRadii[j];
    }
}

// Gather all of the capsules which are potentially casting shadows on a ground segment.
static void CapsuleDataUpdateShadowCapsulesForGroundSegment(CapsuleData* data, Vector3 groundSegmentPosition, Vector3 lightDir, float lightConeAngle)
{
    Vector3 lightRay = Vector3Scale(lightDir, 10.0f);

    data->shadowCapsuleCount = 0;

    for (int i = 0; i < data->capsuleCount; i++)
    {
        Vector3 capsulePosition = data->capsulePositions[i];
        float capsuleHalfLength = data->capsuleHalfLengths[i];
        float capsuleRadius = data->capsuleRadii[i];
        
        float midRayTime = NearestPointBetweenLineSegmentAndGroundPlane(capsulePosition, lightRay);
        Vector3 groundCapsuleMid = Vector3Add(capsulePosition, Vector3Scale(lightRay, midRayTime));
        float maxRatio = 4.0f; // This is a kind of fudge-factor as the soft shadow don't have a fixed falloff

        // Check if the ground segment is more than maxRatio away from the shadow point at the center of the capsule 
        if (Vector3Distance(groundSegmentPosition, groundCapsuleMid) - sqrtf(2.0f) > capsuleHalfLength + maxRatio * capsuleRadius)
        {
            continue;
        }
      
        Quaternion capsuleRotation = data->capsuleRotations[i];
        Vector3 capsuleStart = CapsuleStart(capsulePosition, capsuleRotation, capsuleHalfLength);
        Vector3 capsuleEnd = CapsuleEnd(capsulePosition, capsuleRotation, capsuleHalfLength);
        Vector3 capsuleVector = CapsuleVector(capsulePosition, capsuleRotation, capsuleHalfLength);

        // Find the projected shadow point for the capsule start and end points
        // I think for the ground the darkest part of the shadow is always one of these
        float startRayTime = NearestPointBetweenLineSegmentAndGroundPlane(capsuleStart, lightRay);
        float endRayTime = NearestPointBetweenLineSegmentAndGroundPlane(capsuleEnd, lightRay);

        Vector3 groundCapsuleStart = Vector3Add(capsuleStart, Vector3Scale(lightRay, startRayTime));
        Vector3 groundCapsuleEnd = Vector3Add(capsuleEnd, Vector3Scale(lightRay, endRayTime));
        
        groundCapsuleStart.x = Clamp(groundCapsuleStart.x, groundSegmentPosition.x - 1.0f, groundSegmentPosition.x + 1.0f);
        groundCapsuleStart.z = Clamp(groundCapsuleStart.z, groundSegmentPosition.z - 1.0f, groundSegmentPosition.z + 1.0f);
        groundCapsuleEnd.x = Clamp(groundCapsuleEnd.x, groundSegmentPosition.x - 1.0f, groundSegmentPosition.x + 1.0f);
        groundCapsuleEnd.z = Clamp(groundCapsuleEnd.z, groundSegmentPosition.z - 1.0f, groundSegmentPosition.z + 1.0f);
        
        // Check if both points are more than maxRatio away from the ground segment
        if (Vector3Distance(groundSegmentPosition, groundCapsuleStart) - sqrtf(2.0f) > maxRatio * capsuleRadius &&
            Vector3Distance(groundSegmentPosition, groundCapsuleEnd) - sqrtf(2.0f) > maxRatio * capsuleRadius)
        {
            continue;
        }
        
        // Compute the actual occlusion at both points and take the min
        float capsuleOcclusion = Min(
            CapsuleDirectionalOcclusion(groundCapsuleStart, capsuleStart, capsuleVector, capsuleRadius, lightDir, lightConeAngle),
            CapsuleDirectionalOcclusion(groundCapsuleEnd, capsuleStart, capsuleVector, capsuleRadius, lightDir, lightConeAngle));

        if (capsuleOcclusion < 0.99f)
        {
            data->shadowCapsuleSort[data->shadowCapsuleCount] = (CapsuleSort){ i, capsuleOcclusion };
            data->shadowCapsuleCount++;
        }
    }

    qsort(data->shadowCapsuleSort, data->shadowCapsuleCount, sizeof(CapsuleSort), CapsuleSortCompareGreater);

    for (int i = 0; i < data->shadowCapsuleCount; i++)
    {
        int j = data->shadowCapsuleSort[i].index;
        data->shadowCapsuleStarts[i] = CapsuleStart(data->capsulePositions[j], data->capsuleRotations[j], data->capsuleHalfLengths[j]);
        data->shadowCapsuleVectors[i] = CapsuleVector(data->capsulePositions[j], data->capsuleRotations[j], data->capsuleHalfLengths[j]);
        data->shadowCapsuleRadii[i] = data->capsuleRadii[j];
    }
}

// Gather all of the capsules which are potentially casting shadows on another capsule.
static void CapsuleDataUpdateShadowCapsulesForCapsule(CapsuleData* data, int capsuleIndex, Vector3 lightDir, float lightConeAngle)
{
    Vector3 lightRay = Vector3Scale(lightDir, 10.0f);
    
    Vector3 queryCapsulePosition = data->capsulePositions[capsuleIndex];
    float queryCapsuleHalfLength = data->capsuleHalfLengths[capsuleIndex];
    float queryCapsuleRadius = data->capsuleRadii[capsuleIndex];
    Quaternion queryCapsuleRotation = data->capsuleRotations[capsuleIndex];
    Vector3 queryCapsuleStart = CapsuleStart(queryCapsulePosition, queryCapsuleRotation, queryCapsuleHalfLength);
    Vector3 queryCapsuleEnd = CapsuleEnd(queryCapsulePosition, queryCapsuleRotation, queryCapsuleHalfLength);
    Vector3 queryCapsuleVector = CapsuleVector(queryCapsulePosition, queryCapsuleRotation, queryCapsuleHalfLength);

    data->shadowCapsuleCount = 0;

    for (int i = 0; i < data->capsuleCount; i++)
    {
        if (i == capsuleIndex) { continue; }
        
        Vector3 capsulePosition = data->capsulePositions[i];
        float capsuleHalfLength = data->capsuleHalfLengths[i];
        float capsuleRadius = data->capsuleRadii[i];
        
        // Find the closest point between the capsule and the ray cast from the center of the casting capsule
        float midRayTime = NearestPointOnLineSegment(
            capsulePosition,
            lightRay,
            queryCapsulePosition);
        
        Vector3 capsuleMid = Vector3Add(capsulePosition, Vector3Scale(lightRay, midRayTime));
        float maxRatio = 4.0f;
        
        // Check if this is greater than maxRatio away
        if (Vector3Distance(queryCapsulePosition, capsuleMid) - queryCapsuleHalfLength - queryCapsuleRadius > capsuleHalfLength + maxRatio * capsuleRadius)
        {
            continue;
        }
        
        Quaternion capsuleRotation = data->capsuleRotations[i];
        Vector3 capsuleStart = CapsuleStart(capsulePosition, capsuleRotation, capsuleHalfLength);
        Vector3 capsuleEnd = CapsuleEnd(capsulePosition, capsuleRotation, capsuleHalfLength);
        Vector3 capsuleVector = CapsuleVector(capsulePosition, capsuleRotation, capsuleHalfLength);
        
        // Find the nearest point between the capsule and the swept line of the casting capsule
        float queryCapsuleTime;
        Vector3 nearestRayPoint;
        NearestPointBetweenLineSegmentAndSweptLine(
            &queryCapsuleTime,
            &nearestRayPoint,
            queryCapsuleStart,
            queryCapsuleEnd,
            capsuleStart,
            capsuleEnd,
            lightRay);
        
        Vector3 queryCapsulePoint = Vector3Add(queryCapsuleStart, Vector3Scale(queryCapsuleVector, queryCapsuleTime));
        
        // If this distance is greater than maxRatio away then skip
        if (Vector3Distance(queryCapsulePoint, nearestRayPoint) - queryCapsuleRadius > capsuleHalfLength + maxRatio * capsuleRadius)
        {
            continue;
        }
        
        Vector3 surfaceNormal = Vector3Normalize(Vector3Subtract(nearestRayPoint, queryCapsulePoint));
        Vector3 surfacePoint = Vector3Add(queryCapsulePoint, Vector3Scale(surfaceNormal, queryCapsuleRadius));

        // Find actual occlusion amount
        float capsuleOcclusion = Vector3Distance(queryCapsulePoint, nearestRayPoint) <= queryCapsuleRadius + capsuleRadius ? 0.0f :
            CapsuleDirectionalOcclusion(surfacePoint, capsuleStart, capsuleVector, capsuleRadius, lightDir, lightConeAngle);
        
        if (capsuleOcclusion < 0.99f)
        {
            data->shadowCapsuleSort[data->shadowCapsuleCount] = (CapsuleSort){ i, capsuleOcclusion };
            data->shadowCapsuleCount++;
        }
    }

    qsort(data->shadowCapsuleSort, data->shadowCapsuleCount, sizeof(CapsuleSort), CapsuleSortCompareGreater);

    for (int i = 0; i < data->shadowCapsuleCount; i++)
    {
        int j = data->shadowCapsuleSort[i].index;
        data->shadowCapsuleStarts[i] = CapsuleStart(data->capsulePositions[j], data->capsuleRotations[j], data->capsuleHalfLengths[j]);
        data->shadowCapsuleVectors[i] = CapsuleVector(data->capsulePositions[j], data->capsuleRotations[j], data->capsuleHalfLengths[j]);
        data->shadowCapsuleRadii[i] = data->capsuleRadii[j];
    }
}

// Resize so that we have enough capsules in the buffers for the given set of characters
static inline void CapsuleDataUpdateForCharacters(CapsuleData* capsuleData, CharacterData* characterData)
{
    int totalJointCount = 0;
    for (int i = 0; i < characterData->count; i++)
    {
        totalJointCount += characterData->bvhData[i].jointCount;
    }

    CapsuleDataResize(capsuleData, totalJointCount);
}

//----------------------------------------------------------------------------------
// Shaders
//----------------------------------------------------------------------------------

#define AO_CAPSULES_MAX 32
#define SHADOW_CAPSULES_MAX 64

#define GLSL_DEFINE_VALUE(X) #X
#define GLSL_DEFINE(X) "#define " #X " " GLSL_DEFINE_VALUE(X) " \n"
#define GLSL(X) \
  "#version 300 es\n" \
  GLSL_DEFINE(AO_RATIO_MAX) \
  GLSL_DEFINE(AO_CAPSULES_MAX) \
  GLSL_DEFINE(SHADOW_CAPSULES_MAX) \
  GLSL_DEFINE(PI) \
  #X

// Vertex Shader
static const char* shaderVS = GLSL(

precision highp float;
precision mediump int;

in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;

uniform int isCapsule;
uniform vec3 capsulePosition;
uniform vec4 capsuleRotation;
uniform float capsuleHalfLength;
uniform float capsuleRadius;

uniform mat4 matModel;
uniform mat4 matNormal;
uniform mat4 matView;
uniform mat4 matProjection;

uniform int isCharacter;

out vec3 fragPosition;
out vec2 fragTexCoord;
out vec3 fragNormal;

vec3 Rotate(in vec4 q, vec3 v)
{
    vec3 t = 2.0 * cross(q.xyz, v);
    return v + q.w * t + cross(q.xyz, t);
}

// Stretch capsule according to capsule half length
vec3 CapsuleStretch(vec3 pos, float hlength, float radius)
{
    vec3 scaled = pos * radius;
    scaled.x = scaled.x > 0.0 ? scaled.x + hlength : scaled.x - hlength;
    return scaled;
}

void main()
{
    fragTexCoord = vertexTexCoord;

    if (isCapsule == 1)
    {
        fragPosition = Rotate(capsuleRotation,
            CapsuleStretch(vertexPosition,
            capsuleHalfLength, capsuleRadius)) + capsulePosition;

        fragNormal = Rotate(capsuleRotation, vertexNormal);
    }
    else if (isCharacter == 1)
    {
        fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
        fragNormal = normalize(vec3(matNormal * vec4(vertexNormal, 1.0)));
    }
    else
    {
        fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
        fragNormal = normalize(vec3(matNormal * vec4(vertexNormal, 1.0)));
    }

    gl_Position = matProjection * matView * vec4(fragPosition, 1.0);
}

);

// Fragment Shader
static const char* shaderFS = GLSL(

precision highp float;
precision mediump int;

in vec3 fragPosition;
in vec2 fragTexCoord;
in vec3 fragNormal;

uniform vec3 objectColor;
uniform float objectSpecularity;
uniform float objectGlossiness;
uniform float objectOpacity;

uniform int isCapsule;
uniform vec3 capsulePosition;
uniform vec4 capsuleRotation;
uniform float capsuleHalfLength;
uniform float capsuleRadius;
uniform vec3 capsuleStart;
uniform vec3 capsuleVector;

uniform int shadowCapsuleCount;
uniform vec3 shadowCapsuleStarts[SHADOW_CAPSULES_MAX];
uniform vec3 shadowCapsuleVectors[SHADOW_CAPSULES_MAX];
uniform float shadowCapsuleRadii[SHADOW_CAPSULES_MAX];
uniform sampler2D shadowLookupTable;
uniform vec2 shadowLookupResolution;

uniform int aoCapsuleCount;
uniform vec3 aoCapsuleStarts[AO_CAPSULES_MAX];
uniform vec3 aoCapsuleVectors[AO_CAPSULES_MAX];
uniform float aoCapsuleRadii[AO_CAPSULES_MAX];
uniform sampler2D aoLookupTable;
uniform vec2 aoLookupResolution;

uniform vec3 cameraPosition;

uniform float sunStrength;
uniform vec3 sunDir;
uniform vec3 sunColor;
uniform float skyStrength;
uniform vec3 skyColor;
uniform float ambientStrength;
uniform float groundStrength;
uniform float exposure;

uniform int isCharacter;
uniform sampler2D characterTexture;

out vec4 finalColor;

vec3 ToGamma(in vec3 col)
{
    return vec3(pow(col.x, 2.2), pow(col.y, 2.2), pow(col.z, 2.2));
}

vec3 FromGamma(in vec3 col)
{
    return vec3(pow(col.x, 1.0/2.2), pow(col.y, 1.0/2.2), pow(col.z, 1.0/2.2));
}

float Saturate(in float x)
{
    return clamp(x, 0.0, 1.0);
}

float Square(in float x)
{
    return x * x;
}

float FastAcos(in float x)
{
    float y = abs(x);
    float p = -0.1565827 * y + 1.570796;
    p *= sqrt(max(1.0 - y, 0.0));
    return x >= 0.0 ? p : PI - p;
}

float FastPositiveAcos(in float x)
{
    float p = -0.1565827 * x + 1.570796;
    return p * sqrt(max(1.0 - x, 0.0));
}

vec3 Rotate(in vec4 q, vec3 v)
{
    vec3 t = 2.0 * cross(q.xyz, v);
    return v + q.w * t + cross(q.xyz, t);
}

vec3 Unrotate(in vec4 q, vec3 v)
{
    return Rotate(vec4(-q.x, -q.y, -q.z, q.w), v);
}

float Checker(in vec2 uv)
{
    vec4 uvDDXY = vec4(dFdx(uv), dFdy(uv));
    vec2 w = vec2(length(uvDDXY.xz), length(uvDDXY.yw));
    vec2 i = 2.0*(abs(fract((uv-0.5*w)*0.5)-0.5)-
                  abs(fract((uv+0.5*w)*0.5)-0.5))/w;
    return 0.5 - 0.5*i.x*i.y;
}

float Grid(in vec2 uv, in float lineWidth)
{
    vec4 uvDDXY = vec4(dFdx(uv), dFdy(uv));
    vec2 uvDeriv = vec2(length(uvDDXY.xz), length(uvDDXY.yw));
    float targetWidth = lineWidth > 0.5 ? 1.0 - lineWidth : lineWidth;
    vec2 drawWidth = clamp(
        vec2(targetWidth, targetWidth), uvDeriv, vec2(0.5, 0.5));
    vec2 lineAA = uvDeriv * 1.5;
    vec2 gridUV = abs(fract(uv) * 2.0 - 1.0);
    gridUV = lineWidth > 0.5 ? gridUV : 1.0 - gridUV;
    vec2 g2 = smoothstep(drawWidth + lineAA, drawWidth - lineAA, gridUV);
    g2 *= clamp(targetWidth / drawWidth, 0.0, 1.0);
    g2 = mix(g2, vec2(targetWidth, targetWidth),
        clamp(uvDeriv * 2.0 - 1.0, 0.0, 1.0));
    g2 = lineWidth > 0.5 ? 1.0 - g2 : g2;
    return mix(g2.x, 1.0, g2.y);
}

float SphereOcclusion(in vec3 pos, in vec3 nor, in vec3 sph, in float rad)
{
    vec3 di = sph - pos;
    float l = length(di);
    float nlAngle = FastAcos(dot(nor, di / l));
    float h  = l < rad ? 1.0 : l / rad;
    vec2 uvs = vec2(nlAngle / PI, (h - 1.0) / (AO_RATIO_MAX - 1.0));
    uvs = uvs * (aoLookupResolution - 1.0) / aoLookupResolution + 0.5 / aoLookupResolution;
    return texture(aoLookupTable, uvs).r;
}

float SphereDirectionalOcclusion(
    in vec3 pos, 
    in vec3 sphere, 
    in float radius,
    in vec3 coneDir)
{
    vec3 occluder = sphere - pos;
    float occluderLen2 = dot(occluder, occluder);
    vec3 occluderDir = occluder * inversesqrt(occluderLen2);
    float phi = FastAcos(dot(occluderDir, -coneDir));
    float theta = FastPositiveAcos(sqrt(occluderLen2 / (Square(radius) + occluderLen2)));

    vec2 uvs = vec2(phi / PI, theta / (PI / 2.0));
    uvs = uvs * (shadowLookupResolution - 1.0) / shadowLookupResolution + 0.5 / shadowLookupResolution;
    return texture(shadowLookupTable, uvs).r;
}

float CapsuleOcclusion(
    in vec3 pos, 
    in vec3 nor,
    in vec3 capStart, 
    in vec3 capVec, 
    in float radius)
{
    vec3 ba = capVec;
    vec3 pa = pos - capStart;
    float l = dot(ba, ba);
    float t = abs(l) < 1e-8f ? 0.0 : Saturate(dot(pa, ba) / l);
    return SphereOcclusion(pos, nor, capStart + t * ba, radius);
}

float CapsuleDirectionalOcclusion(
    in vec3 pos, in vec3 capStart, in vec3 capVec,
    in float capRadius, in vec3 coneDir)
{
    vec3 ba = capVec;
    vec3 pa = capStart - pos;
    vec3 cba = dot(-coneDir, ba) * -coneDir - ba;
    float t = Saturate(dot(pa, cba) / dot(cba, cba));

    return SphereDirectionalOcclusion(pos, capStart + t * ba, capRadius, coneDir);
}

vec2 CapsuleUVs(
    in vec3 pos, in vec3 capPos,
    in vec4 capRot, in float capHalfLength,
    in float capRadius, in vec2 scale)
{
    vec3 loc = Unrotate(capRot, pos - capPos);

    vec2 limit = vec2(
        2.0 * capHalfLength + 2.0 * capRadius,
        PI * capRadius);

    vec2 repeat = max(round(scale * limit), 1.0);

    return (repeat / limit) * vec2(loc.x, capRadius * atan(loc.z, loc.y));
}

vec3 CapsuleNormal(
    in vec3 pos, in vec3 capStart,
    in vec3 capVec)
{
    vec3 ba = capVec;
    vec3 pa = pos - capStart;
    float h = Saturate(dot(pa, ba) / dot(ba, ba));
    return normalize(pa - h*ba);
}

void main()
{
    vec3 pos = fragPosition;
    vec3 nor = fragNormal;
    vec2 uvs = fragTexCoord;

    // Recompute uvs and normals if capsule

    if (isCapsule == 1)
    {
        uvs = CapsuleUVs(
            pos,
            capsulePosition,
            capsuleRotation,
            capsuleHalfLength,
            capsuleRadius,
            vec2(4.0, 4.0));
        
        nor = CapsuleNormal(pos, capsuleStart, capsuleVector);
    }

    // Compute sun shadow amount

    float sunShadow = 1.0;
    for (int i = 0; i < shadowCapsuleCount; i++)
    {
        sunShadow = min(sunShadow, CapsuleDirectionalOcclusion(
            pos,
            shadowCapsuleStarts[i],
            shadowCapsuleVectors[i],
            shadowCapsuleRadii[i],
            sunDir));
    }
    
    // Compute ambient shadow amount

    float ambShadow = 1.0;
    for (int i = 0; i < aoCapsuleCount; i++)
    {
        ambShadow = min(ambShadow, CapsuleOcclusion(
            pos, nor,
            aoCapsuleStarts[i],
            aoCapsuleVectors[i],
            aoCapsuleRadii[i]));
    }

    // Compute albedo from grid and checker
    
    vec3 albedo;
    float specularity;
    if (isCharacter == 1)
    {
        albedo = FromGamma(vec3(texture(characterTexture, uvs)));
        specularity = 0.0;
    }
    else
    {
        float gridFine = Grid(20.0 * uvs, 0.025);
        float gridCoarse = Grid(2.0 * uvs, 0.02);
        float check = Checker(2.0 * uvs);
        albedo = FromGamma(objectColor) * mix(mix(mix(0.9, 0.95, check), 0.85, gridFine), 1.0, gridCoarse);
        specularity = objectSpecularity * mix(mix(0.0, 0.75, check), 1.0, gridCoarse);
    }

    // Compute lighting
    
    vec3 eyeDir = normalize(pos - cameraPosition);

    vec3 lightSunColor = FromGamma(sunColor);
    vec3 lightSunHalf = normalize(sunDir + eyeDir);

    vec3 lightSkyColor = FromGamma(skyColor);
    vec3 skyDir = vec3(0.0, -1.0, 0.0);
    vec3 lightSkyHalf = normalize(skyDir + eyeDir);

    float sunFactorDiff = max(dot(nor, -sunDir), 0.0);
    float sunFactorSpec = specularity *
        ((objectGlossiness+2.0) / (8.0 * PI)) *
        pow(max(dot(nor, lightSunHalf), 0.0), objectGlossiness);

    float skyFactorDiff = max(dot(nor, -skyDir), 0.0);
    float skyFactorSpec = specularity *
        ((objectGlossiness+2.0) / (8.0 * PI)) *
        pow(max(dot(nor, lightSkyHalf), 0.0), objectGlossiness);

    float groundFactorDiff = max(dot(nor, skyDir), 0.0);
    
    // Combine
    
    vec3 ambient = ambShadow * ambientStrength * lightSkyColor * albedo;

    vec3 diffuse = sunShadow * sunStrength * lightSunColor * albedo * sunFactorDiff +
        groundStrength * lightSkyColor * albedo * groundFactorDiff;
        skyStrength * lightSkyColor * albedo * skyFactorDiff;

    float specular = sunShadow * sunStrength * sunFactorSpec + skyStrength * skyFactorSpec;

    vec3 final = diffuse + ambient + specular;

    finalColor = vec4(ToGamma(exposure * final), objectOpacity); // Original
    // finalColor = vec4(nor[0], nor[1], nor[2], 1.0); // DEBUG Normals
    // finalColor = vec4(uvs[0], uvs[1], 0.0, 1.0); // DEBUG UV coords
}

);

// Structure containing all the uniform indices for the shader
typedef struct
{
    int isCapsule;
    int capsulePosition;
    int capsuleRotation;
    int capsuleHalfLength;
    int capsuleRadius;
    int capsuleStart;
    int capsuleVector;

    int shadowCapsuleCount;
    int shadowCapsuleStarts;
    int shadowCapsuleVectors;
    int shadowCapsuleRadii;
    int shadowLookupTable;
    int shadowLookupResolution;

    int aoCapsuleCount;
    int aoCapsuleStarts;
    int aoCapsuleVectors;
    int aoCapsuleRadii;
    int aoLookupTable;
    int aoLookupResolution;

    int cameraPosition;

    int objectColor;
    int objectSpecularity;
    int objectGlossiness;
    int objectOpacity;

    int sunStrength;
    int sunDir;
    int sunColor;
    int skyStrength;
    int skyColor;
    int ambientStrength;
    int groundStrength;

    int isCharacter;
    int characterTexture;

    int exposure;

} ShaderUniforms;

// Lookup all shader uniform indices
static void ShaderUniformsInit(ShaderUniforms* uniforms, Shader shader)
{
    uniforms->isCapsule = GetShaderLocation(shader, "isCapsule");
    uniforms->capsulePosition =  GetShaderLocation(shader, "capsulePosition");
    uniforms->capsuleRotation =  GetShaderLocation(shader, "capsuleRotation");
    uniforms->capsuleHalfLength =  GetShaderLocation(shader, "capsuleHalfLength");
    uniforms->capsuleRadius =  GetShaderLocation(shader, "capsuleRadius");
    uniforms->capsuleStart =  GetShaderLocation(shader, "capsuleStart");
    uniforms->capsuleVector =  GetShaderLocation(shader, "capsuleVector");

    uniforms->shadowCapsuleCount = GetShaderLocation(shader, "shadowCapsuleCount");
    uniforms->shadowCapsuleStarts =  GetShaderLocation(shader, "shadowCapsuleStarts");
    uniforms->shadowCapsuleVectors =  GetShaderLocation(shader, "shadowCapsuleVectors");
    uniforms->shadowCapsuleRadii =  GetShaderLocation(shader, "shadowCapsuleRadii");
    uniforms->shadowLookupTable =  GetShaderLocation(shader, "shadowLookupTable");
    uniforms->shadowLookupResolution =  GetShaderLocation(shader, "shadowLookupResolution");

    uniforms->aoCapsuleCount = GetShaderLocation(shader, "aoCapsuleCount");
    uniforms->aoCapsuleStarts =  GetShaderLocation(shader, "aoCapsuleStarts");
    uniforms->aoCapsuleVectors =  GetShaderLocation(shader, "aoCapsuleVectors");
    uniforms->aoCapsuleRadii =  GetShaderLocation(shader, "aoCapsuleRadii");
    uniforms->aoLookupTable =  GetShaderLocation(shader, "aoLookupTable");
    uniforms->aoLookupResolution =  GetShaderLocation(shader, "aoLookupResolution");

    uniforms->cameraPosition = GetShaderLocation(shader, "cameraPosition");

    uniforms->objectColor = GetShaderLocation(shader, "objectColor");
    uniforms->objectSpecularity = GetShaderLocation(shader, "objectSpecularity");
    uniforms->objectGlossiness = GetShaderLocation(shader, "objectGlossiness");
    uniforms->objectOpacity = GetShaderLocation(shader, "objectOpacity");

    uniforms->sunStrength = GetShaderLocation(shader, "sunStrength");
    uniforms->sunDir = GetShaderLocation(shader, "sunDir");
    uniforms->sunColor = GetShaderLocation(shader, "sunColor");
    uniforms->skyStrength = GetShaderLocation(shader, "skyStrength");
    uniforms->skyColor = GetShaderLocation(shader, "skyColor");
    uniforms->ambientStrength = GetShaderLocation(shader, "ambientStrength");
    uniforms->groundStrength = GetShaderLocation(shader, "groundStrength");

    uniforms->isCharacter = GetShaderLocation(shader, "isCharacter");
    uniforms->characterTexture =  GetShaderLocation(shader, "characterTexture");

    uniforms->exposure = GetShaderLocation(shader, "exposure");
}

//----------------------------------------------------------------------------------
// Models
//----------------------------------------------------------------------------------

// Embedded Capsule OBJ file
static const char* capsuleOBJ = "\
v 0.82165808 -0.82165808 -1.0579772e-18\nv 0.82165808 -0.58100000 0.58100000\n\
v 0.82165808 8.7595780e-17 0.82165808\nv 0.82165808 0.58100000 0.58100000\n\
v 0.82165808 0.82165808 9.9566116e-17\nv 0.82165808 0.58100000 -0.58100000\n\
v 0.82165808 2.8884397e-16 -0.82165808\nv 0.82165808 -0.58100000 -0.58100000\n\
v -0.82165808 -0.82165808 -1.0579772e-18\nv -0.82165808 -0.58100000 0.58100000\n\
v -0.82165808 -1.3028313e-17 0.82165808\nv -0.82165808 0.58100000 0.58100000\n\
v -0.82165808 0.82165808 9.9566116e-17\nv -0.82165808 0.58100000 -0.58100000\n\
v -0.82165808 1.8821987e-16 -0.82165808\nv -0.82165808 -0.58100000 -0.58100000\n\
v 1.16200000 1.5874776e-16 -1.0579772e-18\nv -1.16200000 1.6443801e-17 -1.0579772e-18\n\
v -9.1030792e-3 -1.15822938 -1.0579772e-18\nv 9.1030792e-3 -1.15822938 -1.0579772e-18\n\
v 9.1030792e-3 -0.81899185 0.81899185\nv -9.1030792e-3 -0.81899185 0.81899185\n\
v 9.1030792e-3 1.7232088e-17 1.15822938\nv -9.1030792e-3 1.6117282e-17 1.15822938\n\
v 9.1030792e-3 0.81899185 0.81899185\nv -9.1030792e-3 0.81899185 0.81899185\n\
v 9.1030792e-3 1.15822938 1.4078421e-16\nv -9.1030792e-3 1.15822938 1.4078421e-16\n\
v 9.1030792e-3 0.81899185 -0.81899185\nv -9.1030792e-3 0.81899185 -0.81899185\n\
v 9.1030792e-3 3.0091647e-16 -1.15822938\nv -9.1030792e-3 2.9980166e-16 -1.15822938\n\
v 9.1030792e-3 -0.81899185 -0.81899185\nv -9.1030792e-3 -0.81899185 -0.81899185\n\
vn 0.71524683 -0.69887193 -2.5012597e-16\nvn 0.61185516 -0.55930013 0.55930013\n\
vn 0.71524683 0.0000000e+0 0.69887193\nvn 0.61185516 0.55930013 0.55930013\n\
vn 0.71524683 0.69887193 1.5632873e-17\nvn 0.61185516 0.55930013 -0.55930013\n\
vn 0.71524683 6.2531494e-17 -0.69887193\nvn 0.61185516 -0.55930013 -0.55930013\n\
vn -0.71524683 -0.69887193 -2.5012597e-16\nvn -0.61185516 -0.55930013 0.55930013\n\
vn -0.71524683 0.0000000e+0 0.69887193\nvn -0.61185516 0.55930013 0.55930013\n\
vn -0.71524683 0.69887193 4.6898620e-17\nvn -0.61185516 0.55930013 -0.55930013\n\
vn -0.71524683 4.6898620e-17 -0.69887193\nvn -0.61185516 -0.55930013 -0.55930013\n\
vn 1.00000000 1.5208752e-17 -2.6615316e-17\nvn -1.00000000 -1.5208752e-17 2.2813128e-17\n\
vn -0.19614758 -0.98057439 -2.2848712e-16\nvn 0.26047011 -0.96548191 -2.4273177e-16\n\
vn 0.13072302 -0.70103905 0.70103905\nvn -0.19614758 -0.69337080 0.69337080\n\
vn 0.22349711 5.9825845e-2 0.97286685\nvn -0.22349711 -5.9825845e-2 0.97286685\n\
vn 0.15641931 0.75510180 0.63667438\nvn -0.15641931 0.63667438 0.75510180\n\
vn 0.22349711 0.97286685 -5.9825845e-2\nvn -0.22349711 0.97286685 5.9825845e-2\n\
vn 0.15641931 0.63667438 -0.75510180\nvn -0.15641931 0.75510180 -0.63667438\n\
vn 0.22349711 -5.9825845e-2 -0.97286685\nvn -0.22349711 5.9825845e-2 -0.97286685\n\
vn 0.15641931 -0.75510180 -0.63667438\nvn -0.15641931 -0.63667438 -0.75510180\n\
f 1//1 17//17 2//2\nf 1//1 20//20 8//8\nf 2//2 17//17 3//3\nf 2//2 20//20 1//1\n\
f 2//2 23//23 21//21\nf 3//3 17//17 4//4\nf 3//3 23//23 2//2\nf 4//4 17//17 5//5\n\
f 4//4 23//23 3//3\nf 4//4 27//27 25//25\nf 5//5 17//17 6//6\nf 5//5 27//27 4//4\n\
f 6//6 17//17 7//7\nf 6//6 27//27 5//5\nf 6//6 31//31 29//29\nf 7//7 17//17 8//8\n\
f 7//7 31//31 6//6\nf 8//8 17//17 1//1\nf 8//8 20//20 33//33\nf 8//8 31//31 7//7\n\
f 9//9 18//18 16//16\nf 9//9 19//19 10//10\nf 10//10 18//18 9//9\nf 10//10 19//19 22//22\n\
f 10//10 24//24 11//11\nf 11//11 18//18 10//10\nf 11//11 24//24 12//12\nf 12//12 18//18 11//11\n\
f 12//12 24//24 26//26\nf 12//12 28//28 13//13\nf 13//13 18//18 12//12\nf 13//13 28//28 14//14\n\
f 14//14 18//18 13//13\nf 14//14 28//28 30//30\nf 14//14 32//32 15//15\nf 15//15 18//18 14//14\n\
f 15//15 32//32 16//16\nf 16//16 18//18 15//15\nf 16//16 19//19 9//9\nf 16//16 32//32 34//34\n\
f 19//19 33//33 20//20\nf 20//20 21//21 19//19\nf 21//21 20//20 2//2\nf 21//21 24//24 22//22\n\
f 22//22 19//19 21//21\nf 22//22 24//24 10//10\nf 23//23 26//26 24//24\nf 24//24 21//21 23//23\n\
f 25//25 23//23 4//4\nf 25//25 28//28 26//26\nf 26//26 23//23 25//25\nf 26//26 28//28 12//12\n\
f 27//27 30//30 28//28\nf 28//28 25//25 27//27\nf 29//29 27//27 6//6\nf 29//29 32//32 30//30\n\
f 30//30 27//27 29//29\nf 30//30 32//32 14//14\nf 31//31 34//34 32//32\nf 32//32 29//29 31//31\n\
f 33//33 19//19 34//34\nf 33//33 31//31 8//8\nf 34//34 19//19 16//16\nf 34//34 31//31 33//33";

#undef TINYOBJ_LOADER_C_IMPLEMENTATION
#include "external/tinyobj_loader_c.h"

// Extra function for loading OBJ from memory
static Model LoadOBJFromMemory(const char *fileText)
{
    Model model = { 0 };

    tinyobj_attrib_t attrib = { 0 };
    tinyobj_shape_t *meshes = NULL;
    unsigned int meshCount = 0;

    tinyobj_material_t *materials = NULL;
    unsigned int materialCount = 0;

    if (fileText != NULL)
    {
        unsigned int dataSize = (unsigned int)strlen(fileText);

        unsigned int flags = TINYOBJ_FLAG_TRIANGULATE;
        tinyobj_parse_obj(&attrib, &meshes, &meshCount, &materials, &materialCount, fileText, dataSize, flags);
        
        model.meshCount = 1;
        model.meshes = (Mesh *)RL_CALLOC(model.meshCount, sizeof(Mesh));
        model.meshMaterial = (int *)RL_CALLOC(model.meshCount, sizeof(int));

        // Count the faces for each material
        int *matFaces = (int *)RL_CALLOC(model.meshCount, sizeof(int));
        matFaces[0] = attrib.num_faces;

        //--------------------------------------
        // Create the material meshes

        // Running counts/indexes for each material mesh as we are
        // building them at the same time
        int *vCount = (int *)RL_CALLOC(model.meshCount, sizeof(int));
        int *vtCount = (int *)RL_CALLOC(model.meshCount, sizeof(int));
        int *vnCount = (int *)RL_CALLOC(model.meshCount, sizeof(int));
        int *faceCount = (int *)RL_CALLOC(model.meshCount, sizeof(int));

        // Allocate space for each of the material meshes
        for (int mi = 0; mi < model.meshCount; mi++)
        {
            model.meshes[mi].vertexCount = matFaces[mi]*3;
            model.meshes[mi].triangleCount = matFaces[mi];
            model.meshes[mi].vertices = (float *)RL_CALLOC(model.meshes[mi].vertexCount*3, sizeof(float));
            model.meshes[mi].texcoords = (float *)RL_CALLOC(model.meshes[mi].vertexCount*2, sizeof(float));
            model.meshes[mi].normals = (float *)RL_CALLOC(model.meshes[mi].vertexCount*3, sizeof(float));
            model.meshMaterial[mi] = mi;
        }

        // Scan through the combined sub meshes and pick out each material mesh
        for (unsigned int af = 0; af < attrib.num_faces; af++)
        {
            int mm = attrib.material_ids[af];   // mesh material for this face
            if (mm == -1) { mm = 0; }           // no material object..

            // Get indices for the face
            tinyobj_vertex_index_t idx0 = attrib.faces[3*af + 0];
            tinyobj_vertex_index_t idx1 = attrib.faces[3*af + 1];
            tinyobj_vertex_index_t idx2 = attrib.faces[3*af + 2];

            // Fill vertices buffer (float) using vertex index of the face
            for (int v = 0; v < 3; v++) { model.meshes[mm].vertices[vCount[mm] + v] = attrib.vertices[idx0.v_idx*3 + v]; } vCount[mm] +=3;
            for (int v = 0; v < 3; v++) { model.meshes[mm].vertices[vCount[mm] + v] = attrib.vertices[idx1.v_idx*3 + v]; } vCount[mm] +=3;
            for (int v = 0; v < 3; v++) { model.meshes[mm].vertices[vCount[mm] + v] = attrib.vertices[idx2.v_idx*3 + v]; } vCount[mm] +=3;

            if (attrib.num_texcoords > 0)
            {
                // Fill texcoords buffer (float) using vertex index of the face
                // NOTE: Y-coordinate must be flipped upside-down to account for
                // raylib's upside down textures...
                model.meshes[mm].texcoords[vtCount[mm] + 0] = attrib.texcoords[idx0.vt_idx*2 + 0];
                model.meshes[mm].texcoords[vtCount[mm] + 1] = 1.0f - attrib.texcoords[idx0.vt_idx*2 + 1]; vtCount[mm] += 2;
                model.meshes[mm].texcoords[vtCount[mm] + 0] = attrib.texcoords[idx1.vt_idx*2 + 0];
                model.meshes[mm].texcoords[vtCount[mm] + 1] = 1.0f - attrib.texcoords[idx1.vt_idx*2 + 1]; vtCount[mm] += 2;
                model.meshes[mm].texcoords[vtCount[mm] + 0] = attrib.texcoords[idx2.vt_idx*2 + 0];
                model.meshes[mm].texcoords[vtCount[mm] + 1] = 1.0f - attrib.texcoords[idx2.vt_idx*2 + 1]; vtCount[mm] += 2;
            }

            if (attrib.num_normals > 0)
            {
                // Fill normals buffer (float) using vertex index of the face
                for (int v = 0; v < 3; v++) { model.meshes[mm].normals[vnCount[mm] + v] = attrib.normals[idx0.vn_idx*3 + v]; } vnCount[mm] +=3;
                for (int v = 0; v < 3; v++) { model.meshes[mm].normals[vnCount[mm] + v] = attrib.normals[idx1.vn_idx*3 + v]; } vnCount[mm] +=3;
                for (int v = 0; v < 3; v++) { model.meshes[mm].normals[vnCount[mm] + v] = attrib.normals[idx2.vn_idx*3 + v]; } vnCount[mm] +=3;
            }
        }

        model.materialCount = 1;
        model.materials = (Material *)RL_CALLOC(model.materialCount, sizeof(Material));
        model.materials[0] = LoadMaterialDefault();

        tinyobj_attrib_free(&attrib);
        tinyobj_shapes_free(meshes, meshCount);
        tinyobj_materials_free(materials, materialCount);

        RL_FREE(matFaces);
        RL_FREE(vCount);
        RL_FREE(vtCount);
        RL_FREE(vnCount);
        RL_FREE(faceCount);
    }

    // Make sure model transform is set to identity matrix!
    model.transform = MatrixIdentity();

    // Upload vertex data to GPU (static mesh)
    for (int i = 0; i < model.meshCount; i++) { UploadMesh(&model.meshes[i], false); }

    return model;
}

//----------------------------------------------------------------------------------
// Render Settings
//----------------------------------------------------------------------------------

typedef struct {

    Color backgroundColor;

    float sunLightConeAngle;
    float sunLightStrength;
    float sunAzimuth;
    float sunAltitude;
    Color sunColor;

    float skyLightStrength;
    Color skyColor;

    float groundLightStrength;
    float ambientLightStrength;

    float exposure;

    bool drawOrigin;
    bool drawGrid;
    bool drawChecker;
    bool drawCapsules;
    bool drawWireframes;
    bool drawSkeleton;
    bool drawTransforms;
    bool drawAO;
    bool drawShadows;
    bool drawEndSites;
    bool drawFPS;
    bool drawUI;

} RenderSettings;

void RenderSettingsInit(RenderSettings* settings, int argc, char** argv)
{
    settings->backgroundColor = ArgColor(argc, argv, "backgroundColor", WHITE);

    settings->sunLightConeAngle = ArgFloat(argc, argv, "sunLightConeAngle", 0.2f);
    settings->sunLightStrength = ArgFloat(argc, argv, "sunLightStrength", 0.25f);
    settings->sunAzimuth = ArgFloat(argc, argv, "sunAzimuth", PI / 4.0f);
    settings->sunAltitude = ArgFloat(argc, argv, "sunAltitude", 0.8f);
    settings->sunColor = ArgColor(argc, argv, "sunColor", (Color){ 253, 255, 232 });

    settings->skyLightStrength = ArgFloat(argc, argv, "skyLightStrength", 0.2f);
    settings->skyColor = ArgColor(argc, argv, "skyColor", (Color){ 174, 183, 190 });

    settings->groundLightStrength = ArgFloat(argc, argv, "groundLightStrength", 0.1f);
    settings->ambientLightStrength = ArgFloat(argc, argv, "ambientLightStrength", 1.0f);

    settings->exposure = ArgFloat(argc, argv, "exposure", 0.9f);

    settings->drawOrigin = ArgBool(argc, argv, "drawOrigin", true);
    settings->drawGrid = ArgBool(argc, argv, "drawGrid", false);
    settings->drawChecker = ArgBool(argc, argv, "drawChecker", true);
    settings->drawCapsules = ArgBool(argc, argv, "drawCapsules", true);
    settings->drawWireframes = ArgBool(argc, argv, "drawWireframes", false);
    settings->drawSkeleton = ArgBool(argc, argv, "drawSkeleton", true);
    settings->drawTransforms = ArgBool(argc, argv, "drawTransforms", false);
    settings->drawAO = ArgBool(argc, argv, "drawAO", true);
    settings->drawShadows = ArgBool(argc, argv, "drawShadows", true);
    settings->drawEndSites = ArgBool(argc, argv, "drawEndSites", true);
    settings->drawFPS = ArgBool(argc, argv, "drawFPS", false);
    settings->drawUI = ArgBool(argc, argv, "drawUI", true);
}

//--------------------------------------
// Scrubber
//--------------------------------------

typedef struct {

    bool playing;
    bool looping;
    bool inplace;
    float playTime;
    float playSpeed;
    bool frameSnap;
    int sampleMode;

    float timeLimit;
    int frameLimit;
    int frameMin;
    int frameMax;
    int frameMinSelect;
    int frameMaxSelect;
    bool frameMinEdit;
    bool frameMaxEdit;
    float timeMin;
    float timeMax;

} ScrubberSettings;

static inline void ScrubberSettingsInit(ScrubberSettings* settings, int argc, char** argv)
{
    settings->playing = ArgBool(argc, argv, "playing", true);
    settings->looping = ArgBool(argc, argv, "looping", false);
    settings->inplace = ArgBool(argc, argv, "inplace", false);
    settings->playTime = ArgFloat(argc, argv, "playTime", 0.0f);
    settings->playSpeed = ArgFloat(argc, argv, "playSpeed", 1.0f);
    settings->frameSnap = ArgBool(argc, argv, "frameSnap", true);
    settings->sampleMode = ArgEnum(argc, argv, "sampleMode", 3, (const char*[]){ "nearest", "linear", "cubic" }, 1);

    settings->timeLimit = 0.0f;
    settings->frameLimit = 0;
    settings->frameMin = 0;
    settings->frameMax = 0;
    settings->frameMinSelect = 0;
    settings->frameMaxSelect = 0;
    settings->frameMinEdit = false;
    settings->frameMaxEdit = false;
    settings->timeMin = 0.0f;
    settings->timeMax = 0.0f;
}

static inline void ScrubberSettingsRecomputeLimits(ScrubberSettings* settings, CharacterData* characterData)
{
    settings->frameLimit = 0;
    settings->timeLimit = 0.0f;
    for (int i = 0; i < characterData->count; i++)
    {
        settings->frameLimit = MaxInt(settings->frameLimit, characterData->bvhData[i].frameCount - 1);
        settings->timeLimit = Max(settings->timeLimit, (characterData->bvhData[i].frameCount - 1) * characterData->bvhData[i].frameTime);
    }
}

static inline void ScrubberSettingsInitMaxs(ScrubberSettings* settings, CharacterData* characterData)
{
    if (characterData->count == 0) { return; }

    settings->frameMax = characterData->bvhData[characterData->active].frameCount - 1;
    settings->frameMaxSelect = settings->frameMax;
    settings->timeMax = settings->frameMax * characterData->bvhData[characterData->active].frameTime;

    settings->frameMin = 0;
    settings->frameMinSelect = settings->frameMin;
    settings->timeMin = 0.0f;
}

static inline void ScrubberSettingsClamp(ScrubberSettings* settings, CharacterData* characterData)
{
    if (characterData->count == 0) { return; }

    settings->frameMax = ClampInt(settings->frameMax, 0, settings->frameLimit);
    settings->frameMaxSelect = settings->frameMax;
    settings->timeMax = settings->frameMax * characterData->bvhData[characterData->active].frameTime;

    settings->frameMin = ClampInt(settings->frameMin, 0, settings->frameMax);
    settings->frameMinSelect = settings->frameMin;
    settings->timeMin = settings->frameMin * characterData->bvhData[characterData->active].frameTime;

    settings->playTime = Clamp(settings->playTime, settings->timeMin, settings->timeMax);
}

//----------------------------------------------------------------------------------
// Drawing
//----------------------------------------------------------------------------------

static inline void DrawTransform(const Vector3 position, const Quaternion rotation, const float size)
{
    DrawLine3D(position, Vector3Add(position, Vector3RotateByQuaternion((Vector3){ size, 0.0, 0.0 }, rotation)), RED);
    DrawLine3D(position, Vector3Add(position, Vector3RotateByQuaternion((Vector3){ 0.0, size, 0.0 }, rotation)), GREEN);
    DrawLine3D(position, Vector3Add(position, Vector3RotateByQuaternion((Vector3){ 0.0, 0.0, size }, rotation)), BLUE);
}

static inline void DrawSkeleton(TransformData* xformData, bool drawEndSites, Color color, Color endSiteColor)
{
    for (int i = 0; i < xformData->jointCount; i++)
    {
        if (!xformData->endSite[i])
        {
            DrawSphereWires(
                xformData->globalPositions[i],
                0.01f,
                4,
                6,
                color);
        }
        else if (drawEndSites)
        {
            DrawCubeWiresV(
                xformData->globalPositions[i],
                (Vector3){ 0.02f, 0.02f, 0.02f },
                endSiteColor);
        }

        if (xformData->parents[i] != -1)
        {
            if (!xformData->endSite[i])
            {
                DrawLine3D(
                    xformData->globalPositions[i],
                    xformData->globalPositions[xformData->parents[i]],
                    color);
            }
            else if (drawEndSites)
            {
                DrawLine3D(
                    xformData->globalPositions[i],
                    xformData->globalPositions[xformData->parents[i]],
                    endSiteColor);
            }
        }
    }
}

static inline void DrawTransforms(TransformData* xformData)
{
    for (int i = 0; i < xformData->jointCount; i++)
    {
        if (!xformData->endSite[i])
        {
            DrawTransform(
                xformData->globalPositions[i],
                xformData->globalRotations[i],
                0.1f);
        }
    }
}

static inline void DrawWireFrames(CapsuleData* capsuleData, Color color)
{
    for (int i = 0; i < capsuleData->capsuleCount; i++)
    {
        Vector3 capsuleStart = CapsuleStart(capsuleData->capsulePositions[i], capsuleData->capsuleRotations[i], capsuleData->capsuleHalfLengths[i]);
        Vector3 capsuleEnd = CapsuleEnd(capsuleData->capsulePositions[i], capsuleData->capsuleRotations[i], capsuleData->capsuleHalfLengths[i]);
        float capsuleRadius = capsuleData->capsuleRadii[i];

        DrawSphereWires(capsuleStart, capsuleRadius, 4, 6, color);
        DrawSphereWires(capsuleEnd, capsuleRadius, 4, 6, color);
        DrawCylinderWiresEx(capsuleStart, capsuleEnd, capsuleRadius, capsuleRadius, 6, color);
    }
}

//----------------------------------------------------------------------------------
// GUI
//----------------------------------------------------------------------------------

static inline void GuiOrbitCamera(OrbitCamera* camera, CharacterData* characterData, int argc, char** argv)
{
    GuiGroupBox((Rectangle){ 20, 10, 190, 260 }, "Camera");

    GuiLabel((Rectangle){ 30, 20, 150, 20 }, "Ctrl + Left Click - Rotate");
    GuiLabel((Rectangle){ 30, 40, 150, 20 }, "Ctrl + Right Click - Pan");
    GuiLabel((Rectangle){ 30, 60, 150, 20 }, "Mouse Scroll - Zoom");
    GuiLabel((Rectangle){ 30, 80, 150, 20 }, TextFormat("Target: [% 5.3f % 5.3f % 5.3f]", camera->cam3d.target.x, camera->cam3d.target.y, camera->cam3d.target.z));
    GuiLabel((Rectangle){ 30, 100, 150, 20 }, TextFormat("Offset: [% 5.3f % 5.3f % 5.3f]", camera->offset.x, camera->offset.y, camera->offset.z));
    GuiLabel((Rectangle){ 30, 120, 150, 20 }, TextFormat("Azimuth: %5.3f", camera->azimuth));
    GuiLabel((Rectangle){ 30, 140, 150, 20 }, TextFormat("Altitude: %5.3f", camera->altitude));
    GuiLabel((Rectangle){ 30, 160, 150, 20 }, TextFormat("Distance: %5.3f", camera->distance));
    
    if (GuiButton((Rectangle){ 30, 180, 100, 20 }, "Reset"))
    {
        camera->azimuth = ArgFloat(argc, argv, "cameraAzimuth", 0.0f);
        camera->altitude = ArgFloat(argc, argv, "cameraAltitude", 0.4f);
        camera->distance = ArgFloat(argc, argv, "cameraDistance", 4.0f);
        camera->offset = ArgVector3(argc, argv, "cameraOffset", Vector3Zero());
        camera->track = ArgBool(argc, argv, "cameraTrack", true);
        camera->trackBone = ArgInt(argc, argv, "cameraTrackBone", 0);
    }

    if (characterData->count > 0)
    {
        GuiToggle((Rectangle){ 30, 210, 100, 20 }, "Track", &camera->track);
        GuiComboBox((Rectangle){ 30, 240, 150, 20 }, characterData->jointNamesCombo[characterData->active], &camera->trackBone);
    }
}

static inline void GuiRenderSettings(RenderSettings* settings, CapsuleData* capsuleData, int screenWidth, int screenHeight)
{
    GuiGroupBox((Rectangle){ screenWidth - 260, 10, 240, 430 }, "Rendering");

    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 20, 100, 20 },
        "Exposure",
        TextFormat("%5.2f", settings->exposure),
        &settings->exposure,
        0.0f, 3.0f);

    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 50, 100, 20 },
        "Sun Light",
        TextFormat("%5.2f", settings->sunLightStrength),
        &settings->sunLightStrength,
        0.0f, 1.0f);

    if (GuiSliderBar(
        (Rectangle){ screenWidth - 160, 80, 100, 20 },
        "Sun Softness",
        TextFormat("%5.2f", settings->sunLightConeAngle),
        &settings->sunLightConeAngle,
        0.02f, PI / 4.0f))
    {
        CapsuleDataUpdateShadowLookupTable(capsuleData, settings->sunLightConeAngle);
    }

    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 110, 100, 20 },
        "Sky Light",
        TextFormat("%5.2f", settings->skyLightStrength),
        &settings->skyLightStrength,
        0.0f, 1.0f);

    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 140, 100, 20 },
        "Ambient Light",
        TextFormat("%5.2f", settings->ambientLightStrength),
        &settings->ambientLightStrength,
        0.0f, 2.0f);

    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 170, 100, 20 },
        "Ground Light",
        TextFormat("%5.2f", settings->groundLightStrength),
        &settings->groundLightStrength,
        0.0f, 0.5f);

    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 200, 100, 20 },
        "Sun Azimuth",
        TextFormat("%5.2f", settings->sunAzimuth),
        &settings->sunAzimuth,
        -PI, PI);

    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 230, 100, 20 },
        "Sun Altitude",
        TextFormat("%5.2f", settings->sunAltitude),
        &settings->sunAltitude,
        0.0f, 0.49f * PI);

    GuiCheckBox((Rectangle){ screenWidth - 250, 260, 20, 20 }, "Draw Origin", &settings->drawOrigin);
    GuiCheckBox((Rectangle){ screenWidth - 130, 260, 20, 20 }, "Draw Grid", &settings->drawGrid);
    GuiCheckBox((Rectangle){ screenWidth - 250, 290, 20, 20 }, "Draw Checker", &settings->drawChecker);
    GuiCheckBox((Rectangle){ screenWidth - 130, 290, 20, 20 }, "Draw Capsules", &settings->drawCapsules);
    GuiCheckBox((Rectangle){ screenWidth - 250, 320, 20, 20 }, "Draw Wireframes", &settings->drawWireframes);
    GuiCheckBox((Rectangle){ screenWidth - 130, 320, 20, 20 }, "Draw Skeleton", &settings->drawSkeleton);
    GuiCheckBox((Rectangle){ screenWidth - 250, 350, 20, 20 }, "Draw Transforms", &settings->drawTransforms);
    GuiCheckBox((Rectangle){ screenWidth - 130, 350, 20, 20 }, "Draw AO", &settings->drawAO);
    GuiCheckBox((Rectangle){ screenWidth - 250, 380, 20, 20 }, "Draw Shadows", &settings->drawShadows);
    GuiCheckBox((Rectangle){ screenWidth - 130, 380, 20, 20 }, "Draw End Sites", &settings->drawEndSites);
    GuiCheckBox((Rectangle){ screenWidth - 250, 410, 20, 20 }, "Draw FPS", &settings->drawFPS);
    GuiLabel((Rectangle){ screenWidth - 130, 410, 100, 20 }, "H Key - Hide UI");
}

static inline void GuiCharacterData(
    CharacterData* characterData,
    GuiWindowFileDialogState* fileDialogState,
    ScrubberSettings* scrubberSettings,
    char* errMsg,
    int argc,
    char** argv)
{
    int offsetHeight = 280;
  
    GuiGroupBox((Rectangle){ 20, offsetHeight, 190, (CHARACTERS_MAX - 1) * 30 + 150 }, "Characters");

#if !defined(PLATFORM_WEB)
    if (GuiButton((Rectangle){ 30, offsetHeight + 10, 110, 20 }, "Open"))
    {
        fileDialogState->windowActive = true;
    }
#endif

    if (GuiButton((Rectangle){ 150, offsetHeight + 10, 50, 20 }, "Clear"))
    {
        characterData->count = 0;
        errMsg[0] = '\0';
        ScrubberSettingsInit(scrubberSettings, argc, argv);
        SetWindowTitle("BVHView");
   }

    for (int i = 0; i < characterData->count; i++)
    {
        char bvhNameShort[20];
        bvhNameShort[0] = '\0';
        if (strlen(characterData->names[i]) + 1 <= 20)
        {
            strcat(bvhNameShort, characterData->names[i]);
        }
        else
        {
            memcpy(bvhNameShort, characterData->names[i], 16);
            memcpy(bvhNameShort + 16, "...", 4);
        }

        bool bvhSelected = i == characterData->active;
        GuiToggle((Rectangle){ 30, offsetHeight + 40 + i * 30, 140, 20 }, bvhNameShort, &bvhSelected);

        if (bvhSelected && (characterData->active != i))
        {
            characterData->active = i;
            ScrubberSettingsClamp(scrubberSettings, characterData);
            
            char windowTitle[512];
            snprintf(windowTitle, 512, "%s - BVHView", characterData->filePaths[characterData->active]);
            SetWindowTitle(windowTitle);
        }

        DrawRectangleRec((Rectangle){ 180, offsetHeight + 40 + i * 30, 20, 20 }, characterData->colors[i]);
        DrawRectangleLinesEx((Rectangle){ 180, offsetHeight + 40 + i * 30, 20, 20 }, 1, GRAY);

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT))
        {
            Vector2 mousePosition = GetMousePosition();
            if (mousePosition.x > 180 && mousePosition.x < 200 &&
                mousePosition.y > offsetHeight + 40 + i * 30 && mousePosition.y < offsetHeight + 40 + i * 30 + 20)
            {
                characterData->colorPickerActive = !characterData->colorPickerActive;
            }
        }
    }
    
    if (characterData->count > 0)
    {
        bool scaleM = characterData->scales[characterData->active] == 1.0f;
        GuiToggle((Rectangle){ 30, offsetHeight + 60 + (CHARACTERS_MAX - 1) * 30, 30, 20 }, "m", &scaleM);
        if (scaleM) { characterData->scales[characterData->active] = 1.0f; }

        bool scaleCM = characterData->scales[characterData->active] == 0.01f;
        GuiToggle((Rectangle){ 65, offsetHeight + 60 + (CHARACTERS_MAX - 1) * 30, 30, 20 }, "cm", &scaleCM);
        if (scaleCM) { characterData->scales[characterData->active] = 0.01f; }

        bool scaleInches = characterData->scales[characterData->active] == 0.0254f;
        GuiToggle((Rectangle){ 100, offsetHeight + 60 + (CHARACTERS_MAX - 1) * 30, 30, 20 }, "inch", &scaleInches);
        if (scaleInches) { characterData->scales[characterData->active] = 0.0254f; }

        bool scaleFeet = characterData->scales[characterData->active] == 0.3048f;
        GuiToggle((Rectangle){ 135, offsetHeight + 60 + (CHARACTERS_MAX - 1) * 30, 30, 20 }, "feet", &scaleFeet);
        if (scaleFeet) { characterData->scales[characterData->active] = 0.3048f; }

        bool scaleAuto = characterData->scales[characterData->active] == characterData->autoScales[characterData->active];
        GuiToggle((Rectangle){ 170, offsetHeight + 60 + (CHARACTERS_MAX - 1) * 30, 30, 20 }, "auto", &scaleAuto);
        if (scaleAuto) { characterData->scales[characterData->active] = characterData->autoScales[characterData->active]; }

        GuiSliderBar(
            (Rectangle){ 70, offsetHeight + 90 + (CHARACTERS_MAX - 1) * 30, 100, 20 },
            "Radius",
            TextFormat("%5.2f", characterData->radii[characterData->active]),
            &characterData->radii[characterData->active],
            0.01f, 0.1f);

        GuiSliderBar(
            (Rectangle){ 70, offsetHeight + 120 + (CHARACTERS_MAX - 1) * 30, 100, 20 },
            "Opacity",
            TextFormat("%5.2f", characterData->opacities[characterData->active]),
            &characterData->opacities[characterData->active],
            0.0f, 1.0f);
    }
}

static inline void GuiScrubberSettings(
    ScrubberSettings* settings,
    CharacterData* characterData,
    int screenWidth,
    int screenHeight)
{
    if (characterData->count == 0) { return; }

    float frameTime = characterData->bvhData[characterData->active].frameTime;

    GuiGroupBox((Rectangle){ screenWidth / 2 - 600, screenHeight - 100, 1200, 90 }, "Scrubber");

    GuiLabel((Rectangle){ screenWidth / 2 - 480, screenHeight - 80, 150, 20 }, TextFormat("Frame Time: %f", frameTime));
    GuiCheckBox((Rectangle){ screenWidth / 2 - 350, screenHeight - 80, 20, 20 }, "Snap to Frame", &settings->frameSnap);
    GuiComboBox((Rectangle){ screenWidth / 2 - 240, screenHeight - 80, 100, 20 }, "Nearest;Linear;Cubic", &settings->sampleMode);

    GuiToggle((Rectangle){ screenWidth / 2 - 130, screenHeight - 80, 50, 20 }, "Inplace", &settings->inplace);
    GuiToggle((Rectangle){ screenWidth / 2 - 70, screenHeight - 80, 50, 20 }, "Loop", &settings->looping);
    GuiToggle((Rectangle){ screenWidth / 2 - 10, screenHeight - 80, 50, 20 }, "Play", &settings->playing);

    bool speed01x = settings->playSpeed == 0.1f;
    GuiToggle((Rectangle){ screenWidth / 2 + 50, screenHeight - 80, 30, 20 }, "0.1x", &speed01x); if (speed01x) { settings->playSpeed = 0.1f; }
    bool speed05x = settings->playSpeed == 0.5f;
    GuiToggle((Rectangle){ screenWidth / 2 + 90, screenHeight - 80, 30, 20 }, "0.5x", &speed05x); if (speed05x) { settings->playSpeed = 0.5f; }
    bool speed1x = settings->playSpeed == 1.0f;
    GuiToggle((Rectangle){ screenWidth / 2 + 130, screenHeight - 80, 30, 20 }, "1x", &speed1x); if (speed1x) { settings->playSpeed = 1.0f; }
    bool speed2x = settings->playSpeed == 2.0f;
    GuiToggle((Rectangle){ screenWidth / 2 + 170, screenHeight - 80, 30, 20 }, "2x", &speed2x); if (speed2x) { settings->playSpeed = 2.0f; }
    bool speed4x = settings->playSpeed == 4.0f;
    GuiToggle((Rectangle){ screenWidth / 2 + 210, screenHeight - 80, 30, 20 }, "4x", &speed4x); if (speed4x) { settings->playSpeed = 4.0f; }
    GuiSliderBar((Rectangle){ screenWidth / 2 + 250, screenHeight - 80, 70, 20 }, "", TextFormat("%5.2fx", settings->playSpeed), &settings->playSpeed, 0.0f, 4.0f);

    int frame = ClampInt((int)(settings->playTime / frameTime + 0.5f), settings->frameMin, settings->frameMax);

    if (GuiValueBox(
        (Rectangle){ screenWidth / 2 - 540, screenHeight - 80, 50, 20 },
        "Min   ", &settings->frameMinSelect, 0, settings->frameLimit, settings->frameMinEdit))
    {
        settings->frameMinEdit = !settings->frameMinEdit;
        if (!settings->frameMinEdit)
        {
            settings->frameMin = settings->frameMinSelect;
            ScrubberSettingsClamp(settings, characterData);
        }
    }

    if (GuiValueBox(
        (Rectangle){ screenWidth / 2 + 470, screenHeight - 80, 50, 20 },
        "Max   ", &settings->frameMaxSelect, 0, settings->frameLimit, settings->frameMaxEdit))
    {
        settings->frameMaxEdit = !settings->frameMaxEdit;

        if (!settings->frameMaxEdit)
        {
            settings->frameMax = settings->frameMaxSelect;
            ScrubberSettingsClamp(settings, characterData);
        }
    }

    GuiLabel(
        (Rectangle){ screenWidth / 2 + 530, screenHeight - 80, 100, 20 },
        TextFormat("of %i", settings->frameLimit));

    float frameFloatPrev = settings->frameSnap ? (float)frame : settings->playTime / frameTime;
    float frameFloat = frameFloatPrev;

    GuiSliderBar(
        (Rectangle){ screenWidth / 2 - 540, screenHeight - 50, 1080, 20 },
        TextFormat("%5.2f", settings->playTime),
        TextFormat("%i", frame),
        &frameFloat,
        (float)settings->frameMin, (float)settings->frameMax);

    if (frameFloat != frameFloatPrev)
    {
        if (settings->frameSnap)
        {
            frame = ClampInt((int)(frameFloat + 0.5f), settings->frameMin, settings->frameMax);
            settings->playTime = Clamp(frame * frameTime, settings->timeMin, settings->timeMax);
        }
        else
        {
            settings->playTime = Clamp(frameFloat * frameTime, settings->timeMin, settings->timeMax);
        }
    }
}

//----------------------------------------------------------------------------------
// Application
//----------------------------------------------------------------------------------

// Structure containing all of the application state which we can then pass to the Update function
typedef struct {

    int argc;
    char** argv;

    int screenWidth;
    int screenHeight;

    OrbitCamera camera;

    Shader shader;
    ShaderUniforms uniforms;

    Mesh groundPlaneMesh;
    Model groundPlaneModel;
    Model capsuleModel;
    CharacterModel characterModel;

    CharacterData characterData;
    CapsuleData capsuleData;

    ScrubberSettings scrubberSettings;
    RenderSettings renderSettings;

    GuiWindowFileDialogState fileDialogState;

    char errMsg[512];

} ApplicationState;

// Update function - what is called to "tick" the application.
static void ApplicationUpdate(void* voidApplicationState)
{
    ApplicationState* app = voidApplicationState;

    // Process File Dialog

    if (app->fileDialogState.SelectFilePressed)
    {
        if (IsFileExtension(app->fileDialogState.fileNameText, ".bvh"))
        {
            char fileNameToLoad[512];
            snprintf(fileNameToLoad, 512, "%s/%s", app->fileDialogState.dirPathText, app->fileDialogState.fileNameText);

            if (CharacterDataLoadFromFile(&app->characterData, fileNameToLoad, app->errMsg, 512))
            {
                app->characterData.active = app->characterData.count - 1;

                CapsuleDataUpdateForCharacters(&app->capsuleData, &app->characterData);
                ScrubberSettingsRecomputeLimits(&app->scrubberSettings, &app->characterData);
                ScrubberSettingsInitMaxs(&app->scrubberSettings, &app->characterData);
                
                char windowTitle[512];
                snprintf(windowTitle, 512, "%s - BVHView", app->characterData.filePaths[app->characterData.active]);
                SetWindowTitle(windowTitle);
            }
        }
        else
        {
            snprintf(app->errMsg, 512, "Error: File '%s' is not a BVH file.", app->fileDialogState.fileNameText);
        }

        app->fileDialogState.SelectFilePressed = false;
    }

    // Process Dragged and Dropped Files

    if (IsFileDropped())
    {
        FilePathList droppedFiles = LoadDroppedFiles();

        int prevBvhCount = app->characterData.count;

        for (int i = 0; i < droppedFiles.count; i++)
        {
            if (CharacterDataLoadFromFile(&app->characterData, droppedFiles.paths[i], app->errMsg, 512))
            {
                app->characterData.active = app->characterData.count - 1;

                // convert from BVH to RayLib anim data
                app->characterData.animData[app->characterData.active] = BVHToModelAnimation(&app->characterData.bvhData[app->characterData.active], &app->characterModel.model);
            }
        }

        UnloadDroppedFiles(droppedFiles);

        if (app->characterData.count > prevBvhCount)
        {
            CapsuleDataUpdateForCharacters(&app->capsuleData, &app->characterData);
            ScrubberSettingsRecomputeLimits(&app->scrubberSettings, &app->characterData);
            ScrubberSettingsInitMaxs(&app->scrubberSettings, &app->characterData);

            char windowTitle[512];
            snprintf(windowTitle, 512, "%s - BVHView", app->characterData.filePaths[app->characterData.active]);
            SetWindowTitle(windowTitle);
        }
    }

    // Process Key Presses

    if (IsKeyPressed(KEY_H) && !app->fileDialogState.windowActive)
    {
        app->renderSettings.drawUI = !app->renderSettings.drawUI;
    }

    PROFILE_BEGIN(Update);

    // Tick time forward

    if (app->scrubberSettings.playing)
    {
        app->scrubberSettings.playTime += app->scrubberSettings.playSpeed * GetFrameTime();

        if (app->scrubberSettings.playTime >= app->scrubberSettings.timeMax)
        {
            app->scrubberSettings.playTime = (app->scrubberSettings.looping && app->scrubberSettings.timeMax >= 1e-8f) ?
                fmod(app->scrubberSettings.playTime, app->scrubberSettings.timeMax) + app->scrubberSettings.timeMin :
                app->scrubberSettings.timeMax;
        }
    }

    // Sample Animation Data

    for (int i = 0; i < app->characterData.count; i++)
    {
        if (app->scrubberSettings.sampleMode == 0)
        {
            TransformDataSampleFrameNearest(
                &app->characterData.xformData[i],
                &app->characterData.bvhData[i],
                app->scrubberSettings.playTime,
                app->characterData.scales[i]);
        }
        else if (app->scrubberSettings.sampleMode == 1)
        {
            TransformDataSampleFrameLinear(
                &app->characterData.xformData[i],
                &app->characterData.xformTmp0[i],
                &app->characterData.xformTmp1[i],
                &app->characterData.bvhData[i],
                app->scrubberSettings.playTime,
                app->characterData.scales[i]);
        }
        else
        {
            TransformDataSampleFrameCubic(
                &app->characterData.xformData[i],
                &app->characterData.xformTmp0[i],
                &app->characterData.xformTmp1[i],
                &app->characterData.xformTmp2[i],
                &app->characterData.xformTmp3[i],
                &app->characterData.bvhData[i],
                app->scrubberSettings.playTime,
                app->characterData.scales[i]);
        }

        if (app->scrubberSettings.inplace)
        {
            // Remove Translation on ground Plane
          
            app->characterData.xformData[i].localPositions[0].x = 0.0f;
            app->characterData.xformData[i].localPositions[0].z = 0.0f;
            
            // Attempt to extract rotation around vertical axis (this does not work 
            // for all animations but is pretty effective for almost all of them)
            
            Quaternion verticalRotation = QuaternionInvert(QuaternionNormalize((Quaternion){
                0.0f,
                app->characterData.xformData[i].localRotations[0].y,
                0.0f,
                app->characterData.xformData[i].localRotations[0].w,
            }));
            
            // Remove rotation around vertical axis
            
            app->characterData.xformData[i].localRotations[0] = QuaternionMultiply(
                verticalRotation, 
                app->characterData.xformData[i].localRotations[0]);
        }

        TransformDataForwardKinematics(&app->characterData.xformData[i]);
    }

    // Update Camera

    Vector3 cameraTarget = (Vector3){ 0.0f, 1.0f, 0.0f };

    if (app->characterData.count > 0 &&
        app->camera.track &&
        app->camera.trackBone < app->characterData.xformData[app->characterData.active].jointCount)
    {
        cameraTarget = app->characterData.xformData[app->characterData.active].globalPositions[app->camera.trackBone];
    }

    if (!app->fileDialogState.windowActive)
    {
        OrbitCameraUpdate(
            &app->camera,
            cameraTarget,
            (IsKeyDown(KEY_LEFT_CONTROL) && IsMouseButtonDown(0)) ? GetMouseDelta().x : 0.0f,
            (IsKeyDown(KEY_LEFT_CONTROL) && IsMouseButtonDown(0)) ? GetMouseDelta().y : 0.0f,
            (IsKeyDown(KEY_LEFT_CONTROL) && IsMouseButtonDown(1)) ? GetMouseDelta().x : 0.0f,
            (IsKeyDown(KEY_LEFT_CONTROL) && IsMouseButtonDown(1)) ? GetMouseDelta().y : 0.0f,
            GetMouseWheelMove(),
            GetFrameTime());
    }

    // Create Capsules

    CapsuleDataReset(&app->capsuleData);
    for (int i = 0; i < app->characterData.count; i++)
    {
        CapsuleDataAppendFromTransformData(
            &app->capsuleData,
            &app->characterData.xformData[i],
            app->characterData.radii[i],
            app->characterData.colors[i],
            app->characterData.opacities[i],
            !app->renderSettings.drawEndSites);
    }

    PROFILE_END(Update);

    // Rendering

    Frustum frustum = FrustumFromCameraMatrices(
        GetCameraProjectionMatrix(&app->camera.cam3d, app->screenHeight / app->screenWidth),
        GetCameraViewMatrix(&app->camera.cam3d));

    BeginDrawing();

    PROFILE_BEGIN(Rendering);

    ClearBackground(app->renderSettings.backgroundColor);

    BeginMode3D(app->camera.cam3d);

    // Set shader uniforms that don't change based on the object being drawn

    Vector3 sunColorValue = { app->renderSettings.sunColor.r / 255.0f, app->renderSettings.sunColor.g / 255.0f, app->renderSettings.sunColor.b / 255.0f };
    Vector3 skyColorValue = { app->renderSettings.skyColor.r / 255.0f, app->renderSettings.skyColor.g / 255.0f, app->renderSettings.skyColor.b / 255.0f };
    float objectSpecularity = 0.5f;
    float objectGlossiness = 10.0f;
    float objectOpacity = 1.0f;

    Vector3 sunLightPosition = Vector3RotateByQuaternion((Vector3){ 0.0f, 0.0f, 1.0f }, QuaternionFromAxisAngle((Vector3){ 0.0f, 1.0f, 0.0f }, app->renderSettings.sunAzimuth));
    Vector3 sunLightAxis = Vector3Normalize(Vector3CrossProduct(sunLightPosition, (Vector3){ 0.0f, 1.0f, 0.0f }));
    Vector3 sunLightDir = Vector3Negate(Vector3RotateByQuaternion(sunLightPosition, QuaternionFromAxisAngle(sunLightAxis, app->renderSettings.sunAltitude)));

    SetShaderValue(app->shader, app->uniforms.cameraPosition, &app->camera.cam3d.position, SHADER_UNIFORM_VEC3);
    SetShaderValue(app->shader, app->uniforms.exposure, &app->renderSettings.exposure, SHADER_UNIFORM_FLOAT);
    SetShaderValue(app->shader, app->uniforms.sunDir, &sunLightDir, SHADER_UNIFORM_VEC3);
    SetShaderValue(app->shader, app->uniforms.sunStrength, &app->renderSettings.sunLightStrength, SHADER_UNIFORM_FLOAT);
    SetShaderValue(app->shader, app->uniforms.sunColor, &sunColorValue, SHADER_UNIFORM_VEC3);
    SetShaderValue(app->shader, app->uniforms.skyStrength, &app->renderSettings.skyLightStrength, SHADER_UNIFORM_FLOAT);
    SetShaderValue(app->shader, app->uniforms.skyColor, &skyColorValue, SHADER_UNIFORM_VEC3);
    SetShaderValue(app->shader, app->uniforms.ambientStrength, &app->renderSettings.ambientLightStrength, SHADER_UNIFORM_FLOAT);
    SetShaderValue(app->shader, app->uniforms.groundStrength, &app->renderSettings.groundLightStrength, SHADER_UNIFORM_FLOAT);
    SetShaderValue(app->shader, app->uniforms.objectSpecularity, &objectSpecularity, SHADER_UNIFORM_FLOAT);
    SetShaderValue(app->shader, app->uniforms.objectGlossiness, &objectGlossiness, SHADER_UNIFORM_FLOAT);
    SetShaderValue(app->shader, app->uniforms.objectOpacity, &objectOpacity, SHADER_UNIFORM_FLOAT);
    SetShaderValue(app->shader, app->uniforms.aoLookupResolution, &app->capsuleData.aoLookupResolution, SHADER_UNIFORM_VEC2);
    SetShaderValue(app->shader, app->uniforms.shadowLookupResolution, &app->capsuleData.shadowLookupResolution, SHADER_UNIFORM_VEC2);
    SetShaderValueTexture(app->shader, app->uniforms.aoLookupTable, app->capsuleData.aoLookupTable);
    SetShaderValueTexture(app->shader, app->uniforms.shadowLookupTable, app->capsuleData.shadowLookupTable);
    
    // Draw Ground

    PROFILE_BEGIN(RenderingGround);

    if (app->renderSettings.drawChecker)
    {
        int groundIsCapsule = 0;
        int groundIsCharacter = 0;
        Vector3 groundColor = { 1.0f, 1.0f, 1.0f };

        SetShaderValue(app->shader, app->uniforms.isCapsule, &groundIsCapsule, SHADER_UNIFORM_INT);
        SetShaderValue(app->shader, app->uniforms.objectColor, &groundColor, SHADER_UNIFORM_VEC3);
        SetShaderValue(app->shader, app->uniforms.isCharacter, &groundIsCharacter, SHADER_UNIFORM_INT);

        // Draw ground in a grid of 10x10, 2 meter wide segments.
        
        for (int i = 0; i < 11; i++)
        {
            for (int j = 0; j < 11; j++)
            {
                // Check if we can cull ground segment

                Vector3 groundSegmentPosition =
                {
                    (((float)i / 10) - 0.5f) * 20.0f,
                    0.0f,
                    (((float)j / 10) - 0.5f) * 20.0f,
                };                
                
                if (!FrustumContainsSphere(frustum, groundSegmentPosition, sqrtf(2.0f)))
                {
                    continue;
                }

                PROFILE_BEGIN(RenderingGroundSegment);
                
                // Gather all capsules casting AO on this ground segment
                
                PROFILE_BEGIN(RenderingGroundSegmentAO);
                
                app->capsuleData.aoCapsuleCount = 0;
                if (app->renderSettings.drawCapsules && app->renderSettings.drawAO)
                {  
                    CapsuleDataUpdateAOCapsulesForGroundSegment(&app->capsuleData, groundSegmentPosition);
                }
                int aoCapsuleCount = MinInt(app->capsuleData.aoCapsuleCount, AO_CAPSULES_MAX);

                PROFILE_END(RenderingGroundSegmentAO);
                
                SetShaderValue(app->shader, app->uniforms.aoCapsuleCount, &aoCapsuleCount, SHADER_UNIFORM_INT);
                SetShaderValueV(app->shader, app->uniforms.aoCapsuleStarts, app->capsuleData.aoCapsuleStarts, SHADER_UNIFORM_VEC3, aoCapsuleCount);
                SetShaderValueV(app->shader, app->uniforms.aoCapsuleVectors, app->capsuleData.aoCapsuleVectors, SHADER_UNIFORM_VEC3, aoCapsuleCount);
                SetShaderValueV(app->shader, app->uniforms.aoCapsuleRadii, app->capsuleData.aoCapsuleRadii, SHADER_UNIFORM_FLOAT, aoCapsuleCount);
                
                // Gather all capsules casting shadows on this ground segment
                
                PROFILE_BEGIN(RenderingGroundSegmentShadow);

                app->capsuleData.shadowCapsuleCount = 0;
                if (app->renderSettings.drawCapsules && app->renderSettings.drawShadows)
                {
                    CapsuleDataUpdateShadowCapsulesForGroundSegment(&app->capsuleData, groundSegmentPosition, sunLightDir, app->renderSettings.sunLightConeAngle);
                }
                int shadowCapsuleCount = MinInt(app->capsuleData.shadowCapsuleCount, SHADOW_CAPSULES_MAX);

                PROFILE_END(RenderingGroundSegmentShadow);
                
                SetShaderValue(app->shader, app->uniforms.shadowCapsuleCount, &shadowCapsuleCount, SHADER_UNIFORM_INT);
                SetShaderValueV(app->shader, app->uniforms.shadowCapsuleStarts, app->capsuleData.shadowCapsuleStarts, SHADER_UNIFORM_VEC3, shadowCapsuleCount);
                SetShaderValueV(app->shader, app->uniforms.shadowCapsuleVectors, app->capsuleData.shadowCapsuleVectors, SHADER_UNIFORM_VEC3, shadowCapsuleCount);
                SetShaderValueV(app->shader, app->uniforms.shadowCapsuleRadii, app->capsuleData.shadowCapsuleRadii, SHADER_UNIFORM_FLOAT, shadowCapsuleCount);
                
                // Draw

                DrawModel(app->groundPlaneModel, groundSegmentPosition, 1.0f, WHITE);

                PROFILE_END(RenderingGroundSegment);
            }
        }
    }

    PROFILE_END(RenderingGround);

    // Draw Capsules

    PROFILE_BEGIN(RenderingCapsules);

    if (app->renderSettings.drawCapsules)
    {
        // Depth sort back to front for transparency

        for (int i = 0; i < app->capsuleData.capsuleCount; i++)
        {
            app->capsuleData.capsuleSort[i].index = i;
            app->capsuleData.capsuleSort[i].value = Vector3Distance(app->camera.cam3d.position, app->capsuleData.capsulePositions[i]);
        }

        qsort(app->capsuleData.capsuleSort, app->capsuleData.capsuleCount, sizeof(CapsuleSort), CapsuleSortCompareLess);

        // Render

        int capsuleIsCapsule = 1;
        int capsuleIsCharacter = 0;
        SetShaderValue(app->shader, app->uniforms.isCapsule, &capsuleIsCapsule, SHADER_UNIFORM_INT);
        SetShaderValue(app->shader, app->uniforms.isCharacter, &capsuleIsCharacter, SHADER_UNIFORM_INT);

        for (int i = 0; i < app->capsuleData.capsuleCount; i++)
        {
            int j = app->capsuleData.capsuleSort[i].index;
            
            // Check if we can cull capsule
            
            Vector3 capsulePosition = app->capsuleData.capsulePositions[j];
            float capsuleHalfLength = app->capsuleData.capsuleHalfLengths[j];
            float capsuleRadius = app->capsuleData.capsuleRadii[j];

            if (!FrustumContainsSphere(frustum, capsulePosition, capsuleHalfLength + capsuleRadius))
            {
                continue;
            }
            
            PROFILE_BEGIN(RenderingCapsulesCapsule);
            
            // If capsule is semi-transparent disable depth mask
            
            if (app->capsuleData.capsuleOpacities[j] < 1.0f)
            {
                rlDrawRenderBatchActive();
                rlDisableDepthMask();
            }
            
            // Set shader properties
            
            Quaternion capsuleRotation = app->capsuleData.capsuleRotations[j];
            Vector3 capsuleStart = CapsuleStart(capsulePosition, capsuleRotation, capsuleHalfLength);
            Vector3 capsuleVector = CapsuleVector(capsulePosition, capsuleRotation, capsuleHalfLength);

            SetShaderValue(app->shader, app->uniforms.objectColor, &app->capsuleData.capsuleColors[j], SHADER_UNIFORM_VEC3);
            SetShaderValue(app->shader, app->uniforms.objectOpacity, &app->capsuleData.capsuleOpacities[j], SHADER_UNIFORM_FLOAT);
            SetShaderValue(app->shader, app->uniforms.capsulePosition, &app->capsuleData.capsulePositions[j], SHADER_UNIFORM_VEC3);
            SetShaderValue(app->shader, app->uniforms.capsuleRotation, &app->capsuleData.capsuleRotations[j], SHADER_UNIFORM_VEC4);
            SetShaderValue(app->shader, app->uniforms.capsuleHalfLength, &app->capsuleData.capsuleHalfLengths[j], SHADER_UNIFORM_FLOAT);
            SetShaderValue(app->shader, app->uniforms.capsuleRadius, &app->capsuleData.capsuleRadii[j], SHADER_UNIFORM_FLOAT);
            SetShaderValue(app->shader, app->uniforms.capsuleStart, &capsuleStart, SHADER_UNIFORM_VEC3);
            SetShaderValue(app->shader, app->uniforms.capsuleVector, &capsuleVector, SHADER_UNIFORM_VEC3);
            
            // Find all capsules casting AO on this capsule

            PROFILE_BEGIN(RenderingCapsulesCapsuleAO);

            app->capsuleData.aoCapsuleCount = 0;
            if (app->renderSettings.drawAO)
            {
                CapsuleDataUpdateAOCapsulesForCapsule(&app->capsuleData, j);
            }
            int aoCapsuleCount = MinInt(app->capsuleData.aoCapsuleCount, AO_CAPSULES_MAX);
            
            PROFILE_END(RenderingCapsulesCapsuleAO);

            SetShaderValue(app->shader, app->uniforms.aoCapsuleCount, &aoCapsuleCount, SHADER_UNIFORM_INT);
            SetShaderValueV(app->shader, app->uniforms.aoCapsuleStarts, app->capsuleData.aoCapsuleStarts, SHADER_UNIFORM_VEC3, aoCapsuleCount);
            SetShaderValueV(app->shader, app->uniforms.aoCapsuleVectors, app->capsuleData.aoCapsuleVectors, SHADER_UNIFORM_VEC3, aoCapsuleCount);
            SetShaderValueV(app->shader, app->uniforms.aoCapsuleRadii, app->capsuleData.aoCapsuleRadii, SHADER_UNIFORM_FLOAT, aoCapsuleCount);

            // Find all capsules casting shadows on this capsule

            PROFILE_BEGIN(RenderingCapsulesCapsuleShadow);

            app->capsuleData.shadowCapsuleCount = 0;
            if (app->renderSettings.drawShadows)
            {
                CapsuleDataUpdateShadowCapsulesForCapsule(&app->capsuleData, j, sunLightDir, app->renderSettings.sunLightConeAngle);
            }
            int shadowCapsuleCount = MinInt(app->capsuleData.shadowCapsuleCount, SHADOW_CAPSULES_MAX);

            PROFILE_END(RenderingCapsulesCapsuleShadow);

            SetShaderValue(app->shader, app->uniforms.shadowCapsuleCount, &shadowCapsuleCount, SHADER_UNIFORM_INT);
            SetShaderValueV(app->shader, app->uniforms.shadowCapsuleStarts, app->capsuleData.shadowCapsuleStarts, SHADER_UNIFORM_VEC3, shadowCapsuleCount);
            SetShaderValueV(app->shader, app->uniforms.shadowCapsuleVectors, app->capsuleData.shadowCapsuleVectors, SHADER_UNIFORM_VEC3, shadowCapsuleCount);
            SetShaderValueV(app->shader, app->uniforms.shadowCapsuleRadii, app->capsuleData.shadowCapsuleRadii, SHADER_UNIFORM_FLOAT, shadowCapsuleCount);

            // Draw

            DrawModel(app->capsuleModel, Vector3Zero(), 1.0f, WHITE);
            
            // Reset depth mask if rendered semi-transparent
            
            if (app->capsuleData.capsuleOpacities[j] < 1.0f)
            {
                rlDrawRenderBatchActive();
                rlEnableDepthMask();
            }

            PROFILE_END(RenderingCapsulesCapsule);
        }
    }

    // Draw character mesh
    if (app->characterModel.isLoaded)
    {
        int characterIsCapsule = 0;
        int characterIsCharacter = 1;
        SetShaderValue(app->shader, app->uniforms.isCapsule, &characterIsCapsule, SHADER_UNIFORM_INT);
        SetShaderValue(app->shader, app->uniforms.isCharacter, &characterIsCharacter, SHADER_UNIFORM_INT);

        Model* currModel = &app->characterModel;

        // Draw one model instance for each animation
        for (int i = 0; i < app->characterData.count; i++)
        {
            Vector3 position = {1.0f * i, 0.0f, 0.0f};
            int frame = app->scrubberSettings.playTime / app->characterData.bvhData[i].frameTime;
            UpdateModelAnimation(app->characterModel.model, app->characterData.animData[i], frame);
            DrawModel(app->characterModel.model, position, 1.0f, WHITE);
            break;
        }
    }

    PROFILE_END(RenderingCapsules);

    // Grid

    if (app->renderSettings.drawGrid)
    {
        DrawGrid(20, 1.0f);
    }

    // Origin

    if (app->renderSettings.drawOrigin)
    {
        DrawTransform(
            (Vector3){ 0.0f, 0.01f, 0.0f },
            QuaternionIdentity(),
            1.0f);
    }

    // Disable Depth Test

    rlDrawRenderBatchActive();
    rlDisableDepthTest();

    // Draw Capsule Wireframes

    if (app->renderSettings.drawWireframes)
    {
        DrawWireFrames(&app->capsuleData, DARKGRAY);
    }

    // Draw Bones

    if (app->renderSettings.drawSkeleton)
    {
        for (int i = 0; i < app->characterData.count; i++)
        {
            DrawSkeleton(
                &app->characterData.xformData[i],
                app->renderSettings.drawEndSites,
                DARKGRAY,
                GRAY);
        }
    }

    // Draw Joint Transforms

    if (app->renderSettings.drawTransforms)
    {
        for (int i = 0; i < app->characterData.count; i++)
        {
            DrawTransforms(&app->characterData.xformData[i]);
        }
    }

    // Re-Enable Depth Test

    rlDrawRenderBatchActive();
    rlEnableDepthTest();

    // Rendering Done

    EndMode3D();

    PROFILE_END(Rendering);

    // Draw UI

    PROFILE_BEGIN(Gui);

    if (app->renderSettings.drawUI)
    {
        if (app->fileDialogState.windowActive) { GuiLock(); }

        // Error Message

        DrawText(app->errMsg, 250, 20, 15, RED);

        if (app->characterData.count == 0)
        {
            DrawText("Drag and Drop .bvh files to open them.",
              app->screenWidth / 2 - 300, app->screenHeight / 2 - 15, 30, DARKGRAY);
        }

        // Render Settings

        GuiRenderSettings(&app->renderSettings, &app->capsuleData, app->screenWidth, app->screenHeight);

        // FPS

        if (app->renderSettings.drawFPS)
        {
            DrawFPS(230, 10);
        }

        // Camera Settings

        GuiOrbitCamera(&app->camera, &app->characterData, app->argc, app->argv);

        // Characters

        GuiCharacterData(&app->characterData, &app->fileDialogState, &app->scrubberSettings, app->errMsg, app->argc, app->argv);

        // Color Picker

        if (app->characterData.colorPickerActive)
        {
            GuiGroupBox((Rectangle){ app->screenWidth - 180, 450, 160, 140 }, "Color Picker");
            GuiColorPicker((Rectangle){ app->screenWidth - 165, 465, 110, 110 }, NULL, &app->characterData.colors[app->characterData.active]);
        }

        // Scrubber

        GuiScrubberSettings(&app->scrubberSettings, &app->characterData, app->screenWidth, app->screenHeight);

        // File Dialog

        if (app->fileDialogState.windowActive) { GuiUnlock(); }
        
        GuiWindowFileDialog(&app->fileDialogState);
    }

    PROFILE_END(Gui);

#if defined(ENABLE_PROFILE) && defined(_WIN32)

    // Display Profile Records

    PROFILE_TICKERS_UPDATE();
    
    for (int i = 0; i < globalProfileRecords.num; i++)
    {
        GuiLabel((Rectangle){ 260, 10 + (float)i * 20, 200, 20 }, globalProfileRecords.records[i]->name);
        GuiLabel((Rectangle){ 450, 10 + (float)i * 20, 100, 20 }, TextFormat("%6.1f us", globalProfileTickers.times[i]));
        GuiLabel((Rectangle){ 550, 10 + (float)i * 20, 100, 20 }, TextFormat("%i calls", globalProfileTickers.samples[i]));
    }
#endif

    // Done

    EndDrawing();
}

//----------------------------------------------------------------------------------
// Main
//----------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    PROFILE_INIT();
    PROFILE_TICKERS_INIT();
    
    // Init Application State
    
    ApplicationState app;
    app.argc = argc;
    app.argv = argv;    
    app.screenWidth = ArgInt(argc, argv, "screenWidth", 1280);
    app.screenHeight = ArgInt(argc, argv, "screenHeight", 720);
    
    // Init Window

    SetConfigFlags(FLAG_VSYNC_HINT);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(app.screenWidth, app.screenHeight, "BVHView");
    SetTargetFPS(60);

    // Camera

    OrbitCameraInit(&app.camera, argc, argv);

    // Shader

    app.shader = LoadShaderFromMemory(shaderVS, shaderFS);
    ShaderUniformsInit(&app.uniforms, app.shader);

    // Models

    app.groundPlaneMesh = GenMeshPlane(2.0f, 2.0f, 1, 1);
    app.groundPlaneModel = LoadModelFromMesh(app.groundPlaneMesh);
    app.groundPlaneModel.materials[0].shader = app.shader;

    app.capsuleModel = LoadOBJFromMemory(capsuleOBJ);
    app.capsuleModel.materials[0].shader = app.shader;

    CharacterDataInit(&app.characterData, argc, argv);
    CharacterModelInit(&app.characterModel);
    app.characterModel.model.materials[0].shader = app.shader;
    app.characterModel.model.materials[1].shader = app.shader;

    // Capsule Data

    CapsuleDataInit(&app.capsuleData);

    // Scrubber Settings

    ScrubberSettingsInit(&app.scrubberSettings, argc, argv);

    // Render Settings

    RenderSettingsInit(&app.renderSettings, argc, argv);
    CapsuleDataUpdateShadowLookupTable(&app.capsuleData, app.renderSettings.sunLightConeAngle);

    // File Dialog

    app.fileDialogState = InitGuiWindowFileDialog(GetWorkingDirectory());

    // Reset Error Message
    
    app.errMsg[0] = '\0';

    // Load any files given as command line arguments

    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-') { continue; }

        CharacterDataLoadFromFile(&app.characterData, argv[i], app.errMsg, 512);
    }

    // If any characters loaded, update capsules and scrubber

    if (app.characterData.count > 0)
    {
        app.characterData.active = app.characterData.count - 1;

        CapsuleDataUpdateForCharacters(&app.capsuleData, &app.characterData);
        ScrubberSettingsRecomputeLimits(&app.scrubberSettings, &app.characterData);
        ScrubberSettingsInitMaxs(&app.scrubberSettings, &app.characterData);
        
        char windowTitle[512];
        snprintf(windowTitle, 512, "%s - BVHView", app.characterData.filePaths[app.characterData.active]);
        SetWindowTitle(windowTitle);
    }

    // Game Loop

#if defined(PLATFORM_WEB)
    emscripten_set_main_loop_arg(ApplicationUpdate, &app, 0, 1);
#else
    while (!WindowShouldClose())
    {
        ApplicationUpdate(&app);
    }
#endif

    // Unload and finish

    CapsuleDataFree(&app.capsuleData);
    CharacterDataFree(&app.characterData);
    CharacterModelFree(&app.characterModel);

    UnloadModel(app.capsuleModel);
    UnloadModel(app.groundPlaneModel);
    UnloadShader(app.shader);

    CloseWindow();

    return 0;
}
