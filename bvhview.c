#include <assert.h>
#include <stdlib.h>

#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#define GUI_WINDOW_FILE_DIALOG_IMPLEMENTATION
#include "../examples/custom_file_dialog/gui_window_file_dialog.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

//--------------------------------------
// Additional C Functions
//--------------------------------------

#define PIf 3.14159265358979323846f

static inline int max(int x, int y)
{
    return x > y ? x : y;
}

static inline int min(int x, int y)
{
    return x < y ? x : y;
}

static inline int clamp(int x, int min, int max)
{
    return x < min ? min : x > max ? max : x;
}

static inline float maxf(float x, float y)
{
    return x > y ? x : y;
}

static inline float minf(float x, float y)
{
    return x < y ? x : y;
}

static inline float squaref(float x)
{
    return x * x;
}

static inline float clampf(float x, float min, float max)
{
    return x > max ? max : x < min ? min : x;
}

static inline float deg2radf(const float x)
{
    return x * (PIf / 180.0f);
}

static inline float fast_acosf(float x)
{
    float y = fabs(x);
    float p = -0.1565827f * y + 1.570796f;
    p *= sqrtf(maxf(1.0f - y, 0.0f));
    return x >= 0.0f ? p : PIf - p;
}

static inline float fast_positive_acosf(float x)
{
    float p = -0.1565827f * x + 1.570796f;
    return p * sqrtf(maxf(1.0f - x, 0.0f));
}  

static inline float fast_positive_atanf(float x)
{
    float w = x > 1.0f ? 1.0f / x : x;
    float y = (PIf / 4.0f)*w - w*(w - 1.0f)*(0.2447f + 0.0663f*w);
    return fabs(x > 1.0f ? PIf / 2.0f - y : y);
}

//--------------------------------------
// Additional Raylib Functions
//--------------------------------------

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
        float halfangle = acosf(clampf(q.w, -1.0f, 1.0f));
        return Vector3Scale((Vector3){ q.x, q.y, q.z }, halfangle / length);
    }
}

static inline Vector3 QuaternionToScaleAngleAxis(Quaternion q)
{
    return Vector3Scale(QuaternionLog(q), 2.0f);
}

static inline Quaternion QuaternionFromScaledAngleAxis(Vector3 v)
{
    return QuaternionExp(Vector3Scale(v, 0.5f));
}

static inline Quaternion QuaternionHermite(Quaternion r0, Quaternion r1, Vector3 v0, Vector3 v1, float alpha)
{
    float x = alpha;
    float w1 = 3*x*x - 2*x*x*x;
    float w2 = x*x*x - 2*x*x + x;
    float w3 = x*x*x - x*x;
    
    Vector3 r1_sub_r0 = QuaternionToScaleAngleAxis(QuaternionAbsolute(QuaternionMultiply(r1, QuaternionInvert(r0))));   
    
    return QuaternionMultiply(QuaternionFromScaledAngleAxis(
        Vector3Add(Vector3Add(Vector3Scale(r1_sub_r0, w1), Vector3Scale(v0, w2)), Vector3Scale(v1, w3))), r0);
}

static inline Quaternion QuaternionInterpolateCubic(Quaternion r0, Quaternion r1, Quaternion r2, Quaternion r3, float alpha)
{
    Vector3 r1_sub_r0 = QuaternionToScaleAngleAxis(QuaternionAbsolute(QuaternionMultiply(r1, QuaternionInvert(r0))));
    Vector3 r2_sub_r1 = QuaternionToScaleAngleAxis(QuaternionAbsolute(QuaternionMultiply(r2, QuaternionInvert(r1))));
    Vector3 r3_sub_r2 = QuaternionToScaleAngleAxis(QuaternionAbsolute(QuaternionMultiply(r3, QuaternionInvert(r2))));
  
    Vector3 v1 = Vector3Scale(Vector3Add(r1_sub_r0, r2_sub_r1), 0.5f);
    Vector3 v2 = Vector3Scale(Vector3Add(r2_sub_r1, r3_sub_r2), 0.5f);
    return QuaternionHermite(r1, r2, v1, v2, alpha);
}

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
        QuaternionFromAxisAngle((Vector3){ 1.0f, 0.0f, 0.0f }, PIf) :
        QuaternionNormalize(o);
}

//--------------------------------------
// Camera
//--------------------------------------

static inline float OrbitCameraUpdateAzimuth(float azimuth, float mouseDx, float dt)
{
    return azimuth + 1.0f * dt * -mouseDx;
}

static inline float OrbitCameraUpdateAltitude(float altitude, float mouseDy, float dt)
{
    return clampf(altitude + 1.0f * dt * mouseDy, 0.0, 0.4f * PIf);
}

static inline float OrbitCameraUpdateDistance(float distance, float dt)
{
    return clampf(distance +  20.0f * dt * -GetMouseWheelMove(), 0.1f, 100.0f);
}

static inline void OrbitCameraUpdate(
    Camera3D* cam, 
    float* cameraAzimuth,
    float* cameraAltitude,
    float* cameraDistance,
    Vector3 target,
    float mouse_dx,
    float mouse_dy,
    float dt)
{
    *cameraAzimuth = OrbitCameraUpdateAzimuth(*cameraAzimuth, mouse_dx, dt);
    *cameraAltitude = OrbitCameraUpdateAltitude(*cameraAltitude, mouse_dy, dt);
    *cameraDistance = OrbitCameraUpdateDistance(*cameraDistance, dt);
    
    Quaternion rotationAzimuth = QuaternionFromAxisAngle((Vector3){0, 1, 0}, *cameraAzimuth);
    Vector3 position = Vector3RotateByQuaternion((Vector3){0, 0, *cameraDistance}, rotationAzimuth);
    Vector3 axis = Vector3Normalize(Vector3CrossProduct(position, (Vector3){0, 1, 0}));
    
    Quaternion rotationAltitude = QuaternionFromAxisAngle(axis, *cameraAltitude);
    
    Vector3 eye = Vector3Add(target, Vector3RotateByQuaternion(position, rotationAltitude));

    cam->target = target;
    cam->position = eye;
}

//--------------------------------------
// BVH File Data
//--------------------------------------

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

typedef struct 
{
    int parent;
    char* name;
    Vector3 offset;
    int channelCount;
    char channels[CHANNELS_MAX];
    bool endSite;
    
} BVHJointData;

static void BVHJointDataInit(BVHJointData* data)
{
    data->parent = -1;
    data->name = NULL;
    data->offset = (Vector3){ 0.0f, 0.0f, 0.0f };
    data->channelCount = 0;
    data->endSite = false;
}

static void BVHJointDataRename(BVHJointData* data, const char* name)
{
    data->name = realloc(data->name, strlen(name) + 1);
    strcpy(data->name, name);
}

static void BVHJointDataFree(BVHJointData* data)
{
    free(data->name);
}

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
    
    // Extra Data
   
    char* jointNamesCombo;

} BVHData;

static void BVHDataInit(BVHData* bvh)
{
    bvh->jointCount = 0;
    bvh->joints = NULL;
    bvh->frameCount = 0;
    bvh->channelCount = 0;
    bvh->frameTime = 0.0f;
    bvh->motionData = NULL;
    bvh->jointNamesCombo = NULL;
}

static void BVHDataFree(BVHData* bvh)
{
    for (int i = 0; i < bvh->jointCount; i++)
    {
        BVHJointDataFree(&bvh->joints[i]);
    }
    free(bvh->joints);
    
    free(bvh->motionData);
    free(bvh->jointNamesCombo);
}

static int BVHDataAddJoint(BVHData* bvh)
{
    bvh->joints = (BVHJointData*)realloc(bvh->joints, (bvh->jointCount + 1) * sizeof(BVHJointData));
    bvh->jointCount++;
    BVHJointDataInit(&bvh->joints[bvh->jointCount - 1]);
    return bvh->jointCount - 1;
}

static void BVHDataComputeJointNamesCombo(BVHData* bvh)
{
    int total_size = 0;
    for (int i = 0; i < bvh->jointCount; i++)
    {
        total_size += (i > 0 ? 1 : 0) + strlen(bvh->joints[i].name);   
    }
    total_size++;
  
    bvh->jointNamesCombo = malloc(total_size);
    bvh->jointNamesCombo[0] = '\0';
    for (int i = 0; i < bvh->jointCount; i++)
    {
        if (i > 0)
        {
            strcat(bvh->jointNamesCombo, ";");                  
        }
        strcat(bvh->jointNamesCombo, bvh->joints[i].name);
    }
}

//--------------------------------------
// Parser
//--------------------------------------

enum
{
    PARSER_ERR_MAX = 512,
};

typedef struct {
  
    const char* filename;
    int offset;
    const char* data;
    int row;
    int col;
    char err[PARSER_ERR_MAX];
  
} Parser;

static void ParserInit(Parser* par, const char* filename, const char* data)
{
    par->filename = filename;
    par->offset = 0;
    par->data = data;
    par->row = 0;
    par->col = 0;
    par->err[0] = '\0';
}

static char ParserPeek(const Parser* par)
{
    return *(par->data + par->offset);
}

static char ParserPeekForward(const Parser* par, int steps)
{
    return *(par->data + par->offset + steps);
}

static bool ParserMatch(const Parser* par, char match)
{
    return match == *(par->data + par->offset);
}

static bool ParserOneOf(const Parser* par, const char* matches)
{
    return strchr(matches, *(par->data + par->offset));
}

static bool ParserStartsWith(const Parser* par, const char* prefix)
{
    const char* start = par->data + par->offset;
    while (*prefix)
    {
        if (*prefix != *start) { return false; }
        prefix++;
        start++;
    }

    return true;
}

static void ParserInc(Parser* par)
{
    if (*(par->data + par->offset) == '\n')
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

static void ParserAdvance(Parser* par, int num)
{
    for (int i = 0; i < num; i++) { ParserInc(par); }
}

#define ParserError(par, fmt, ...) \
    snprintf(par->err, PARSER_ERR_MAX, "%s:%i:%i: error: " fmt, par->filename, par->row, par->col, ##__VA_ARGS__)

//--------------------------------------
// BVH Parser
//--------------------------------------

// TODO: name of char for error

static void BVHParseWhitespace(Parser* par)
{
    while (ParserOneOf(par, " \r\t\v")) { ParserInc(par); }
}

static bool BVHParseString(Parser* par, const char* string)
{
    if (ParserStartsWith(par, string))
    {
        ParserAdvance(par, strlen(string));
        return true;
    }
    else
    {
        ParserError(par, "expected '%s' at '%c'", string, ParserPeek(par));
        return false;
    }
}

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
        ParserError(par, "expected newline at '%c'", ParserPeek(par));
        return false;
    }    
}

static bool BVHParseJointName(BVHJointData* jnt, Parser* par)
{
    BVHParseWhitespace(par);

    char buffer[256];
    int chrnum = 0;
    while (chrnum < 255 && ParserOneOf(par, 
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"))
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
        ParserError(par, "expected joint name at '%c'", ParserPeek(par));
        return false;
    }   
}

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
        ParserError(par, "expected float at '%c'", ParserPeek(par));
        return false;
    }
}

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
        ParserError(par, "expected integer at '%c'", ParserPeek(par));
        return false;
    }
}

static bool BVHParseJointOffset(BVHJointData* jnt, Parser* par)
{
    if (!BVHParseString(par, "OFFSET")) { return false; }
    if (!BVHParseFloat(&jnt->offset.x, par)) { return false; }
    if (!BVHParseFloat(&jnt->offset.y, par)) { return false; }
    if (!BVHParseFloat(&jnt->offset.z, par)) { return false; }
    if (!BVHParseNewline(par)) { return false; }
    return true;
}

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

static bool BVHParseChannel(char* channel, Parser* par)
{
    BVHParseWhitespace(par);

    if (ParserPeek(par) == '\0')
    {
        ParserError(par, "expected channel at end of file");
        return false;
    }
  
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

static bool BVHParseJoints(BVHData* bvh, int parent, Parser* par)
{
    while (ParserOneOf(par, "JE")) // Either "JOINT" or "End Site"
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

static bool BVHParseFrames(BVHData* bvh, Parser* par)
{
    if (!BVHParseString(par, "Frames:")) { return false; }
    if (!BVHParseInt(&bvh->frameCount, par)) { return false; }
    if (!BVHParseNewline(par)) { return false; }
    return true;
}

static bool BVHParseFrameTime(BVHData* bvh, Parser* par)
{
    if (!BVHParseString(par, "Frame Time:")) { return false; }
    if (!BVHParseFloat(&bvh->frameTime, par)) { return false; }
    if (!BVHParseNewline(par)) { return false; }
    if (bvh->frameTime == 0.0f) { bvh->frameTime = 1.0f / 60.0f; }
    return true;
}

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

static bool BVHDataLoad(BVHData* bvh, const char* filename, char* errMsg, int errMsgSize)
{
    // Read file Contents
  
    FILE* f = fopen(filename, "r");
    
    if (f == NULL)
    {
        snprintf(errMsg, errMsgSize, "Error: Could not find file '%s'\n", filename);
        return false;
    }
    
    fseek(f, 0, SEEK_END);
    long int length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buffer = malloc(length + 1);
    fread(buffer, 1, length, f);
    buffer[length] = '\n';
    fclose(f);
    
    // Create Parser
    
    Parser par;
    ParserInit(&par, filename, buffer);
    
    // Parse BVH
    
    BVHDataFree(bvh);
    BVHDataInit(bvh);
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
        
        // Compute some additional info
        BVHDataComputeJointNamesCombo(bvh);
    }
    
    return result;
}

//--------------------------------------
// Transform Data
//--------------------------------------

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

static void TransformDataInit(TransformData* data)
{
    data->jointCount = 0;
    data->parents = NULL;
    data->endSite = NULL;
    data->localPositions = NULL;
    data->localRotations = NULL;
    data->globalPositions = NULL;
    data->globalRotations = NULL;
}

static void TransformDataResize(TransformData* data, BVHData* bvh)
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

static void TransformDataFree(TransformData* data)
{
    free(data->parents);
    free(data->endSite);
    free(data->localPositions);
    free(data->localRotations);
    free(data->globalPositions);
    free(data->globalRotations);
}

static void TransformDataSampleFrame(TransformData* data, BVHData* bvh, int frame, float scale)
{
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
                        (Vector3){1, 0, 0}, deg2radf(bvh->motionData[frame * bvh->channelCount + offset])));
                    offset++;
                    break;
                
                case CHANNEL_Y_ROTATION:
                    rotation = QuaternionMultiply(rotation, QuaternionFromAxisAngle(
                        (Vector3){0, 1, 0}, deg2radf(bvh->motionData[frame * bvh->channelCount + offset])));
                    offset++;
                    break;
                    
                case CHANNEL_Z_ROTATION:
                    rotation = QuaternionMultiply(rotation, QuaternionFromAxisAngle(
                        (Vector3){0, 0, 1}, deg2radf(bvh->motionData[frame * bvh->channelCount + offset])));
                    offset++;
                    break;
            }
        }
        
        data->localPositions[i] = position;
        data->localRotations[i] = rotation;
    }
    
    assert(offset == bvh->channelCount);
}

static void TransformDataSampleFrameNearest(TransformData* data, BVHData* bvh, float time, float scale)
{
    int frame = clamp((int)(time / bvh->frameTime + 0.5f), 0, bvh->frameCount - 1);
    
    TransformDataSampleFrame(data, bvh, frame, scale);
}

static void TransformDataSampleFrameLinear(
    TransformData* data, 
    TransformData* tmp0, 
    TransformData* tmp1, 
    BVHData* bvh, 
    float time, 
    float scale)
{
    const float alpha = fmod(time / bvh->frameTime, 1.0f);
    int frame0 = clamp((int)(time / bvh->frameTime) + 0, 0, bvh->frameCount - 1);
    int frame1 = clamp((int)(time / bvh->frameTime) + 1, 0, bvh->frameCount - 1);
    
    TransformDataSampleFrame(tmp0, bvh, frame0, scale);
    TransformDataSampleFrame(tmp1, bvh, frame1, scale);
  
    for (int i = 0; i < data->jointCount; i++)
    {
        data->localPositions[i] = Vector3Lerp(tmp0->localPositions[i], tmp1->localPositions[i], alpha);
        data->localRotations[i] = QuaternionSlerp(tmp0->localRotations[i], tmp1->localRotations[i], alpha);
    }
}

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
    int frame0 = clamp((int)(time / bvh->frameTime) - 1, 0, bvh->frameCount - 1);
    int frame1 = clamp((int)(time / bvh->frameTime) + 0, 0, bvh->frameCount - 1);
    int frame2 = clamp((int)(time / bvh->frameTime) + 1, 0, bvh->frameCount - 1);
    int frame3 = clamp((int)(time / bvh->frameTime) + 2, 0, bvh->frameCount - 1);
    
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

//--------------------------------------
// Capsule Functions
//--------------------------------------

static Vector3 CapsuleStart(Vector3 capsulePosition, Quaternion capsuleRotation, float capsuleHalfLength)
{
    return Vector3Add(capsulePosition, 
        Vector3RotateByQuaternion((Vector3){+capsuleHalfLength, 0.0f, 0.0f}, capsuleRotation));
}

static Vector3 CapsuleEnd(Vector3 capsulePosition, Quaternion capsuleRotation, float capsuleHalfLength)
{
    return Vector3Add(capsulePosition, 
        Vector3RotateByQuaternion((Vector3){-capsuleHalfLength, 0.0f, 0.0f}, capsuleRotation)); 
}

static Vector3 CapsuleVector(Vector3 capsulePosition, Quaternion capsuleRotation, float capsuleHalfLength)
{
    Vector3 capsuleStart = CapsuleStart(capsulePosition, capsuleRotation, capsuleHalfLength);
  
    return Vector3Subtract(Vector3Add(capsulePosition, 
        Vector3RotateByQuaternion((Vector3){-capsuleHalfLength, 0.0f, 0.0f}, capsuleRotation)), capsuleStart);   
}

static Vector3 NearestPointOnLineSegment(
    Vector3 lineStart, 
    Vector3 lineEnd,
    Vector3 point)
{
    Vector3 ab = Vector3Subtract(lineEnd, lineStart);
    Vector3 ap = Vector3Subtract(point, lineStart);
    float lengthsq = Vector3LengthSqr(ab);
    
    if (lengthsq < 1e-8f)
    {
        return lineStart;
    }
    else
    {
        float t = clampf(Vector3DotProduct(ab, ap) / lengthsq, 0.0f, 1.0f);
        return Vector3Add(lineStart, Vector3Scale(ab, t));
    }
}

static void NearestPointBetweenLineSegments(
    Vector3* nearestLine0, 
    Vector3* nearestLine1, 
    Vector3 line0Start, 
    Vector3 line0End,
    Vector3 line1Start,
    Vector3 line1End)
{
    float d0 = Vector3LengthSqr(Vector3Subtract(line1Start, line0Start)); 
    float d1 = Vector3LengthSqr(Vector3Subtract(line1End, line0Start)); 
    float d2 = Vector3LengthSqr(Vector3Subtract(line1Start, line0End)); 
    float d3 = Vector3LengthSqr(Vector3Subtract(line1End, line0End));

    if (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1)
    {
        *nearestLine0 = line0End;
    }
    else
    {
        *nearestLine0 = line0Start;
    }
    
    *nearestLine1 = NearestPointOnLineSegment(line1Start, line1End, *nearestLine0);
    *nearestLine0 = NearestPointOnLineSegment(line0Start, line0End, *nearestLine1);
}

static void NearestPointBetweenLineSegmentAndGroundPlane(
    Vector3* nearestPointOnLine,
    Vector3* nearestPointOnGroundPlane,
    Vector3 lineStart, 
    Vector3 lineEnd)
{
    Vector3 lineVector = Vector3Subtract(lineEnd, lineStart);
    
    if (fabs(lineVector.y) < 1e-8f)
    {
        *nearestPointOnLine = lineStart;
        *nearestPointOnGroundPlane = (Vector3){ nearestPointOnLine->x, 0.0f, nearestPointOnLine->z };
        return;
    }
    
    float alpha = clampf((-lineStart.y) / lineVector.y, 0.0f, 1.0f);
    
    *nearestPointOnLine = Vector3Add(lineStart, Vector3Scale(lineVector, alpha));
    *nearestPointOnGroundPlane = (Vector3){ nearestPointOnLine->x, 0.0f, nearestPointOnLine->z };
}

static void NearestPointBetweenLineSegmentAndGroundSegment(
    Vector3* nearestPointOnLine,
    Vector3* nearestPointOnGround,
    Vector3 lineStart, 
    Vector3 lineEnd,
    Vector3 groundMins,
    Vector3 groundMaxs)
{
    // Check Against Plane
  
    NearestPointBetweenLineSegmentAndGroundPlane(
        nearestPointOnLine, 
        nearestPointOnGround,
        lineStart,
        lineEnd);
        
    // If point is inside plane bounds it must be the nearest
  
    if (nearestPointOnGround->x >= groundMins.x &&
        nearestPointOnGround->x <= groundMaxs.x &&
        nearestPointOnGround->z <= groundMaxs.z &&
        nearestPointOnGround->z <= groundMaxs.z)
    {
        return;
    }
    
    // Check against four edges

    Vector3 nearestPointOnLine0, nearestPointOnLine1, nearestPointOnLine2, nearestPointOnLine3;
    Vector3 nearestPointOnGround0, nearestPointOnGround1, nearestPointOnGround2, nearestPointOnGround3;
    float distance0, distance1, distance2, distance3;
  
    NearestPointBetweenLineSegments(
        &nearestPointOnLine0,
        &nearestPointOnGround0,
        lineStart, lineEnd,
        (Vector3){ groundMins.x, 0.0f, groundMins.z },
        (Vector3){ groundMins.x, 0.0f, groundMaxs.z });
  
    NearestPointBetweenLineSegments(
        &nearestPointOnLine1,
        &nearestPointOnGround1,
        lineStart, lineEnd,
        (Vector3){ groundMins.x, 0.0f, groundMaxs.z },
        (Vector3){ groundMaxs.x, 0.0f, groundMaxs.z });
  
    NearestPointBetweenLineSegments(
        &nearestPointOnLine2,
        &nearestPointOnGround2,
        lineStart, lineEnd,
        (Vector3){ groundMaxs.x, 0.0f, groundMaxs.z },
        (Vector3){ groundMaxs.x, 0.0f, groundMins.z });
  
    NearestPointBetweenLineSegments(
        &nearestPointOnLine3,
        &nearestPointOnGround3,
        lineStart, lineEnd,
        (Vector3){ groundMaxs.x, 0.0f, groundMins.z },
        (Vector3){ groundMins.x, 0.0f, groundMins.z });
    
    distance0 = Vector3Distance(nearestPointOnLine0, nearestPointOnGround0);
    distance1 = Vector3Distance(nearestPointOnLine1, nearestPointOnGround1);
    distance2 = Vector3Distance(nearestPointOnLine2, nearestPointOnGround2);
    distance3 = Vector3Distance(nearestPointOnLine3, nearestPointOnGround3);
    
    if (distance0 <= distance1 && distance0 <= distance2 && distance0 <= distance3)
    {
        *nearestPointOnLine = nearestPointOnLine0;
        *nearestPointOnGround = nearestPointOnGround0;
        return;
    }
    
    if (distance1 <= distance0 && distance1 <= distance2 && distance1 <= distance3)
    {
        *nearestPointOnLine = nearestPointOnLine1;
        *nearestPointOnGround = nearestPointOnGround1;
        return;
    }
    
    if (distance2 <= distance0 && distance2 <= distance1 && distance2 <= distance3)
    {
        *nearestPointOnLine = nearestPointOnLine2;
        *nearestPointOnGround = nearestPointOnGround2;
        return;
    }
    
    if (distance3 <= distance0 && distance3 <= distance1 && distance3 <= distance2)
    {
        *nearestPointOnLine = nearestPointOnLine3;
        *nearestPointOnGround = nearestPointOnGround3;
        return;
    }
}

static float SphereOcclusion(Vector3 pos, Vector3 nor, Vector3 sph, float rad)
{
    Vector3 di = Vector3Subtract(sph, pos);
    float l = maxf(Vector3Length(di), rad);
    float nl = Vector3DotProduct(nor, Vector3Normalize(di));
    float h  = l / rad;
    float h2 = h*h;
    
    float res;
    float k2 = 1.0 - h2*nl*nl;
    if (k2 > 0.001)                                                        
    {                                                                      
        res = nl * fast_acosf(-nl*sqrtf((h2 - 1.0f) / (1.0f - nl*nl))) - sqrtf(k2*(h2 - 1.0f));   
        res = (res / h2 + fast_positive_atanf(sqrt(k2 / (h2 - 1.0f)))) / PIf; 
    }                                                                      
    else                                                                   
    {                                                                      
        res = maxf(0.0, nl) / h2;                                
    }      

    return 1.0f - clampf(res, 0.0f, 1.0f);
}

static float CapsuleSphereIntersectionArea(
    float cosCap1, float cosCap2,
    float cap2, float cosCist)
{
    float r1 = fast_positive_acosf(cosCap1);
    float r2 = cap2;
    float d  = fast_acosf(cosCist);

    if (minf(r1, r2) <= maxf(r1, r2) - d)
    {
        return 1.0f - maxf(cosCap1, cosCap2);
    }
    else if (r1 + r2 <= d)
    {
        return 0.0f; 
    }

    float delta = fabs(r1 - r2);
    float x = 1.0f - clampf((d - delta) / maxf(r1 + r2 - delta, 0.0001f), 0.0f, 1.0f);
    float area = squaref(x) * (-2.0f * x + 3.0f);

    return area * (1.0f - maxf(cosCap1, cosCap2));
}

static float SphereDirectionalOcclusion(
    Vector3 pos, Vector3 sphere, float radius,
    Vector3 coneDir, float coneAngle)
{
    Vector3 occluder = Vector3Subtract(sphere, pos);
    float occluderLen2 = Vector3DotProduct(occluder, occluder);
    Vector3 occluderDir = Vector3Scale(occluder, 1.0f / maxf(sqrtf(occluderLen2), 1e-8f));

    float cosPhi = Vector3DotProduct(occluderDir, Vector3Negate(coneDir));
    float cosTheta = sqrtf(occluderLen2 / (squaref(radius) + occluderLen2));
    float cosCone = cosf(coneAngle / 2.0f) / (1.0f - 1e-8f);

    return 1.0f - CapsuleSphereIntersectionArea(
        cosTheta, cosCone, coneAngle / 2.0f, cosPhi) / (1.0f - cosCone);
}

static float CapsuleDirectionalOcclusion(
    Vector3 pos, Vector3 capStart, Vector3 capVec,
    float capRadius, Vector3 coneDir, float coneAngle)
{
    Vector3 ba = capVec;
    Vector3 pa = Vector3Subtract(capStart, pos);
    float a = Vector3DotProduct(Vector3Negate(coneDir), ba);
    float t = clampf(
        Vector3DotProduct(pa, Vector3Subtract(Vector3Scale(Vector3Negate(coneDir), a), ba)) / 
        (Vector3DotProduct(ba, ba) - a * a), 
        0.0f, 1.0f);

    return SphereDirectionalOcclusion(
        pos, Vector3Add(capStart, Vector3Scale(ba, t)), capRadius, coneDir, coneAngle);
}

float CapsuleRayIntersect(
    Vector3 capStart, Vector3 capVec, float capRadius,
    Vector3 rayStart, Vector3 rayDir)
{
    Vector3 ba = capVec;
    Vector3 oa = Vector3Subtract(rayStart, capStart);

    float baba = Vector3DotProduct(ba, ba);
    float bard = Vector3DotProduct(ba, rayDir);
    float baoa = Vector3DotProduct(ba, oa);
    float rdoa = Vector3DotProduct(rayDir, oa);
    float oaoa = Vector3DotProduct(oa, oa);

    float a = baba      - bard*bard;
    float b = baba*rdoa - baoa*bard;
    float c = baba*oaoa - baoa*baoa - capRadius*capRadius*baba;
    float h = b*b - a*c;
    
    if (h >= 0.0)
    {
        float t = (-b - sqrtf(h))/a;
        float y = baoa + t*bard;
        
        // body
        if (y > 0.0 && y < baba)
        {
            return t > 1e-4 ? 0.0 : 1.0;
        }
        
        // caps
        Vector3 oc = y <= 0.0 ? oa : Vector3Subtract(rayStart, Vector3Add(capStart, capVec));
        b = Vector3DotProduct(rayDir, oc);
        c = Vector3DotProduct(oc, oc) - capRadius*capRadius;
        h = b*b - c;
        
        if (h > 0.0)
        {
            t = -b - sqrtf(h);
            return t > 1e-4 ? 0.0 : 1.0;
        }
    }
    
    return 1.0;
}

//--------------------------------------
// Capsule Data
//--------------------------------------

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

typedef struct
{
    int capsuleCount;
    Vector3* capsulePositions;
    Quaternion* capsuleRotations;
    float* capsuleRadii;
    float* capsuleHalfLengths;
    Vector3* capsuleColors;
    CapsuleSort* capsuleSort;
    
    int aoCapsuleCount;
    Vector3* aoCapsuleStarts;
    Vector3* aoCapsuleVectors;
    float* aoCapsuleRadii;
    CapsuleSort* aoCapsuleSort;

    int shadowCapsuleCount;
    Vector3* shadowCapsuleStarts;
    Vector3* shadowCapsuleVectors;
    float* shadowCapsuleRadii;
    CapsuleSort* shadowCapsuleSort;

} CapsuleData;

static void CapsuleDataInit(CapsuleData* data)
{
    data->capsuleCount = 0;
    data->capsulePositions = NULL;
    data->capsuleRotations = NULL;
    data->capsuleRadii = NULL;
    data->capsuleHalfLengths = NULL;
    data->capsuleColors = NULL;
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
}

static void CapsuleDataResize(CapsuleData* data, int maxCapsuleCount)
{
    data->capsuleCount = 0;
    data->capsulePositions = realloc(data->capsulePositions, maxCapsuleCount * sizeof(Vector3));
    data->capsuleRotations = realloc(data->capsuleRotations, maxCapsuleCount * sizeof(Quaternion));
    data->capsuleRadii = realloc(data->capsuleRadii, maxCapsuleCount * sizeof(float));
    data->capsuleHalfLengths = realloc(data->capsuleHalfLengths, maxCapsuleCount * sizeof(float));
    data->capsuleColors = realloc(data->capsuleColors, maxCapsuleCount * sizeof(Vector3));
    data->capsuleSort = realloc(data->capsuleSort, maxCapsuleCount * sizeof(CapsuleSort));

    data->aoCapsuleCount = 0;
    data->aoCapsuleStarts = realloc(data->aoCapsuleStarts, maxCapsuleCount * sizeof(Vector3));
    data->aoCapsuleVectors = realloc(data->aoCapsuleVectors, maxCapsuleCount * sizeof(Vector3));
    data->aoCapsuleRadii = realloc(data->aoCapsuleRadii, maxCapsuleCount * sizeof(float));
    data->aoCapsuleSort = realloc(data->aoCapsuleSort, maxCapsuleCount * sizeof(CapsuleSort));

    data->shadowCapsuleCount = 0;
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
    free(data->capsuleSort);
    
    free(data->aoCapsuleStarts);
    free(data->aoCapsuleVectors);
    free(data->aoCapsuleRadii);
    free(data->aoCapsuleSort);
    
    free(data->shadowCapsuleStarts);
    free(data->shadowCapsuleVectors);
    free(data->shadowCapsuleRadii);
    free(data->shadowCapsuleSort);
}

static void CapsuleDataReset(CapsuleData* data)
{
    data->capsuleCount = 0;
    data->aoCapsuleCount = 0;
    data->shadowCapsuleCount = 0;
}

static void CapsuleDataAppendFromTransformData(CapsuleData* data, TransformData* xforms, float maxCapsuleRadius, Color color, bool ignoreEndSite)
{
    for (int i = 0; i < xforms->jointCount; i++)
    {
        int p = xforms->parents[i];
        if (p == -1) { continue; }
        if (ignoreEndSite && xforms->endSite[i]) { continue; }
        
        float capsuleHalfLength = Vector3Length(xforms->localPositions[i]) / 2.0f;
        float capsuleRadius = minf(maxCapsuleRadius, capsuleHalfLength) + (i % 2) * 0.001f;

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
        data->capsuleCount++;
    }
}

static void CapsuleDataUpdateAOCapsulesForGroundSegment(CapsuleData* data, Vector3 groundSegmentPosition)
{
    data->aoCapsuleCount = 0;
  
    for (int i = 0; i < data->capsuleCount; i++)
    {      
        Vector3 capsuleStart = CapsuleStart(data->capsulePositions[i], data->capsuleRotations[i], data->capsuleHalfLengths[i]);
        Vector3 capsuleEnd = CapsuleEnd(data->capsulePositions[i], data->capsuleRotations[i], data->capsuleHalfLengths[i]);
        float capsuleRadius = data->capsuleRadii[i];
        
        Vector3 capsulePoint, groundPoint;
        NearestPointBetweenLineSegmentAndGroundSegment(
            &capsulePoint,
            &groundPoint,
            capsuleStart,
            capsuleEnd,
            (Vector3){ groundSegmentPosition.x - 1.0f, 0.0f, groundSegmentPosition.z - 1.0f },
            (Vector3){ groundSegmentPosition.x + 1.0f, 0.0f, groundSegmentPosition.z + 1.0f });
            
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

static void CapsuleDataUpdateAOCapsulesForCapsule(CapsuleData* data, int capsuleIndex)
{
    Vector3 queryCapsuleStart = CapsuleStart(data->capsulePositions[capsuleIndex], data->capsuleRotations[capsuleIndex], data->capsuleHalfLengths[capsuleIndex]);
    Vector3 queryCapsuleEnd = CapsuleEnd(data->capsulePositions[capsuleIndex], data->capsuleRotations[capsuleIndex], data->capsuleHalfLengths[capsuleIndex]);
    float queryCapsuleRadius = data->capsuleRadii[capsuleIndex];
    
    data->aoCapsuleCount = 0;
  
    for (int i = 0; i < data->capsuleCount; i++)
    {
        if (i == capsuleIndex) { continue; }
        
        Vector3 capsuleStart = CapsuleStart(data->capsulePositions[i], data->capsuleRotations[i], data->capsuleHalfLengths[i]);
        Vector3 capsuleEnd = CapsuleEnd(data->capsulePositions[i], data->capsuleRotations[i], data->capsuleHalfLengths[i]);
        float capsuleRadius = data->capsuleRadii[i];
        
        Vector3 capsulePoint, queryPoint;
        NearestPointBetweenLineSegments(
            &capsulePoint,
            &queryPoint,
            capsuleStart,
            capsuleEnd,
            queryCapsuleStart,
            queryCapsuleEnd);
        
        Vector3 surfaceNormal = Vector3Normalize(Vector3Subtract(capsulePoint, queryPoint));
        Vector3 surfacePoint = Vector3Add(queryPoint, Vector3Scale(surfaceNormal, queryCapsuleRadius));
        
        float capsuleOcclusion = Vector3Distance(capsulePoint, capsulePoint) <= queryCapsuleRadius + capsuleRadius ? 0.0f :
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

static void CapsuleDataUpdateShadowCapsulesForGroundSegment(CapsuleData* data, Vector3 groundSegmentPosition, Vector3 lightDir, float lightConeAngle, int shadowMode)
{
    data->shadowCapsuleCount = 0;
  
    for (int i = 0; i < data->capsuleCount; i++)
    {      
        Vector3 capsulePosition = data->capsulePositions[i];
        Vector3 capsuleStart = CapsuleStart(data->capsulePositions[i], data->capsuleRotations[i], data->capsuleHalfLengths[i]);
        Vector3 capsuleEnd = CapsuleEnd(data->capsulePositions[i], data->capsuleRotations[i], data->capsuleHalfLengths[i]);
        Vector3 capsuleVector = Vector3Subtract(capsuleEnd, capsuleStart);
        float capsuleRadius = data->capsuleRadii[i];
        
        Vector3 rayPoint, groundPlanePoint;
        NearestPointBetweenLineSegmentAndGroundPlane(
            &rayPoint, 
            &groundPlanePoint,
            capsulePosition,
            Vector3Add(capsulePosition, Vector3Scale(lightDir, 10.0f)));
        
        Vector3 groundPoint =
        {
            clampf(groundPlanePoint.x, groundSegmentPosition.x - 1.0f, groundSegmentPosition.x + 1.0f),
            0.0f,
            clampf(groundPlanePoint.z, groundSegmentPosition.z - 1.0f, groundSegmentPosition.z + 1.0f),
        };
        
        float capsuleOcclusion = shadowMode == 2 ? 
                CapsuleDirectionalOcclusion(groundPoint, capsuleStart, capsuleVector, capsuleRadius, lightDir, lightConeAngle) :
                CapsuleRayIntersect(capsuleStart, capsuleVector, capsuleRadius, groundPoint, Vector3Negate(lightDir));
        
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

static void CapsuleDataUpdateShadowCapsulesForCapsule(CapsuleData* data, int capsuleIndex, Vector3 lightDir, float lightConeAngle, int shadowMode)
{
    Vector3 queryCapsuleStart = CapsuleStart(data->capsulePositions[capsuleIndex], data->capsuleRotations[capsuleIndex], data->capsuleHalfLengths[capsuleIndex]);
    Vector3 queryCapsuleEnd = CapsuleEnd(data->capsulePositions[capsuleIndex], data->capsuleRotations[capsuleIndex], data->capsuleHalfLengths[capsuleIndex]);
    float queryCapsuleRadius = data->capsuleRadii[capsuleIndex];
  
    data->shadowCapsuleCount = 0;
  
    for (int i = 0; i < data->capsuleCount; i++)
    {      
        if (i == capsuleIndex) { continue; }
        
        Vector3 capsulePosition = data->capsulePositions[i];
        Vector3 capsuleStart = CapsuleStart(data->capsulePositions[i], data->capsuleRotations[i], data->capsuleHalfLengths[i]);
        Vector3 capsuleEnd = CapsuleEnd(data->capsulePositions[i], data->capsuleRotations[i], data->capsuleHalfLengths[i]);
        Vector3 capsuleVector = Vector3Subtract(capsuleEnd, capsuleStart);
        float capsuleRadius = data->capsuleRadii[i];
        
        // TODO: This isn't right but it will do for now.
        Vector3 rayPoint, queryPoint;
        NearestPointBetweenLineSegments(
            &rayPoint,
            &queryPoint,
            capsulePosition,
            Vector3Add(capsulePosition, Vector3Scale(lightDir, 10.0f)),
            queryCapsuleStart,
            queryCapsuleEnd);
        
        Vector3 surfaceNormal = Vector3Normalize(Vector3Subtract(rayPoint, queryPoint));
        Vector3 surfacePoint = Vector3Add(queryPoint, Vector3Scale(surfaceNormal, queryCapsuleRadius));
        
        float capsuleOcclusion = 
            Vector3Distance(surfacePoint, capsulePosition) <= capsuleRadius ? 0.0f : 
            shadowMode == 2 ? 
                CapsuleDirectionalOcclusion(surfacePoint, capsuleStart, capsuleVector, capsuleRadius, lightDir, lightConeAngle) :
                CapsuleRayIntersect(capsuleStart, capsuleVector, capsuleRadius, surfacePoint, Vector3Negate(lightDir));
        
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

//--------------------------------------
// Shaders
//--------------------------------------

enum
{
    AO_CAPSULES_MAX = 32,
    SHADOW_CAPSULES_MAX = 32,
};

#define GLSL(X) \
  "#version 300 es\n" \
  "#define AO_CAPSULES_MAX 32\n" \
  "#define SHADOW_CAPSULES_MAX 32\n" \
  "#define PI 3.14159265359\n" \
  #X

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
                                                                           
out vec3 fragPosition;                                                     
out vec2 fragTexCoord;                                                     
out vec3 fragNormal;                                                       
                                                                           
vec3 rotate(in vec4 q, vec3 v)                                             
{                                                                          
    vec3 t = 2.0 * cross(q.xyz, v);                                        
    return v + q.w * t + cross(q.xyz, t);                                  
}                                                                          
                                                                           
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
        fragPosition = rotate(capsuleRotation,                             
            CapsuleStretch(vertexPosition,                                
            capsuleHalfLength, capsuleRadius)) + capsulePosition;          
                                                                           
        fragNormal = rotate(capsuleRotation, vertexNormal);                
    }                                                                      
    else                                                                   
    {                                                                      
        fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));         
        fragNormal = normalize(vec3(matNormal * vec4(vertexNormal, 1.0))); 
    }                                                                      
                                                                           
    gl_Position = matProjection * matView * vec4(fragPosition, 1.0);       
}                                          
                                
);

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
                                                                           
uniform int aoCapsuleCount;                                                  
uniform vec3 aoCapsuleStarts[AO_CAPSULES_MAX];                              
uniform vec3 aoCapsuleVectors[AO_CAPSULES_MAX];                             
uniform float aoCapsuleRadii[AO_CAPSULES_MAX];                              
                                                                           
uniform vec3 cameraPosition;                                               
                                                                           
uniform float sunStrength;                                                 
uniform float sunConeAngle;                                                
uniform vec3 sunDir;                                                       
uniform vec3 sunColor;                                                     
uniform float skyStrength;                                                 
uniform vec3 skyColor;                                                     
uniform float ambientStrength;                                             
uniform float groundStrength;                                              
uniform float exposure;                                                    
uniform int shadowMode;

out vec4 finalColor;                                                       
                                                                           
float saturate(in float x)                                                 
{                                                                          
    return clamp(x, 0.0, 1.0);                                             
}                                                                          
                                                                           
float square(in float x)                                                   
{                                                                          
    return x * x;                                                          
}                                                                          
                                                                           
float fast_acos(in float x)                                                
{                                                                          
    float y = abs(x);                                                      
    float p = -0.1565827 * y + 1.570796;                                   
    p *= sqrt(max(1.0 - y, 0.0));                                                    
    return x >= 0.0 ? p : PI - p;                                          
}                                                                          
                                                                           
float fast_positive_acos(in float x)                                       
{                                                                          
    float p = -0.1565827 * x + 1.570796;                                   
    return p * sqrt(max(1.0 - x, 0.0));                                              
}                                                                                                                                             
                                                                           
float fast_positive_atan(in float x)                                       
{                                                                          
    float w = x > 1.0 ? 1.0 / x : x;                                       
    float y = (PI / 4.0)*w - w*(w - 1.0)*(0.2447 + 0.0663*w);              
    return abs(x > 1.0 ? PI / 2.0 - y : y);                                
}                                                                          
                                                                           
float SphereOcclusion(in vec3 pos, in vec3 nor, in vec3 sph, in float rad)
{                                                                          
    vec3 di = sph - pos;                                                  
    float l = max(length(di), rad);
    float nl = dot(nor, normalize(di));                                           
    float h  = l / rad;                                                    
    float h2 = h*h;                                                        
    
    float res;                                                             
    float k2 = 1.0 - h2*nl*nl;                                             
    if (k2 > 0.001)                                                        
    {                                                                      
        res = nl * fast_acos(-nl*sqrt((h2 - 1.0) / (1.0 - nl*nl))) - sqrt(k2*(h2 - 1.0));   
        res = (res / h2 + fast_positive_atan(sqrt(k2 / (h2 - 1.0)))) / PI; 
    }                                                                      
    else                                                                   
    {                                                                      
        res = max(0.0, nl) / h2;                                
    }                                                                      
                                                                           
    return 1.0 - saturate(res);                                            
}                                                                          
                                                                           
float CapsuleOcclusion(                                                   
    in vec3 pos, in vec3 nor,                                              
    in vec3 capStart, in vec3 capVec, in float radius)                   
{                                                                          
    vec3 ba = capVec;                                                     
    vec3 pa = pos - capStart;                 
    float t = saturate(dot(pa, ba) / dot(ba, ba));                          
    
    return SphereOcclusion(pos, nor, capStart + t * ba, radius);         
}                                                                          
                                                                           
float CapsuleSphereIntersectionArea(                                         
    in float cosCap1, in float cosCap2,                                  
    in float cap2, in float cosDist)                                      
{                                                                          
    float r1 = fast_positive_acos(cosCap1);                               
    float r2 = cap2;                                                       
    float d  = fast_acos(cosDist);                                        
                                                                           
    if (min(r1, r2) <= max(r1, r2) - d)                                    
    {                                                                      
        return 1.0 - max(cosCap1, cosCap2);                              
    }                                                                      
    else if (r1 + r2 <= d)                                                 
    {                                                                      
        return 0.0;                                                        
    }                                                                      
                                                                           
    float delta = abs(r1 - r2);                                            
    float x = 1.0 - saturate((d - delta) / max(r1 + r2 - delta, 0.0001));  
    float area = square(x) * (-2.0 * x + 3.0);                             
                                                                           
    return area * (1.0 - max(cosCap1, cosCap2));                         
}                                                                          
                                                                           
float SphereDirectionalOcclusion(                                        
    in vec3 pos, in vec3 sphere, in float radius,                          
    in vec3 coneDir, in float coneAngle)                                 
{                                                                          
    vec3 occluder = sphere - pos;                                          
    float occluder_len2 = dot(occluder, occluder);                         
    vec3 occluder_dir = occluder * inversesqrt(occluder_len2);             
                                                                           
    float cos_phi = dot(occluder_dir, -coneDir);                          
    float cos_theta = sqrt(occluder_len2 /                                 
        (square(radius) + occluder_len2));                                 
    float cos_cone = cos(coneAngle / 2.0) / (1.0 - 1e-8);                                
                                                                           
    return 1.0 - CapsuleSphereIntersectionArea(                              
        cos_theta, cos_cone, coneAngle / 2.0, cos_phi) / (1.0 - cos_cone);
}                                                                          
                                                                           
float CapsuleDirectionalOcclusion(                                       
    in vec3 pos, in vec3 capStart, in vec3 capVec,                       
    in float capRadius, in vec3 coneDir, in float coneAngle)            
{                                                                          
    vec3 ba = capVec;                                                     
    vec3 pa = capStart - pos;                                             
    float a = dot(-coneDir, ba);                                          
    float t = saturate(                                                    
        dot(pa, a * -coneDir - ba) / (dot(ba, ba) - a * a));              
        
    return SphereDirectionalOcclusion(                                   
        pos, capStart + t * ba, capRadius, coneDir, coneAngle);        
}                                                                          

float CapsuleRayIntersect(
    in vec3 capStart, in vec3 capVec, in float capRadius,
    in vec3 rayStart, in vec3 rayDir)
{
    vec3 ba = capVec;
    vec3 oa = rayStart - capStart;

    float baba = dot(ba, ba);
    float bard = dot(ba, rayDir);
    float baoa = dot(ba, oa);
    float rdoa = dot(rayDir, oa);
    float oaoa = dot(oa, oa);

    float a = baba      - bard*bard;
    float b = baba*rdoa - baoa*bard;
    float c = baba*oaoa - baoa*baoa - capRadius*capRadius*baba;
    float h = b*b - a*c;
    
    if (h >= 0.0)
    {
        float t = (-b - sqrt(h))/a;
        float y = baoa + t*bard;
        
        // body
        if (y > 0.0 && y < baba)
        {
            return t > 1e-4 ? 0.0 : 1.0;
        }
        
        // caps
        vec3 oc = y <= 0.0 ? oa : rayStart - (capStart + capVec);
        b = dot(rayDir, oc);
        c = dot(oc, oc) - capRadius*capRadius;
        h = b*b - c;
        
        if (h > 0.0)
        {
            t = -b - sqrt(h);
            return t > 1e-4 ? 0.0 : 1.0;
        }
    }
    
    return 1.0;
}

float checker(in vec2 uv)                                                  
{                                                                          
    vec4 uvDDXY = vec4(dFdx(uv), dFdy(uv));                                 
    vec2 w = vec2(length(uvDDXY.xz), length(uvDDXY.yw));                    
    vec2 i = 2.0*(abs(fract((uv-0.5*w)*0.5)-0.5)-                           
                  abs(fract((uv+0.5*w)*0.5)-0.5))/w;                        
    return 0.5 - 0.5*i.x*i.y;                                                     
}                                                                          
                                                                           
float grid(in vec2 uv, in float lineWidth)                                 
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
                                                                           
vec3 to_gamma(in vec3 col)                                                 
{                                                                          
    return vec3(pow(col.x, 2.2), pow(col.y, 2.2), pow(col.z, 2.2));        
}                                                                          
                                                                           
vec3 from_gamma(in vec3 col)                                               
{                                                                          
    return vec3(                                                           
        pow(col.x, 1.0/2.2), pow(col.y, 1.0/2.2), pow(col.z, 1.0/2.2));    
}                                                                          
                                                                           
vec3 rotate(in vec4 q, vec3 v)                                             
{                                                                          
    vec3 t = 2.0 * cross(q.xyz, v);                                        
    return v + q.w * t + cross(q.xyz, t);                                  
}                                                                          
                                                                           
vec4 inv(in vec4 q)                                                        
{                                                                          
    return vec4(-q.x, -q.y, -q.z, q.w);                                    
}                                                                          
                                                                           
vec3 unrotate(in vec4 q, vec3 v)                                           
{                                                                          
    return rotate(inv(q), v);                                              
}                                                                          
                                                                           
vec2 CapsuleUVs(                                                          
    in vec3 pos, in vec3 capPos,                                          
    in vec4 capRot, in float capHalfLength,                                 
    in float capRadius, in vec2 scale)                                    
{                                                                          
    vec3 loc = unrotate(capRot, pos - capPos);                           
                                                                           
    vec2 limit = vec2(                                                     
        2.0 * capHalfLength + 2.0 * capRadius,                              
        PI * capRadius);                                                  
                                                                           
    vec2 repeat = max(round(scale * limit), 1.0);                          
                                                                           
    return (repeat / limit) * vec2(                                        
        loc.x, capRadius * atan(loc.z, loc.y));                           
}                                                                          
                                                                           
vec3 CapsuleNormal(                                                       
    in vec3 pos, in vec3 capStart,                                        
    in vec3 capVec)                                      
{                                                                          
    vec3 ba = capVec;                                                     
    vec3 pa = pos - capStart;                                             
    float h = saturate(dot(pa, ba) / dot(ba, ba));
    
    return normalize(pa - h*ba);                                          
}                                                                                            

void main()                                                                
{                                                                          
    vec3 pos = fragPosition;                                               
    vec3 nor = fragNormal;                                                 
    vec2 uvs = fragTexCoord;                                               
                                                                           
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
                                                                           
    vec3 eye_dir = normalize(fragPosition - cameraPosition);               
                                                                           
    vec3 light_sunColor = from_gamma(sunColor);                           
    vec3 light_sun_half = normalize(sunDir + eye_dir);                     
                                                                           
    vec3 light_skyColor = from_gamma(skyColor);                           
    vec3 skyDir = vec3(0.0, -1.0, 0.0);                                    
    vec3 light_sky_half = normalize(skyDir + eye_dir);                     
                                                                           
    float sunShadow = 1.0;
    
    if (shadowMode == 2)
    {
        for (int i = 0; i < shadowCapsuleCount; i++)                             
        {                                                                      
            sunShadow = min(sunShadow, CapsuleDirectionalOcclusion(                       
                pos,                                                           
                shadowCapsuleStarts[i],                                        
                shadowCapsuleVectors[i],                                       
                shadowCapsuleRadii[i],                                         
                sunDir,                                                        
                sunConeAngle));                                        
        }                  
    }
    else if (shadowMode == 1)
    {
        for (int i = 0; i < shadowCapsuleCount; i++)                             
        {
            sunShadow = min(sunShadow, CapsuleRayIntersect(                       
                shadowCapsuleStarts[i],                                        
                shadowCapsuleVectors[i],                                       
                shadowCapsuleRadii[i],
                pos,
                -sunDir));
        }  
    }
    
    float ambShadow = 1.0;                                                
    for (int i = 0; i < aoCapsuleCount; i++)                                 
    {                                                                      
        ambShadow = min(ambShadow, CapsuleOcclusion(                                   
            pos, nor,                                                      
            aoCapsuleStarts[i],                                            
            aoCapsuleVectors[i],                                           
            aoCapsuleRadii[i]));                                            
    }                                                                      
                                                                           
    float grid_fine = grid(20.0 * uvs, 0.025);                             
    float grid_coarse = grid(2.0 * uvs, 0.02);                             
    float check = checker(2.0 * uvs);                                      
                                                                           
    vec3 color = from_gamma(objectColor);                                  
    vec3 albedo = color * mix(mix(mix(                                     
        0.9, 0.95, check),                                                 
        0.85, grid_fine),                                                  
        1.0, grid_coarse);                                                 
    float specularity =                                                    
       objectSpecularity * mix(mix(0.0, 0.75, check), 1.0, grid_coarse);   
                                                                           
    float sun_factor_diff = max(dot(nor, -sunDir), 0.0);                   
    float sun_factor_spec = specularity *                                  
        ((objectGlossiness+2.0) / (8.0 * PI)) *                            
        max(pow(dot(nor, light_sun_half), objectGlossiness), 0.0);         
                                                                           
    float sky_factor_diff = max(dot(nor, -skyDir), 0.0);                   
    float sky_factor_spec = specularity *                                  
        ((objectGlossiness+2.0) / (8.0 * PI)) *                            
        max(pow(dot(nor, light_sky_half), objectGlossiness), 0.0);         
                                                                           
    float ground_factor_diff = max(dot(nor, skyDir), 0.0);                 
                                                                           
    vec3 ambient =                                                         
        ambShadow * ambientStrength * light_skyColor * albedo;           
                                                                           
    vec3 diffuse =                                                         
        sunShadow * sunStrength * light_sunColor *                       
                     albedo * sun_factor_diff +                            
                     groundStrength * light_skyColor *                    
                     albedo * ground_factor_diff;                          
                     skyStrength * light_skyColor *                       
                     albedo * sky_factor_diff;                             
                                                                           
    float specular = sunShadow * sunStrength * sun_factor_spec +          
                                  skyStrength * sky_factor_spec;           
                                                                           
    vec3 final = diffuse + ambient + specular;                             
                                                                           
    finalColor = vec4(to_gamma(exposure * final), objectOpacity);                    
}    
                                                                      
);

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
    
    int aoCapsuleCount;
    int aoCapsuleStarts;
    int aoCapsuleVectors;        
    int aoCapsuleRadii;
    
    int cameraPosition;
    
    int objectColor;
    int objectSpecularity;
    int objectGlossiness;
    int objectOpacity;
    
    int sunStrength;
    int sunConeAngle;
    int sunDir;
    int sunColor;
    int skyStrength;
    int skyColor;
    int ambientStrength;
    int groundStrength;
    int shadowMode;
    
    int exposure;
  
} ShaderUniforms;

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
    
    uniforms->aoCapsuleCount = GetShaderLocation(shader, "aoCapsuleCount");
    uniforms->aoCapsuleStarts =  GetShaderLocation(shader, "aoCapsuleStarts");
    uniforms->aoCapsuleVectors =  GetShaderLocation(shader, "aoCapsuleVectors");        
    uniforms->aoCapsuleRadii =  GetShaderLocation(shader, "aoCapsuleRadii");
    
    uniforms->cameraPosition = GetShaderLocation(shader, "cameraPosition");
    
    uniforms->objectColor = GetShaderLocation(shader, "objectColor");
    uniforms->objectSpecularity = GetShaderLocation(shader, "objectSpecularity");
    uniforms->objectGlossiness = GetShaderLocation(shader, "objectGlossiness");
    uniforms->objectOpacity = GetShaderLocation(shader, "objectOpacity");
    
    uniforms->sunStrength = GetShaderLocation(shader, "sunStrength");
    uniforms->sunConeAngle = GetShaderLocation(shader, "sunConeAngle");
    uniforms->sunDir = GetShaderLocation(shader, "sunDir");
    uniforms->sunColor = GetShaderLocation(shader, "sunColor");
    uniforms->skyStrength = GetShaderLocation(shader, "skyStrength");
    uniforms->skyColor = GetShaderLocation(shader, "skyColor");
    uniforms->ambientStrength = GetShaderLocation(shader, "ambientStrength");
    uniforms->groundStrength = GetShaderLocation(shader, "groundStrength");
    uniforms->shadowMode = GetShaderLocation(shader, "shadowMode");
    
    uniforms->exposure = GetShaderLocation(shader, "exposure");
}

//--------------------------------------
// Models
//--------------------------------------

#define OBJDATA(X) #X

static const char* capsuleOBJ =
"v 0.82165808 -0.82165808 -1.0579772e-18     \n"
"v 0.82165808 -0.58100000 0.58100000         \n"
"v 0.82165808 8.7595780e-17 0.82165808       \n"
"v 0.82165808 0.58100000 0.58100000          \n"
"v 0.82165808 0.82165808 9.9566116e-17       \n"
"v 0.82165808 0.58100000 -0.58100000         \n"
"v 0.82165808 2.8884397e-16 -0.82165808      \n"
"v 0.82165808 -0.58100000 -0.58100000        \n"
"v -0.82165808 -0.82165808 -1.0579772e-18    \n"
"v -0.82165808 -0.58100000 0.58100000        \n"
"v -0.82165808 -1.3028313e-17 0.82165808     \n"
"v -0.82165808 0.58100000 0.58100000         \n"
"v -0.82165808 0.82165808 9.9566116e-17      \n"
"v -0.82165808 0.58100000 -0.58100000        \n"
"v -0.82165808 1.8821987e-16 -0.82165808     \n"
"v -0.82165808 -0.58100000 -0.58100000       \n"
"v 1.16200000 1.5874776e-16 -1.0579772e-18   \n"
"v -1.16200000 1.6443801e-17 -1.0579772e-18  \n"
"v -9.1030792e-3 -1.15822938 -1.0579772e-18  \n"
"v 9.1030792e-3 -1.15822938 -1.0579772e-18   \n"
"v 9.1030792e-3 -0.81899185 0.81899185       \n"
"v -9.1030792e-3 -0.81899185 0.81899185      \n"
"v 9.1030792e-3 1.7232088e-17 1.15822938     \n"
"v -9.1030792e-3 1.6117282e-17 1.15822938    \n"
"v 9.1030792e-3 0.81899185 0.81899185        \n"
"v -9.1030792e-3 0.81899185 0.81899185       \n"
"v 9.1030792e-3 1.15822938 1.4078421e-16     \n"
"v -9.1030792e-3 1.15822938 1.4078421e-16    \n"
"v 9.1030792e-3 0.81899185 -0.81899185       \n"
"v -9.1030792e-3 0.81899185 -0.81899185      \n"
"v 9.1030792e-3 3.0091647e-16 -1.15822938    \n"
"v -9.1030792e-3 2.9980166e-16 -1.15822938   \n"
"v 9.1030792e-3 -0.81899185 -0.81899185      \n"
"v -9.1030792e-3 -0.81899185 -0.81899185     \n"
"vn 0.71524683 -0.69887193 -2.5012597e-16    \n"
"vn 0.61185516 -0.55930013 0.55930013        \n"
"vn 0.71524683 0.0000000e+0 0.69887193       \n"
"vn 0.61185516 0.55930013 0.55930013         \n"
"vn 0.71524683 0.69887193 1.5632873e-17      \n"
"vn 0.61185516 0.55930013 -0.55930013        \n"
"vn 0.71524683 6.2531494e-17 -0.69887193     \n"
"vn 0.61185516 -0.55930013 -0.55930013       \n"
"vn -0.71524683 -0.69887193 -2.5012597e-16   \n"
"vn -0.61185516 -0.55930013 0.55930013       \n"
"vn -0.71524683 0.0000000e+0 0.69887193      \n"
"vn -0.61185516 0.55930013 0.55930013        \n"
"vn -0.71524683 0.69887193 4.6898620e-17     \n"
"vn -0.61185516 0.55930013 -0.55930013       \n"
"vn -0.71524683 4.6898620e-17 -0.69887193    \n"
"vn -0.61185516 -0.55930013 -0.55930013      \n"
"vn 1.00000000 1.5208752e-17 -2.6615316e-17  \n"
"vn -1.00000000 -1.5208752e-17 2.2813128e-17 \n"
"vn -0.19614758 -0.98057439 -2.2848712e-16   \n"
"vn 0.26047011 -0.96548191 -2.4273177e-16    \n"
"vn 0.13072302 -0.70103905 0.70103905        \n"
"vn -0.19614758 -0.69337080 0.69337080       \n"
"vn 0.22349711 5.9825845e-2 0.97286685       \n"
"vn -0.22349711 -5.9825845e-2 0.97286685     \n"
"vn 0.15641931 0.75510180 0.63667438         \n"
"vn -0.15641931 0.63667438 0.75510180        \n"
"vn 0.22349711 0.97286685 -5.9825845e-2      \n"
"vn -0.22349711 0.97286685 5.9825845e-2      \n"
"vn 0.15641931 0.63667438 -0.75510180        \n"
"vn -0.15641931 0.75510180 -0.63667438       \n"
"vn 0.22349711 -5.9825845e-2 -0.97286685     \n"
"vn -0.22349711 5.9825845e-2 -0.97286685     \n"
"vn 0.15641931 -0.75510180 -0.63667438       \n"
"vn -0.15641931 -0.63667438 -0.75510180      \n"
"f 1//1 17//17 2//2                          \n"
"f 1//1 20//20 8//8                          \n"
"f 2//2 17//17 3//3                          \n"
"f 2//2 20//20 1//1                          \n"
"f 2//2 23//23 21//21                        \n"
"f 3//3 17//17 4//4                          \n"
"f 3//3 23//23 2//2                          \n"
"f 4//4 17//17 5//5                          \n"
"f 4//4 23//23 3//3                          \n"
"f 4//4 27//27 25//25                        \n"
"f 5//5 17//17 6//6                          \n"
"f 5//5 27//27 4//4                          \n"
"f 6//6 17//17 7//7                          \n"
"f 6//6 27//27 5//5                          \n"
"f 6//6 31//31 29//29                        \n"
"f 7//7 17//17 8//8                          \n"
"f 7//7 31//31 6//6                          \n"
"f 8//8 17//17 1//1                          \n"
"f 8//8 20//20 33//33                        \n"
"f 8//8 31//31 7//7                          \n"
"f 9//9 18//18 16//16                        \n"
"f 9//9 19//19 10//10                        \n"
"f 10//10 18//18 9//9                        \n"
"f 10//10 19//19 22//22                      \n"
"f 10//10 24//24 11//11                      \n"
"f 11//11 18//18 10//10                      \n"
"f 11//11 24//24 12//12                      \n"
"f 12//12 18//18 11//11                      \n"
"f 12//12 24//24 26//26                      \n"
"f 12//12 28//28 13//13                      \n"
"f 13//13 18//18 12//12                      \n"
"f 13//13 28//28 14//14                      \n"
"f 14//14 18//18 13//13                      \n"
"f 14//14 28//28 30//30                      \n"
"f 14//14 32//32 15//15                      \n"
"f 15//15 18//18 14//14                      \n"
"f 15//15 32//32 16//16                      \n"
"f 16//16 18//18 15//15                      \n"
"f 16//16 19//19 9//9                        \n"
"f 16//16 32//32 34//34                      \n"
"f 19//19 33//33 20//20                      \n"
"f 20//20 21//21 19//19                      \n"
"f 21//21 20//20 2//2                        \n"
"f 21//21 24//24 22//22                      \n"
"f 22//22 19//19 21//21                      \n"
"f 22//22 24//24 10//10                      \n"
"f 23//23 26//26 24//24                      \n"
"f 24//24 21//21 23//23                      \n"
"f 25//25 23//23 4//4                        \n"
"f 25//25 28//28 26//26                      \n"
"f 26//26 23//23 25//25                      \n"
"f 26//26 28//28 12//12                      \n"
"f 27//27 30//30 28//28                      \n"
"f 28//28 25//25 27//27                      \n"
"f 29//29 27//27 6//6                        \n"
"f 29//29 32//32 30//30                      \n"
"f 30//30 27//27 29//29                      \n"
"f 30//30 32//32 14//14                      \n"
"f 31//31 34//34 32//32                      \n"
"f 32//32 29//29 31//31                      \n"
"f 33//33 19//19 34//34                      \n"
"f 33//33 31//31 8//8                        \n"
"f 34//34 19//19 16//16                      \n"
"f 34//34 31//31 33//33                      \n"
"                                            \n";

#undef TINYOBJ_LOADER_C_IMPLEMENTATION
#include "external/tinyobj_loader_c.h"      // OBJ/MTL file formats loading

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
        int ret = tinyobj_parse_obj(&attrib, &meshes, &meshCount, &materials, &materialCount, fileText, dataSize, flags);
        assert(ret == 0);

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

//--------------------------------------
// Command Line Args
//--------------------------------------

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

static inline int ArgBool(int argc, char** argv, const char* name, bool defaultValue)
{ 
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }
    if (strcmp(value, "true") == 0) { printf("INFO: Parsed option '%s' as '%s'\n", name, value); return true; }
    if (strcmp(value, "false") == 0) { printf("INFO: Parsed option '%s' as '%s'\n", name, value); return false; }
    
    printf("ERROR: Could not parse value '%s' given for option '%s' as bool\n", value, name);
    return defaultValue;
}

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

static inline const char* ArgStr(int argc, char** argv, const char* name, const char* defaultValue)
{ 
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }

    printf("INFO: Parsed option '%s' as '%s'\n", name, value);
    return value;
}

static inline Color ArgColor(int argc, char** argv, const char* name, Color defaultValue)
{ 
    char* value = ArgFind(argc, argv, name);
    if (!value) { return defaultValue; }

    int cx, cy, cz;
    if (sscanf(value, "%i,%i,%i", &cx, &cy, &cz) == 3)
    {
        printf("INFO: Parsed option '%s' as '%s'\n", name, value);
        return (Color){ clamp(cx, 0, 255), clamp(cy, 0, 255), clamp(cz, 0, 255) };
    }

    printf("ERROR: Could not parse value '%s' given for option '%s' as color\n", value, name);
    return defaultValue;
}

//--------------------------------------
// Rendering
//--------------------------------------

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
    bool drawEndSites;
    bool drawFPS;
    bool drawUI;
    
    int shadowMode;
    
    float maxCapsuleRadius;
    float capsuleOpacity;
  
} RenderSettings;

void RenderSettingsInit(RenderSettings* settings, int argc, char** argv)
{
    settings->backgroundColor = ArgColor(argc, argv, "backgroundColor", WHITE);
    
    settings->sunLightConeAngle = ArgFloat(argc, argv, "sunLightConeAngle", 0.2f);
    settings->sunLightStrength = ArgFloat(argc, argv, "sunLightStrength", 0.25f);
    settings->sunAzimuth = ArgFloat(argc, argv, "sunAzimuth", 0.0f);
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
    settings->drawEndSites = ArgBool(argc, argv, "drawEndSites", true);
    settings->drawFPS = ArgBool(argc, argv, "drawFPS", false);
    settings->drawUI = ArgBool(argc, argv, "drawUI", true);
    
    settings->shadowMode = ArgEnum(argc, argv, "shadowMode", 3, (const char*[]){"none", "hard", "soft"}, 2);
    
    settings->maxCapsuleRadius = ArgFloat(argc, argv, "maxCapsuleRadius", 0.04f);
    settings->capsuleOpacity = ArgFloat(argc, argv, "capsuleOpacity", 1.0f);
}

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

//--------------------------------------
// GUI
//--------------------------------------

static inline void GuiRenderSettings(RenderSettings* renderSettings, int screenWidth, int screenHeight)
{
    GuiGroupBox((Rectangle){ screenWidth - 260, 10, 240, 490 }, "Rendering");
    
    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 20, 100, 20 }, 
        "Exposure", 
        TextFormat("%5.2f", renderSettings->exposure),
        &renderSettings->exposure, 
        0.0f, 3.0f);
    
    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 50, 100, 20 }, 
        "Sun Light", 
        TextFormat("%5.2f", renderSettings->sunLightStrength),
        &renderSettings->sunLightStrength, 
        0.0f, 1.0f);
    
    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 80, 100, 20 }, 
        "Sun Softness", 
        TextFormat("%5.2f", renderSettings->sunLightConeAngle),
        &renderSettings->sunLightConeAngle, 
        0.02f, PIf / 4.0f);
    
    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 110, 100, 20 }, 
        "Sky Light", 
        TextFormat("%5.2f", renderSettings->skyLightStrength),
        &renderSettings->skyLightStrength, 
        0.0f, 1.0f);
    
    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 140, 100, 20 }, 
        "Ambient Light", 
        TextFormat("%5.2f", renderSettings->ambientLightStrength),
        &renderSettings->ambientLightStrength, 
        0.0f, 2.0f);
    
    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 170, 100, 20 }, 
        "Ground Light", 
        TextFormat("%5.2f", renderSettings->groundLightStrength),
        &renderSettings->groundLightStrength, 
        0.0f, 0.5f);
    
    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 200, 100, 20 }, 
        "Sun Azimuth", 
        TextFormat("%5.2f", renderSettings->sunAzimuth),
        &renderSettings->sunAzimuth, 
        -PIf, PIf);
    
    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 230, 100, 20 }, 
        "Sun Altitude", 
        TextFormat("%5.2f", renderSettings->sunAltitude),
        &renderSettings->sunAltitude, 
        0.0f, 0.49f * PIf);
    
    GuiCheckBox((Rectangle){ screenWidth - 250, 260, 20, 20 }, "Draw Origin", &renderSettings->drawOrigin);
    GuiCheckBox((Rectangle){ screenWidth - 130, 260, 20, 20 }, "Draw Grid", &renderSettings->drawGrid);
    GuiCheckBox((Rectangle){ screenWidth - 250, 290, 20, 20 }, "Draw Checker", &renderSettings->drawChecker);
    GuiCheckBox((Rectangle){ screenWidth - 130, 290, 20, 20 }, "Draw Capsules", &renderSettings->drawCapsules);
    GuiCheckBox((Rectangle){ screenWidth - 250, 320, 20, 20 }, "Draw Wireframes", &renderSettings->drawWireframes);
    GuiCheckBox((Rectangle){ screenWidth - 130, 320, 20, 20 }, "Draw Skeleton", &renderSettings->drawSkeleton);
    GuiCheckBox((Rectangle){ screenWidth - 250, 350, 20, 20 }, "Draw Transforms", &renderSettings->drawTransforms);
    GuiCheckBox((Rectangle){ screenWidth - 130, 350, 20, 20 }, "Draw AO", &renderSettings->drawAO);
    GuiCheckBox((Rectangle){ screenWidth - 250, 380, 20, 20 }, "Draw End Sites", &renderSettings->drawEndSites);
    GuiCheckBox((Rectangle){ screenWidth - 130, 380, 20, 20 }, "Draw FPS", &renderSettings->drawFPS);
    GuiLabel((Rectangle){ screenWidth - 240, 410, 100, 20 }, "Shadows");
    GuiComboBox((Rectangle){ screenWidth - 180, 410, 150, 20 }, "None;Hard;Soft", &renderSettings->shadowMode);

    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 440, 100, 20 }, 
        "Capsule Radius", 
        TextFormat("%5.2f", renderSettings->maxCapsuleRadius),
        &renderSettings->maxCapsuleRadius, 
        0.01f, 0.1f);  

    GuiSliderBar(
        (Rectangle){ screenWidth - 160, 470, 100, 20 }, 
        "Capsule Opacity", 
        TextFormat("%5.2f", renderSettings->capsuleOpacity),
        &renderSettings->capsuleOpacity, 
        0.0f, 1.0f);  
}

//--------------------------------------
// Loading
//--------------------------------------

enum
{
    MAX_BVH_FILES = 6,
};

static bool LoadBVHFile(
    const char* path, 
    char* errMsg,
    int errMsgSize,
    BVHData* bvhData, 
    TransformData* xformData,
    TransformData* xformTmp0,
    TransformData* xformTmp1,
    TransformData* xformTmp2,
    TransformData* xformTmp3,
    char* bvhName, 
    int bvhNameSize, 
    float* bvhScale,
    float* bvhAutoScale)
{
    printf("INFO: Loading '%s'\n", path);

    if (BVHDataLoad(bvhData, path, errMsg, errMsgSize))
    {
        TransformDataResize(xformData, bvhData);
        TransformDataResize(xformTmp0, bvhData);
        TransformDataResize(xformTmp1, bvhData);
        TransformDataResize(xformTmp2, bvhData);
        TransformDataResize(xformTmp3, bvhData);
      
        const char* filename = path;
        while (strchr(filename, '/')) { filename = strchr(filename, '/') + 1; }
        while (strchr(filename, '\\')) { filename = strchr(filename, '\\') + 1; }
      
        snprintf(bvhName, bvhNameSize, "%s", filename);
        *bvhScale = 1.0f;
        
        // Set Window Title
        
        char windowTitle[512];
        snprintf(windowTitle, 512, "%s - BVHView", path);
        SetWindowTitle(windowTitle);
        
        // Auto-Scaling and unit detection
        
        if (bvhData->frameCount > 0)
        {
            TransformDataSampleFrame(xformData, bvhData, 0, 1.0f);
            TransformDataForwardKinematics(xformData);
            
            float height = 1e-8f;
            for (int j = 0; j < xformData->jointCount; j++)
            {
                height = maxf(height, xformData->globalPositions[j].y);
            } 
            
            *bvhScale = height > 10.0f ? 0.01f : 1.0f;
            *bvhAutoScale = 1.8 / height;
        }
        else
        {
            *bvhAutoScale = 1.0f;
        }
        
        return true;
    }
    else
    {
        printf("INFO: Failed to Load '%s'\n", path);

        return false;
    }
}

static inline void UpdateCapsuleBuffersForBVHs(BVHData bvhData[], int bvhCount, CapsuleData* capsuleData)
{
    int totalJointCount = 0;
    for (int i = 0; i < bvhCount; i++)
    {
        totalJointCount += bvhData[i].jointCount;
    }
    
    CapsuleDataResize(capsuleData, totalJointCount);
}

static inline void UpdateScrubberForBVHs(
    int* frameLimit, 
    float* timeLimit, 
    int* frameMax,
    int* frameMaxSelect,
    float* timeMax,
    BVHData bvhData[], 
    int bvhCount,
    int bvhActive)
{
    for (int i = 0; i < bvhCount; i++)
    {
        int bvhFrameLimit = bvhData[i].frameCount - 1;
        float bvhTimeLimit = bvhFrameLimit * bvhData[i].frameTime;
        *frameLimit = max(*frameLimit, bvhFrameLimit);
        *timeLimit = maxf(*timeLimit, bvhTimeLimit);
    }
    
    if (bvhCount > 0)
    {
        *frameMax = bvhData[bvhActive].frameCount - 1;
        *frameMaxSelect = bvhData[bvhActive].frameCount - 1;
        *timeMax = (bvhData[bvhActive].frameCount - 1) * bvhData[bvhActive].frameTime;
    }
}

//--------------------------------------
// Main
//--------------------------------------

int main(int argc, char** argv)
{
    // Init Window
    
    const int screenWidth = ArgInt(argc, argv, "screenWidth", 1280);
    const int screenHeight = ArgInt(argc, argv, "screenHeight", 720);
    
    SetConfigFlags(FLAG_VSYNC_HINT);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(screenWidth, screenHeight, "BVHView");
    SetTargetFPS(60);
        
    // Camera

    Camera3D camera = { 0 };
    camera.position = (Vector3){ 2.0f, 3.0f, 5.0f };
    camera.target = (Vector3){ -0.5f, 1.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = ArgFloat(argc, argv, "cameraFOV", 45.0f);
    camera.projection = CAMERA_PERSPECTIVE;
    
    float cameraAzimuth = ArgFloat(argc, argv, "cameraAzimuth", 0.0f);
    float cameraAltitude = ArgFloat(argc, argv, "cameraAltitude", 0.4f);
    float cameraDistance = ArgFloat(argc, argv, "cameraDistance", 4.0f);
    bool cameraTrack = ArgBool(argc, argv, "cameraTrack", true);
    int cameraTrackBone = ArgInt(argc, argv, "cameraTrackBone", 0);
    
    // Shader
    
    Shader shader = LoadShaderFromMemory(shaderVS, shaderFS);
    ShaderUniforms uniforms;
    ShaderUniformsInit(&uniforms, shader);
    
    // Meshes
    
    Mesh groundPlaneMesh = GenMeshPlane(2.0f, 2.0f, 1, 1);
    Model groundPlaneModel = LoadModelFromMesh(groundPlaneMesh);
    groundPlaneModel.materials[0].shader = shader;
    
    Model capsuleModel = LoadOBJFromMemory(capsuleOBJ);
    capsuleModel.materials[0].shader = shader;

    // BVH, Transform, and Capsule Data
    
    int bvhCount = 0;
    int bvhActive = 0;

    BVHData bvhData[MAX_BVH_FILES];
    float bvhScales[MAX_BVH_FILES];
    char bvhNames[MAX_BVH_FILES][128];
    float bvhAutoScales[MAX_BVH_FILES];
    Color bvhColors[MAX_BVH_FILES] = 
    {
        ORANGE,
        VIOLET,
        PINK,
        PURPLE,
        MAGENTA,
        GREEN,
    };
    
    TransformData xformData[MAX_BVH_FILES];
    TransformData xformTmp0[MAX_BVH_FILES];
    TransformData xformTmp1[MAX_BVH_FILES];
    TransformData xformTmp2[MAX_BVH_FILES];
    TransformData xformTmp3[MAX_BVH_FILES];
    
    for (int i = 0; i < MAX_BVH_FILES; i++)
    {
        BVHDataInit(&bvhData[i]);
        bvhScales[i] = 1.0f;
        bvhNames[i][0] = '\0';
        bvhAutoScales[i] = 1.0f;
        TransformDataInit(&xformData[i]);
        TransformDataInit(&xformTmp0[i]);
        TransformDataInit(&xformTmp1[i]);
        TransformDataInit(&xformTmp2[i]);
        TransformDataInit(&xformTmp3[i]);    
    }
    
    CapsuleData capsuleData;
    CapsuleDataInit(&capsuleData);
    
    // Playback and Scrubber Settings
    
    bool playing = true;
    bool looping = false;
    float playTime = 0.0f;
    float playSpeed = 1.0f;
    bool frameSnap = true;
    int sampleMode = 1;
    float timeLimit = 0.0f;
    int frameLimit = 0;
    int frame = 0;
    int frameMin = 0;
    int frameMax = 0;
    int frameMinSelect = 0;
    int frameMaxSelect = 0;
    bool frameMinEdit = false;
    bool frameMaxEdit = false;
    float timeMin = 0.0f;
    float timeMax = 0.0f;
    
    
    // Render Settings
    
    RenderSettings renderSettings;
    RenderSettingsInit(&renderSettings, argc, argv);
    
    // File Dialog

    GuiWindowFileDialogState fileDialogState = InitGuiWindowFileDialog(GetWorkingDirectory());
    
    // Load File(s)
    
    char errMsg[512];
    errMsg[0] = '\0';
    
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-') { continue; }
      
        if (bvhCount == MAX_BVH_FILES)
        {
            snprintf(errMsg, 512, "Error: Maximum number of BVH files loaded (%i)", MAX_BVH_FILES);
            break;
        }
        
        if (LoadBVHFile(
            argv[i], 
            errMsg,
            512,
            &bvhData[bvhCount], 
            &xformData[bvhCount],
            &xformTmp0[bvhCount],
            &xformTmp1[bvhCount],
            &xformTmp2[bvhCount],
            &xformTmp3[bvhCount],
            bvhNames[bvhCount], 
            128, 
            &bvhScales[bvhCount],
            &bvhAutoScales[bvhCount]))
        {
            bvhCount++;
        }
    }
    
    // Allocate capsule data
    
    UpdateCapsuleBuffersForBVHs(bvhData, bvhCount, &capsuleData);
    
    // Init Scrubber Frame Max
    
    UpdateScrubberForBVHs(
        &frameLimit, 
        &timeLimit, 
        &frameMax,
        &frameMaxSelect,
        &timeMax,
        bvhData, 
        bvhCount,
        bvhActive);
            
    // Go
    
    while (!WindowShouldClose())
    {        
        // File Dialog
        
        if (fileDialogState.SelectFilePressed)
        {
            if (bvhCount == MAX_BVH_FILES)
            {
                snprintf(errMsg, 512, "Error: Maximum number of BVH files loaded (%i)", MAX_BVH_FILES);
            }
            else
            {
                if (IsFileExtension(fileDialogState.fileNameText, ".bvh"))
                {
                    char fileNameToLoad[512];
                    snprintf(fileNameToLoad, 512, "%s/%s", fileDialogState.dirPathText, fileDialogState.fileNameText);
                    
                    if (LoadBVHFile(
                        fileNameToLoad, 
                        errMsg,
                        512,
                        &bvhData[bvhCount], 
                        &xformData[bvhCount],
                        &xformTmp0[bvhCount],
                        &xformTmp1[bvhCount],
                        &xformTmp2[bvhCount],
                        &xformTmp3[bvhCount],
                        bvhNames[bvhCount], 
                        128, 
                        &bvhScales[bvhCount],
                        &bvhAutoScales[bvhCount]))
                    {
                        bvhCount++;
                        bvhActive = bvhCount - 1;
                    }
                    
                    UpdateCapsuleBuffersForBVHs(bvhData, bvhCount, &capsuleData);
                    UpdateScrubberForBVHs(
                        &frameLimit, 
                        &timeLimit, 
                        &frameMax,
                        &frameMaxSelect,
                        &timeMax,
                        bvhData, 
                        bvhCount,
                        bvhActive);
                }
                else
                {
                    snprintf(errMsg, 512, "Error: File '%s' is not a BVH file.", fileDialogState.fileNameText);
                }
            }
            
            fileDialogState.SelectFilePressed = false;
        }
        
        // Drag and Drop Files
        
        if (IsFileDropped())
        {
            FilePathList droppedFiles = LoadDroppedFiles();
            
            for (int i = 0; i < droppedFiles.count; i++)
            {
                if (bvhCount == MAX_BVH_FILES)
                {
                    snprintf(errMsg, 512, "Error: Maximum number of BVH files loaded (%i)", MAX_BVH_FILES);
                    break;
                }
              
                if (LoadBVHFile(
                    droppedFiles.paths[i], 
                    errMsg,
                    512,
                    &bvhData[bvhCount], 
                    &xformData[bvhCount],
                    &xformTmp0[bvhCount],
                    &xformTmp1[bvhCount],
                    &xformTmp2[bvhCount],
                    &xformTmp3[bvhCount],
                    bvhNames[bvhCount], 
                    128, 
                    &bvhScales[bvhCount],
                    &bvhAutoScales[bvhCount]))
                {
                    bvhCount++;
                    bvhActive = bvhCount - 1;
                }
            }

            UpdateCapsuleBuffersForBVHs(bvhData, bvhCount, &capsuleData);
            UpdateScrubberForBVHs(
                &frameLimit, 
                &timeLimit, 
                &frameMax,
                &frameMaxSelect,
                &timeMax,
                bvhData, 
                bvhCount,
                bvhActive);

            UnloadDroppedFiles(droppedFiles);
        }
  
        // Move time forward
        
        if (playing)
        {
            playTime += playSpeed * GetFrameTime();
            
            if (playTime >= timeMax)
            {
                playTime = looping ? fmod(playTime, timeMax) + timeMin : timeMax;
            }
        }
        
        // Sample Animation Data
        
        for (int i = 0; i < bvhCount; i++)
        {
            if (sampleMode == 0)
            {
                TransformDataSampleFrameNearest(
                    &xformData[i],
                    &bvhData[i],
                    playTime,
                    bvhScales[i]);
            }
            else if (sampleMode == 1)
            {
                TransformDataSampleFrameLinear(
                    &xformData[i],
                    &xformTmp0[i],
                    &xformTmp1[i],
                    &bvhData[i],
                    playTime,
                    bvhScales[i]);
            }
            else
            {
                TransformDataSampleFrameCubic(
                    &xformData[i],
                    &xformTmp0[i],
                    &xformTmp1[i],
                    &xformTmp2[i],
                    &xformTmp3[i],
                    &bvhData[i],
                    playTime,
                    bvhScales[i]);
            }
                
            TransformDataForwardKinematics(&xformData[i]);
        }
        
        // Update Camera
        
        Vector3 cameraTarget = (Vector3){ 0.0f, 1.0f, 0.0f };

        if (bvhCount > 0 && cameraTrack && cameraTrackBone < xformData[bvhActive].jointCount)
        {
            cameraTarget = xformData[bvhActive].globalPositions[cameraTrackBone];
        }

        if (!fileDialogState.windowActive)
        {
            OrbitCameraUpdate(
                &camera, 
                &cameraAzimuth,
                &cameraAltitude,
                &cameraDistance,
                cameraTarget,
                (IsKeyDown(KEY_LEFT_CONTROL) && IsMouseButtonDown(0)) ? GetMouseDelta().x : 0.0f,
                (IsKeyDown(KEY_LEFT_CONTROL) && IsMouseButtonDown(0)) ? GetMouseDelta().y : 0.0f,
                GetFrameTime());
        }
        
        // Create Capsules
        
        CapsuleDataReset(&capsuleData);
        for (int i = 0; i < bvhCount; i++)
        {
            CapsuleDataAppendFromTransformData(
                &capsuleData, 
                &xformData[i], 
                renderSettings.maxCapsuleRadius, 
                bvhColors[i], 
                !renderSettings.drawEndSites);
        }
        
        // Render
        
        BeginDrawing();
        ClearBackground(renderSettings.backgroundColor);
        
        BeginMode3D(camera);
        
        // Set Global Uniforms
        
        Vector3 sunColorValue = { renderSettings.sunColor.r / 255.0f, renderSettings.sunColor.g / 255.0f, renderSettings.sunColor.b / 255.0f };
        Vector3 skyColorValue = { renderSettings.skyColor.r / 255.0f, renderSettings.skyColor.g / 255.0f, renderSettings.skyColor.b / 255.0f };
        float objectSpecularity = 0.5f;
        float objectGlossiness = 10.0f;
        float objectOpacity = 1.0f;
        
        Vector3 sunLightPosition = Vector3RotateByQuaternion((Vector3){ 0.0f, 0.0f, 1.0f }, QuaternionFromAxisAngle((Vector3){ 0.0f, 1.0f, 0.0f }, renderSettings.sunAzimuth));
        Vector3 sunLightAxis = Vector3Normalize(Vector3CrossProduct(sunLightPosition, (Vector3){ 0.0f, 1.0f, 0.0f }));
        Vector3 sunLightDir = Vector3Negate(Vector3RotateByQuaternion(sunLightPosition, QuaternionFromAxisAngle(sunLightAxis, renderSettings.sunAltitude)));
        
        SetShaderValue(shader, uniforms.cameraPosition, &camera.position, SHADER_UNIFORM_VEC3);
        SetShaderValue(shader, uniforms.exposure, &renderSettings.exposure, SHADER_UNIFORM_FLOAT);
        SetShaderValue(shader, uniforms.sunConeAngle, &renderSettings.sunLightConeAngle, SHADER_UNIFORM_FLOAT);
        SetShaderValue(shader, uniforms.sunDir, &sunLightDir, SHADER_UNIFORM_VEC3);
        SetShaderValue(shader, uniforms.sunStrength, &renderSettings.sunLightStrength, SHADER_UNIFORM_FLOAT);
        SetShaderValue(shader, uniforms.sunColor, &sunColorValue, SHADER_UNIFORM_VEC3);
        SetShaderValue(shader, uniforms.skyStrength, &renderSettings.skyLightStrength, SHADER_UNIFORM_FLOAT);
        SetShaderValue(shader, uniforms.skyColor, &skyColorValue, SHADER_UNIFORM_VEC3);
        SetShaderValue(shader, uniforms.ambientStrength, &renderSettings.ambientLightStrength, SHADER_UNIFORM_FLOAT);
        SetShaderValue(shader, uniforms.groundStrength, &renderSettings.groundLightStrength, SHADER_UNIFORM_FLOAT);
        SetShaderValue(shader, uniforms.objectSpecularity, &objectSpecularity, SHADER_UNIFORM_FLOAT);
        SetShaderValue(shader, uniforms.objectGlossiness, &objectGlossiness, SHADER_UNIFORM_FLOAT);
        SetShaderValue(shader, uniforms.objectOpacity, &objectOpacity, SHADER_UNIFORM_FLOAT);
        SetShaderValue(shader, uniforms.shadowMode, &renderSettings.shadowMode, SHADER_UNIFORM_INT);

        // Draw Ground        
        
        if (renderSettings.drawChecker)
        {
            int groundIsCapsule = 0;
            Vector3 groundColor = { 0.9f, 0.9f, 0.9f };
            
            SetShaderValue(shader, uniforms.isCapsule, &groundIsCapsule, SHADER_UNIFORM_INT);
            SetShaderValue(shader, uniforms.objectColor, &groundColor, SHADER_UNIFORM_VEC3);

            for (int i = 0; i < 11; i++)
            {
                for (int j = 0; j < 11; j++)
                {
                    Vector3 groundSegmentPosition = {
                        (((float)i / 10) - 0.5f) * 20.0f,
                        -0.001f,
                        (((float)j / 10) - 0.5f) * 20.0f, 
                    };
                    
                    capsuleData.aoCapsuleCount = 0;
                    capsuleData.shadowCapsuleCount = 0;
                    
                    if (renderSettings.drawCapsules && renderSettings.drawAO)
                    {
                        CapsuleDataUpdateAOCapsulesForGroundSegment(&capsuleData, groundSegmentPosition);
                    }
                    
                    if (renderSettings.drawCapsules && renderSettings.shadowMode != 0)
                    {
                        CapsuleDataUpdateShadowCapsulesForGroundSegment(&capsuleData, groundSegmentPosition, sunLightDir, renderSettings.sunLightConeAngle, renderSettings.shadowMode);
                    }
                    
                    int aoCapsuleCount = min(capsuleData.aoCapsuleCount, AO_CAPSULES_MAX);
                    int shadowCapsuleCount = min(capsuleData.shadowCapsuleCount, SHADOW_CAPSULES_MAX);
                    
                    SetShaderValue(shader, uniforms.shadowCapsuleCount, &shadowCapsuleCount, SHADER_UNIFORM_INT);
                    SetShaderValueV(shader, uniforms.shadowCapsuleStarts, capsuleData.shadowCapsuleStarts, SHADER_UNIFORM_VEC3, shadowCapsuleCount);
                    SetShaderValueV(shader, uniforms.shadowCapsuleVectors, capsuleData.shadowCapsuleVectors, SHADER_UNIFORM_VEC3, shadowCapsuleCount);        
                    SetShaderValueV(shader, uniforms.shadowCapsuleRadii, capsuleData.shadowCapsuleRadii, SHADER_UNIFORM_FLOAT, shadowCapsuleCount);
                    
                    SetShaderValue(shader, uniforms.aoCapsuleCount, &aoCapsuleCount, SHADER_UNIFORM_INT);
                    SetShaderValueV(shader, uniforms.aoCapsuleStarts, capsuleData.aoCapsuleStarts, SHADER_UNIFORM_VEC3, aoCapsuleCount);
                    SetShaderValueV(shader, uniforms.aoCapsuleVectors, capsuleData.aoCapsuleVectors, SHADER_UNIFORM_VEC3, aoCapsuleCount);        
                    SetShaderValueV(shader, uniforms.aoCapsuleRadii, capsuleData.aoCapsuleRadii, SHADER_UNIFORM_FLOAT, aoCapsuleCount);
                    
                    DrawModel(groundPlaneModel, groundSegmentPosition, 1.0f, WHITE);
                }
            }
        }
        
        // Draw Capsules
        
        
        if (renderSettings.drawCapsules)
        {
            int capsuleIsCapsule = 1;
            SetShaderValue(shader, uniforms.isCapsule, &capsuleIsCapsule, SHADER_UNIFORM_INT);
            SetShaderValue(shader, uniforms.objectOpacity, &renderSettings.capsuleOpacity, SHADER_UNIFORM_FLOAT);
            
            // Depth Sort
            
            for (int i = 0; i < capsuleData.capsuleCount; i++)
            {
                capsuleData.capsuleSort[i].index = i;
                capsuleData.capsuleSort[i].value = Vector3Distance(camera.position, capsuleData.capsulePositions[i]);
            }
            
            if (renderSettings.capsuleOpacity < 1.0f)
            {
                qsort(capsuleData.capsuleSort, capsuleData.capsuleCount, sizeof(CapsuleSort), CapsuleSortCompareLess);
                rlDrawRenderBatchActive();        
                rlDisableDepthTest(); 
            }
            else
            {
                qsort(capsuleData.capsuleSort, capsuleData.capsuleCount, sizeof(CapsuleSort), CapsuleSortCompareGreater);
            }
            
            // Render
            
            for (int i = 0; i < capsuleData.capsuleCount; i++)
            {
                int j = capsuleData.capsuleSort[i].index;
              
                Vector3 capsuleStart = CapsuleStart(capsuleData.capsulePositions[j], capsuleData.capsuleRotations[j], capsuleData.capsuleHalfLengths[j]);
                Vector3 capsuleVector = CapsuleVector(capsuleData.capsulePositions[j], capsuleData.capsuleRotations[j], capsuleData.capsuleHalfLengths[j]);

                capsuleData.aoCapsuleCount = 0;
                capsuleData.shadowCapsuleCount = 0;

                if (renderSettings.drawAO)
                {
                    CapsuleDataUpdateAOCapsulesForCapsule(&capsuleData, j);                  
                }
                
                if (renderSettings.shadowMode != 0)
                {
                    CapsuleDataUpdateShadowCapsulesForCapsule(&capsuleData, j, sunLightDir, renderSettings.sunLightConeAngle, renderSettings.shadowMode);
                }
                
                int shadowCapsuleCount = min(capsuleData.shadowCapsuleCount, SHADOW_CAPSULES_MAX);
                int aoCapsuleCount = min(capsuleData.aoCapsuleCount, AO_CAPSULES_MAX);                
                
                SetShaderValue(shader, uniforms.objectColor, &capsuleData.capsuleColors[j], SHADER_UNIFORM_VEC3);

                SetShaderValue(shader, uniforms.shadowCapsuleCount, &shadowCapsuleCount, SHADER_UNIFORM_INT);
                SetShaderValueV(shader, uniforms.shadowCapsuleStarts, capsuleData.shadowCapsuleStarts, SHADER_UNIFORM_VEC3, shadowCapsuleCount);
                SetShaderValueV(shader, uniforms.shadowCapsuleVectors, capsuleData.shadowCapsuleVectors, SHADER_UNIFORM_VEC3, shadowCapsuleCount);        
                SetShaderValueV(shader, uniforms.shadowCapsuleRadii, capsuleData.shadowCapsuleRadii, SHADER_UNIFORM_FLOAT, shadowCapsuleCount);
                
                SetShaderValue(shader, uniforms.aoCapsuleCount, &aoCapsuleCount, SHADER_UNIFORM_INT);
                SetShaderValueV(shader, uniforms.aoCapsuleStarts, capsuleData.aoCapsuleStarts, SHADER_UNIFORM_VEC3, aoCapsuleCount);
                SetShaderValueV(shader, uniforms.aoCapsuleVectors, capsuleData.aoCapsuleVectors, SHADER_UNIFORM_VEC3, aoCapsuleCount);        
                SetShaderValueV(shader, uniforms.aoCapsuleRadii, capsuleData.aoCapsuleRadii, SHADER_UNIFORM_FLOAT, aoCapsuleCount);
                
                SetShaderValue(shader, uniforms.capsulePosition, &capsuleData.capsulePositions[j], SHADER_UNIFORM_VEC3);
                SetShaderValue(shader, uniforms.capsuleRotation, &capsuleData.capsuleRotations[j], SHADER_UNIFORM_VEC4);
                SetShaderValue(shader, uniforms.capsuleHalfLength, &capsuleData.capsuleHalfLengths[j], SHADER_UNIFORM_FLOAT);
                SetShaderValue(shader, uniforms.capsuleRadius, &capsuleData.capsuleRadii[j], SHADER_UNIFORM_FLOAT);
                SetShaderValue(shader, uniforms.capsuleStart, &capsuleStart, SHADER_UNIFORM_VEC3);
                SetShaderValue(shader, uniforms.capsuleVector, &capsuleVector, SHADER_UNIFORM_VEC3);
                
                DrawModel(capsuleModel, Vector3Zero(), 1.0f, WHITE);
            }
        }
        
        if (renderSettings.capsuleOpacity < 1.0f)
        {
            rlDrawRenderBatchActive();
            rlEnableDepthTest();
        }
        
        // Draw Debug
        
        if (renderSettings.drawGrid)
        {
            DrawGrid(20, 1.0f);
        }
        
        if (renderSettings.drawOrigin)
        {
            DrawTransform(
                (Vector3){ 0.0f, 0.01f, 0.0f }, 
                QuaternionIdentity(),
                1.0f);
        }  
        
        rlDrawRenderBatchActive();        
        rlDisableDepthTest(); 
        
        // Draw Capsule Wireframes
        
        if (renderSettings.drawWireframes)
        {
            DrawWireFrames(&capsuleData, DARKGRAY);
        }
        
        // Draw Bones
        
        if (renderSettings.drawSkeleton)
        {
            for (int i = 0; i < bvhCount; i++)
            {
                DrawSkeleton(
                    &xformData[i],
                    renderSettings.drawEndSites,
                    DARKGRAY,
                    GRAY);
            }
        }
        
        // Draw Joint Transforms
        
        if (renderSettings.drawTransforms)
        {
            for (int i = 0; i < bvhCount; i++)
            {
                DrawTransforms(&xformData[i]);
            }
        }
        
        rlDrawRenderBatchActive();
        rlEnableDepthTest();
        
        EndMode3D();
        
        // Draw UI
        
        if (IsKeyPressed(KEY_H) && !fileDialogState.windowActive)
        {
            renderSettings.drawUI = !renderSettings.drawUI;
        }
        
        if (renderSettings.drawUI)
        {
            if (fileDialogState.windowActive) { GuiLock(); }

            // Error Message
            
            DrawText(errMsg, 250, 20, 15, RED);

            // Rendering            
            
            GuiRenderSettings(&renderSettings, screenWidth, screenHeight);
        
            // FPS
            
            if (renderSettings.drawFPS)
            {
                DrawFPS(230, 10);              
            }
            
            // Camera
            
            GuiGroupBox((Rectangle){ 20, 10, 190, 200 }, "Camera");
            
            GuiLabel((Rectangle){ 30, 20, 150, 20 }, "Ctrl + Left Click - Rotate");
            GuiLabel((Rectangle){ 30, 40, 150, 20 }, "Mouse Scroll - Zoom");
            GuiLabel((Rectangle){ 30, 60, 150, 20 }, TextFormat("Target: [% 5.3f % 5.3f % 5.3f]", cameraTarget.x, cameraTarget.y, cameraTarget.z));
            GuiLabel((Rectangle){ 30, 80, 150, 20 }, TextFormat("Azimuth: %5.3f", cameraAzimuth));
            GuiLabel((Rectangle){ 30, 100, 150, 20 }, TextFormat("Altitude: %5.3f", cameraAltitude));
            GuiLabel((Rectangle){ 30, 120, 150, 20 }, TextFormat("Distance: %5.3f", cameraDistance));
            
            if (bvhCount > 0)
            {
                GuiToggle((Rectangle){ 30, 150, 100, 20 }, "Track", &cameraTrack);
                GuiComboBox((Rectangle){ 30, 180, 150, 20 }, bvhData[bvhActive].jointNamesCombo, &cameraTrackBone);          
            }
            
            // Files
            
            GuiGroupBox((Rectangle){ 20, 230, 190, 40 + bvhCount * 60 }, "Files");
            
            if (GuiButton((Rectangle){ 30, 240, 110, 20 }, "Open"))
            {
                fileDialogState.windowActive = true;
            }
            
            if (GuiButton((Rectangle){ 150, 240, 50, 20 }, "Clear"))
            {
                bvhCount = 0;
                errMsg[0] = '\0';
            }
            
            for (int i = 0; i < bvhCount; i++)
            {
                char bvhNameShort[20];
                bvhNameShort[0] = '\0';
                if (strlen(bvhNames[i]) + 1 <= 20)
                {
                    strcat(bvhNameShort, bvhNames[i]);
                }
                else
                {
                    memcpy(bvhNameShort, bvhNames[i], 16);
                    memcpy(bvhNameShort + 16, "...", 4);
                }
                
                bool bvhSelected = i == bvhActive;
                GuiToggle((Rectangle){ 30, 270 + i * 60, 140, 20 }, bvhNameShort, &bvhSelected);
                
                if (bvhSelected)
                {
                    bvhActive = i;
                    cameraTrackBone = min(cameraTrackBone, xformData[bvhActive].jointCount - 1);
                }
                
                DrawRectangleRec((Rectangle){ 180, 270 + i * 60, 20, 20 }, bvhColors[i]);
                DrawRectangleLinesEx((Rectangle){ 180, 270 + i * 60, 20, 20 }, 1, GRAY);
                
                bool scaleM = bvhScales[i] == 1.0f;
                GuiToggle((Rectangle){ 30, 300 + i * 60, 30, 20 }, "m", &scaleM); if (scaleM) { bvhScales[i] = 1.0f; }
                bool scaleCM = bvhScales[i] == 0.01f;
                GuiToggle((Rectangle){ 65, 300 + i * 60, 30, 20 }, "cm", &scaleCM); if (scaleCM) { bvhScales[i] = 0.01f; }
                bool scaleInches = bvhScales[i] == 0.0254f;
                GuiToggle((Rectangle){ 100, 300 + i * 60, 30, 20 }, "inch", &scaleInches); if (scaleInches) { bvhScales[i] = 0.0254f; }
                bool scaleFeet = bvhScales[i] == 0.3048f;
                GuiToggle((Rectangle){ 135, 300 + i * 60, 30, 20 }, "feet", &scaleFeet); if (scaleFeet) { bvhScales[i] = 0.3048f; }
                bool scaleAuto = bvhScales[i] == bvhAutoScales[i];
                GuiToggle((Rectangle){ 170, 300 + i * 60, 30, 20 }, "auto", &scaleAuto); if (scaleAuto) { bvhScales[i] = bvhAutoScales[i]; }
            }
            
            // Scrubber
            
            if (bvhCount > 0)
            {
                GuiLabel((Rectangle){ 160, screenHeight - 80, 150, 20 }, TextFormat("Frame Time: %f", bvhData[bvhActive].frameTime));
                GuiCheckBox((Rectangle){ 290, screenHeight - 80, 20, 20 }, "Snap to Frame", &frameSnap);
                GuiComboBox((Rectangle){ 400, screenHeight - 80, 100, 20 }, "Nearest;Linear;Cubic", &sampleMode);
            
                GuiToggle((Rectangle){ screenWidth / 2 - 25, screenHeight - 80, 50, 20 }, "Play", &playing);
                GuiToggle((Rectangle){ screenWidth / 2 - 85, screenHeight - 80, 50, 20 }, "Loop", &looping);
                
                bool speed01x = playSpeed == 0.1f;
                GuiToggle((Rectangle){ screenWidth / 2 + 40, screenHeight - 80, 30, 20 }, "0.1x", &speed01x); if (speed01x) { playSpeed = 0.1f; }
                bool speed05x = playSpeed == 0.5f;
                GuiToggle((Rectangle){ screenWidth / 2 + 80, screenHeight - 80, 30, 20 }, "0.5x", &speed05x); if (speed05x) { playSpeed = 0.5f; }
                bool speed1x = playSpeed == 1.0f;
                GuiToggle((Rectangle){ screenWidth / 2 + 120, screenHeight - 80, 30, 20 }, "1x", &speed1x); if (speed1x) { playSpeed = 1.0f; }
                bool speed2x = playSpeed == 2.0f;
                GuiToggle((Rectangle){ screenWidth / 2 + 160, screenHeight - 80, 30, 20 }, "2x", &speed2x); if (speed2x) { playSpeed = 2.0f; }
                bool speed4x = playSpeed == 4.0f;
                GuiToggle((Rectangle){ screenWidth / 2 + 200, screenHeight - 80, 30, 20 }, "4x", &speed4x); if (speed4x) { playSpeed = 4.0f; }
                GuiSliderBar((Rectangle){ screenWidth / 2 + 240, screenHeight - 80, 70, 20 }, "", TextFormat("%5.2fx", playSpeed), &playSpeed, 0.0f, 4.0f);
              
                frame = clamp((int)(playTime / bvhData[bvhActive].frameTime + 0.5f), frameMin, frameMax);

                if (GuiValueBox(
                    (Rectangle){ 100, screenHeight - 80, 50, 20 }, 
                    "Min   ", &frameMinSelect, 0, frameLimit, frameMinEdit))
                {
                    frameMinEdit = !frameMinEdit;
                    if (!frameMinEdit)
                    {
                        frameMin = frameMinSelect;
                        frameMaxSelect = frameMaxSelect < frameMin ? frameMin : frameMaxSelect;
                        frameMax = frameMaxSelect;
                        frame = frame < frameMin ? frameMin : frame;
                        playTime = frame * bvhData[bvhActive].frameTime;
                        timeMin = frameMin * bvhData[bvhActive].frameTime;
                        timeMax = frameMax * bvhData[bvhActive].frameTime;
                    }
                }
                
                if (GuiValueBox(
                    (Rectangle){ screenWidth - 170, screenHeight - 80, 50, 20 }, 
                    "Max   ", &frameMaxSelect, 0, frameLimit, frameMaxEdit))
                {
                    frameMaxEdit = !frameMaxEdit;
                    
                    if (!frameMaxEdit)
                    {
                        frameMax = frameMaxSelect;
                        frameMinSelect = frameMinSelect > frameMax ? frameMax : frameMinSelect;
                        frameMin = frameMinSelect;
                        frame = frame > frameMax ? frameMax : frame;
                        playTime = frame * bvhData[bvhActive].frameTime;
                        timeMin = frameMin * bvhData[bvhActive].frameTime;
                        timeMax = frameMax * bvhData[bvhActive].frameTime;
                    }
                }
                
                GuiLabel(
                    (Rectangle){ screenWidth - 110, screenHeight - 80, 100, 20 }, 
                    TextFormat("of %i", frameLimit));
                
                float frameFloatPrev = frameSnap ? (float)frame : playTime / bvhData[bvhActive].frameTime;
                float frameFloat = frameFloatPrev;
                
                GuiSliderBar(
                    (Rectangle){ 100, screenHeight - 50, screenWidth - 220, 20 }, 
                    TextFormat("%5.2f", playTime), 
                    TextFormat("%i", frame),
                    &frameFloat, 
                    (float)frameMin, (float)frameMax);
                
                if (frameFloat != frameFloatPrev)
                {
                    if (frameSnap)
                    {
                        frame = clamp((int)(frameFloat + 0.5f), frameMin, frameMax);
                        playTime = frame * bvhData[bvhActive].frameTime;
                    }
                    else
                    {
                        playTime = frameFloat * bvhData[bvhActive].frameTime;
                    }
                }
            }
        
            GuiUnlock();
            
            GuiWindowFileDialog(&fileDialogState);
        }
        
        // Done

        EndDrawing();
    }
  
    // Unload stuff and finish
  
    CapsuleDataFree(&capsuleData);
  
    for (int i = 0; i < bvhCount; i++)
    {
        TransformDataFree(&xformData[i]);
        TransformDataFree(&xformTmp0[i]);
        TransformDataFree(&xformTmp1[i]);
        TransformDataFree(&xformTmp2[i]);
        TransformDataFree(&xformTmp3[i]);
        BVHDataFree(&bvhData[i]);
    }
    
    UnloadModel(capsuleModel);    
    UnloadModel(groundPlaneModel);
    UnloadShader(shader);

    CloseWindow();

    return 0;
}