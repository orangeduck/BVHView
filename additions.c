#include "additions.h"
#include "raylib.h"
#include "raudio.c"
#include <ctype.h>
#include <errno.h>
#include <stdio.h>

bool BVHALoadCharacterModelFromFile(CharacterModel* characterModel, const char *fileName)
{
    Model loadedModel = LoadModel(fileName);
    characterModel->model = loadedModel;
    characterModel->isLoaded = true;
    return true;
}

void BVHAUnloadCharacterModel(CharacterModel* characterModel)
{
    UnloadModel(characterModel->model);
    characterModel->isLoaded = false;
}

bool boneExistsInModel(const char* boneName, const Model* model) {
    for (int i = 0; i < model->boneCount; ++i) {
        if (strcmp(boneName, model->bones[i].name) == 0) {
            return true;
        }
    }
    return false;
}

void SetAudioTimeInSeconds(Sound* audio, float seconds)
{
    audio->stream.buffer->frameCursorPos = audio->stream.sampleRate * seconds;
}

bool OpenFFmpegPipe(FFmpegPipe* pipe)
{
    if (!pipe)
    {
        fprintf(stderr, "[ERROR - BVHVIEW] FFmpegPipe is null\n");
        return false;
    }

    char ffmpegCommand[FFMPEG_COMMAND_BUFFER_SIZE];
    if (pipe->audioPath[0] == '\0') {
        snprintf(ffmpegCommand, FFMPEG_COMMAND_BUFFER_SIZE,
            "ffmpeg -y -f rawvideo -pixel_format rgba -video_size %dx%d -framerate %d -i - -vf format=yuv420p \"%s\"",
            pipe->width,
            pipe->height,
            pipe->framerate,
            pipe->outputPath
        );
    } else {
        snprintf(ffmpegCommand, FFMPEG_COMMAND_BUFFER_SIZE,
            "ffmpeg -y -f rawvideo -pixel_format rgba -video_size %dx%d -framerate %d -i - -i \"%s\" -c:v libx264 -vf format=yuv420p -c:a aac -shortest \"%s\"",
            pipe->width,
            pipe->height,
            pipe->framerate,
            pipe->audioPath,
            pipe->outputPath
        );
    }

#if defined(_WIN32)
    pipe->pipeHandle = _popen(ffmpegCommand, "wb");
#else
    pipe->pipeHandle = popen(ffmpegCommand, "w");
#endif

    if (!pipe->pipeHandle) {
        fprintf(stderr, "[ERROR - BVHVIEW] Failed to open pipe to FFmpeg: %s\n", strerror(errno));
        return false;
    }

    return true;
}

bool WriteImageToFFmpegPipe(FFmpegPipe* pipe, Image* image)
{
    if (!pipe) {
        fprintf(stderr, "[ERROR - BVHVIEW] FFmpeg pipe is null.\n");
        return false;
    }

    if (!pipe->pipeHandle) {
        fprintf(stderr, "[ERROR - BVHVIEW] FFMpeg pipe handle is null.\n");
        return false;
    }

    if (!image) {
        fprintf(stderr, "[ERROR - BVHVIEW] Image is null.\n");
        return false;
    }

    if (pipe->width != image->width || pipe->height != image->height) {
        fprintf(stderr, "[ERROR - BVHVIEW] Image dimensions do not match pipe dimensions (%dx%d vs %dx%d).\n",
                pipe->width, pipe->height, image->width, image->height);
        return false;
    }

    if (image->format != PIXELFORMAT_UNCOMPRESSED_R8G8B8A8) {
        fprintf(stderr, "[ERROR - BVHVIEW] Image data format must be %d but got %d.\n",
                PIXELFORMAT_UNCOMPRESSED_R8G8B8A8, image->format);
        return false;
    }

    const int pixel_size = 4; // bytes
    const int pixel_count = image->width * image->height;
    fwrite(image->data, pixel_size, pixel_count, pipe->pipeHandle);
    return true;
}

void CloseFFmpegPipe(FFmpegPipe* pipe)
{
    if (pipe && pipe->pipeHandle)
    {
        pclose(pipe->pipeHandle);
        pipe->pipeHandle = NULL;
    }
}

void FreeFFmpegPipe(FFmpegPipe* pipe) {
    if (pipe == NULL) {
        return;
    }

    if (pipe->pipeHandle != NULL) {
        CloseFFmpegPipe(pipe);
    }

    pipe = NULL;
}

void NormalizePath(char *path) {
    char target_separator = PATH_SEPARATOR;
    char other_separator = (target_separator == '/') ? '\\' : '/';
    for (size_t i = 0; i < strlen(path); i++) {
        if (path[i] == other_separator) {
            path[i] = target_separator;
        }
    }
}

void TrimTrailingSpaces(char *str) {
    size_t len = strlen(str);
    while (len > 0 && isspace((unsigned char)str[len - 1])) {
        len--;
    }
    str[len] = '\0';
}

int CreateDirectories(const char *path) {
    char temp[PATH_MAX];

    // Copy path to temp and ensure null termination
    snprintf(temp, sizeof(temp), "%s", path);
    TrimTrailingSpaces(temp);
    size_t len = strlen(temp);
    size_t start = 0;

    // Normalize the path to use the platform-specific separator
    NormalizePath(temp);

#ifdef _WIN32
    // Handle Windows volume (drive letter) if present, e.g., "C:/"
    if (len > 2 && temp[1] == ':' && (temp[2] == PATH_SEPARATOR || temp[2] == '\0')) {
        start = 3; // Skip the "C:/" part of the path
    }
#endif

    // Iterate over each part of the path and create directories
    for (size_t i = start; i < len; i++) {
        if (temp[i] == PATH_SEPARATOR || temp[i] == '\0') {
            if (i == 0) {
                continue;
            }

            // Temporarily end the string to create the directory
            temp[i] = '\0';

            // Check if the directory exists
            struct stat st = {0};
            if (stat(temp, &st) == -1) {
                mkdir_p(temp);
                printf("[INFO - BVHVIEW] Created directory: %s\n", temp);
            }

            // Restore the path separator after creating the directory
            temp[i] = PATH_SEPARATOR;
        }
    }

    // Create the final directory if it's not handled in the loop
    struct stat st = {0};
    if (stat(temp, &st) == -1) {
        mkdir_p(temp);
        printf("[INFO - BVHVIEW] Created directory: %s\n", temp);
    }

    return 0;
}
