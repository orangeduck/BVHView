#include "additions.h"
#include "raylib.h"
#include "raudio.c"
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
        fprintf(stderr, "Error: FFmpegPipe is null\n");
        return false;
    }

    char ffmpegCommand[FFMPEG_COMMAND_BUFFER_SIZE];
    snprintf(ffmpegCommand, FFMPEG_COMMAND_BUFFER_SIZE,
        "ffmpeg -y -f rawvideo -pixel_format rgba -video_size %dx%d -framerate %d -i - -vf format=yuv420p \"%s\"",
        pipe->width,
        pipe->height,
        pipe->framerate,
        pipe->outputPath
    );

#if defined(_WIN32)
    pipe->pipeHandle = _popen(ffmpegCommand, "wb");
#else
    pipe->pipeHandle = popen(ffmpegCommand, "w");
#endif

    if (!pipe->pipeHandle) {
        fprintf(stderr, "Error: Failed to open pipe to FFmpeg: %s\n", strerror(errno));
        return false;
    }

    return true;
}

bool WriteImageToFFmpegPipe(FFmpegPipe* pipe, Image* image)
{
    if (!pipe) {
        fprintf(stderr, "Error: FFmpeg pipe is null.\n");
        return false;
    }

    if (!pipe->pipeHandle) {
        fprintf(stderr, "Error: FFMpeg pipe handle is null.\n");
        return false;
    }

    if (!image) {
        fprintf(stderr, "Error: Image is null.\n");
        return false;
    }

    if (pipe->width != image->width || pipe->height != image->height) {
        fprintf(stderr, "Error: Image dimensions do not match pipe dimensions (%dx%d vs %dx%d).\n",
                pipe->width, pipe->height, image->width, image->height);
        return false;
    }

    if (image->format != PIXELFORMAT_UNCOMPRESSED_R8G8B8A8) {
        fprintf(stderr, "Error: Image data format must be %d but got %d.\n",
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
