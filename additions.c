#include "additions.h"
#include "raylib.h"
#include "raudio.c"

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
        fprintf(stderr, "FFmpegPipe is null\n");
        return false;
    }

    if (!pipe->outputPath)
    {
        fprintf(stderr, "FFmpegPipe output path is null/missing\n");
        return false;
    }

    char ffmpegCommand[FFMPEG_COMMAND_BUFFER_SIZE];
    snprintf(ffmpegCommand, FFMPEG_COMMAND_BUFFER_SIZE,
        "ffmpeg -y -f rawvideo -pixel_format rgb24 -video_size %dx%d -framerate %d -i - -vf format=yuv420p %s",
        pipe->width,
        pipe->height,
        pipe->framerate,
        pipe->outputPath
    );

    pipe->pipeHandle = _popen(ffmpegCommand, "w");

    if (!pipe->pipeHandle) {
        fprintf(stderr, "Failed to open pipe to FFmpeg\n");
        return false;
    }

    return true;
}

void WriteImageToFFmpegPipe(FFmpegPipe* pipe, Image* image)
{
    if (pipe && pipe->pipeHandle && image)
    {
        if (pipe->width == image->width && pipe->height == image->height)
        {
            fwrite(image->data, 1, image->width * image->height * 3, pipe->pipeHandle);
        }
    }
}

void CloseFFmpegPipe(FFmpegPipe* pipe)
{
    if (pipe && pipe->pipeHandle)
    {
        _pclose(pipe->pipeHandle);
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
}
