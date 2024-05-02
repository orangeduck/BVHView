#include "additions.h"
#include "raylib.h"

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
