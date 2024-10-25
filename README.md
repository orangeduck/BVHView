# BVHView (fork)

This repo is a fork of [BVHView](https://github.com/orangeduck/BVHView) developed by [Daniel Holden](https://theorangeduck.com/). This BVHView version has been developed for the [GENEA Leaderboard](https://genea-workshop.github.io/leaderboard/) research project.

## Roadmap

- [x] Add textured mesh for GENEA avatar
- [X] Playback BVH animation onto GENEA avatar
- [X] Add scrubber-synchronized WAV audio
- [X] Load BVH and WAV via command line args
- [ ] Add SMPL-X meshes for the BEAT dataset
- [X] Add .mp4 video recordings using FFMPEG (for automated pipelines)
- [ ] Add audio recording to .mp4 videos
- [ ] Add support for multiple animated models in the same scene

(The list could change with time.)

## Dependencies

### Raylib

Used to render window, GUI, and OpenGL graphics.

1. Prepare `raylib` (refer to the [original BVHView instructions](https://github.com/orangeduck/BVHView)).
2. Uncomment `#define SUPPORT_FILEFORMAT_JPG` in `raylib/raylib/src/config.h`
   - This allows loading `.jpg` textures for user meshes
   - Other formats have not been tested, but feel free to try
3. Build `raylib` by running `make -B` in `raylib/raylib/src`
4. Build `raygui` (refer to the [original BVHView instructions](https://github.com/orangeduck/BVHView))

### cwalk

Convenience library for manipulating system paths.

1. Download `include/cwalk.h` and `src/cwalk.c` from [here](https://github.com/likle/cwalk).
2. Copy `cwalk.h` and `cwalk.c` to the `external` folder in this repository.
3. The `Makefile` will handle the rest when you start building.

### FFmpeg

Used to stream framebuffer data from raylib for recording videos.

1. Refer to the [original FFmpeg instructions](https://www.ffmpeg.org/download.html).

## Building

Refer to the [original BVHView instructions](https://github.com/orangeduck/BVHView).

## Usage

It is recommended that you run the software via the *command line interface (CLI)*, as some features (such as `.wav` loading) are only supported via command line arguments.

CLI example:
- `cd "[...]/BVHView/"`
- `./bvhview.exe --bvh="./samples/genea/trn_2023_v0_000_main-agent.bvh" --wav="./samples/genea/trn_2023_v0_000_main-agent.wav"` *(Remove `.exe` on Linux)*

You can find example files in `/samples`.

### Args

Most args are optional.

- `--bvh` : Path to a `.bvh` animation data file
- `--wav` : Path to a `.wav` audio file
- `--record` : Toggle recording of the 3D scene (this hides the window and UI). *Default is disabled.*
- `--record_fps` : Toggle the FPS at which to record. *Default is 30.*
- `--record_out_dir` : Specify an absolute or relative path where recordings will be saved. *Default is `output/video` in the working directory.*
- `--record_out_name` : Specify the name of the recording. *Default is the name of the `.bvh` file.*

### Executable

Windows:
- For releases at `BVHView/bvhview.exe`
- For your own builds at `build/BVHView/bvhview.exe`

Linux:
- There are no releases currently
- For your own builds at `build/BVHView/bvhview`

## Importing your mesh

Loading your mesh in BVHView is trickier than other DCC and game engines which have been developed for a long time. For example, there is currently no `.fbx` support to import a mesh (but we can use `.gltf ` instead).

The steps below guide you through the process of preparing your mesh for BVHView. Note that these steps may differ depending on your character -- you may need to play with the import/export settings, textures, bone orientation.

1. Setup
   1. Install `Blender` ([Official website](https://www.blender.org/download/), [Steam](https://store.steampowered.com/app/365670/Blender/))
   2. Start `Blender`
   3. Delete all scene objects (Camera, Cube, Light...)
2. Import
   1. `File > Import > FBX (.fbx)`
   2. Select your `.fbx` mesh
   3. Import settings:
      1. Reset settings via `Restore Operator Defaults`
      2. Disable `Animation`
   4. Press `Import FBX`
3. Configure model
   1. Select your armature and then `Object > Apply > All transforms`
   2. Select your mesh and then `Object > Apply > All transforms`
4. Export
   1. `File > Export > glTF 2.0 (.glb/.gltf)`
   2. Choose export directory
      - *Example: `[...]/BVHView/assets/models/<your_model_folder/>`*
   3. Choose a file name
      - *Warning: you should update the contents of the exported `.gltf` file if you rename any files after the export*
   4. Export settings:
      1. Reset settings via `Restore Operator Defaults`
      2. `Format > glTF Separate (.gltf + .bin + textures)`
      3. Disable `Animation`
      4. *(You may need to set `Textures` to export the textures as files.)*
   5. Press `Export glTF 2.0`
   6. Test export by running `./bvhview --mesh="<path_to_gltf>"`

### Texture issues

If your textures are not loaded correctly or visible, it could be due to:
- Texture file is not `.png` or `.jpg`
  - Only `.png` and `.jpg` are supported currently. For other formats, you should uncomment the corresponding lines in `/raylib/raylib/src/config.h`, rebuild `raylib`, then rebuild `BVHView`.

---
---

# BVHView (original description)

BVHView is a simple [.bvh animation file format](https://research.cs.wisc.edu/graphics/Courses/cs-838-1999/Jeff/BVH.html) viewer built using [raylib](https://www.raylib.com/).



https://github.com/orangeduck/BVHView/assets/177299/9f976284-02c1-4a14-bca4-b8b3d7c9a774



* [Download (Windows)](https://theorangeduck.com/media/uploads/BVHView/bvhview.zip)
* [Project Page](https://theorangeduck.com/page/bvhview)
* [Web Demo](https://theorangeduck.com/media/uploads/BVHView/bvhview.html)

# Building

## Windows

Download and install [MinGW](https://www.mingw-w64.org/) in some form. Perhaps [w64devkit](https://www.mingw-w64.org/downloads/#w64devkit) or [MSYS2](https://www.mingw-w64.org/downloads/#msys2).

Download [raylib](https://github.com/raysan5/raylib) into `C:/raylib/raylib`.

Download [raygui](https://github.com/raysan5/raygui) into `C:/raylib/raygui`.

Build raylib by going to `C:/raylib/raylib/src` and running `make`.

Download this repo and run `make` in the main directory to build `bvhview.exe`.

To build a release version with optimizations enabled and no console window run `make BUILD_MODE=RELEASE` in the main directory instead.

## Linux

Download [raylib](https://github.com/raysan5/raylib) into `~/raylib/raylib`.

Download [raygui](https://github.com/raysan5/raygui) into `~/raylib/raygui`.

Follow [this guide](https://github.com/raysan5/raylib/wiki/Working-on-GNU-Linux) to install any dependencies. 

Build raylib by going to `~/raylib/raylib/src` and running `make`.

Download this repo and run `make` in the main directory to build `bvhview`.

To build a release version with optimizations enabled run `make BUILD_MODE=RELEASE` in the main directory instead.

## Other

For other platforms you should be able to build BVHView by hacking the `Makefile` a bit. Contributions here welcome.
