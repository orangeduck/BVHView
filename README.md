# BVHView (fork)

This repo is a fork of [BVHView](https://github.com/orangeduck/BVHView) developed by [Daniel Holden](https://theorangeduck.com/). This BVHView version has been developed for the [GENEA Leaderboard](https://genea-workshop.github.io/leaderboard/) research project.

## Roadmap

- [x] Add textured mesh for GENEA avatar
- [X] Playback BVH animation onto GENEA avatar
- [X] Add scrubber-synchronized WAV audio
- [X] Load BVH and WAV via command line args
- [ ] Add .mp4 recordings using FFMPEG (for automated pipelines)
- [ ] Add SMPL-X meshes for the BEAT dataset

(The list could change with time.)

## Building

Please follow the original BVHView build instructions. I work on Windows so Linux has not been tested. If there are build issues, or you have recommendations, please make a GitHub issue or contact me by email: tnikolov@hotmail.com.

## Using

### Windows

Locate the executable in `BVHView/bvhview.exe` for Releases, and `build/BVHView/bvhview.exe` when building yourself. It is recommended that you run the software via a **command line interface**, as some features (such as `.wav` loading) are only supported via command line arguments.

You can find example files in `BVHView/samples/genea/`.

Args:
- `--bvh` : Path to GENEA-compatible `.bvh` file
- `--wav` : Path to `.wav` file

CLI example:
- `cd ".../BVHView/"`
- `./bvhview.exe --bvh="./samples/genea/trn_2023_v0_000_main-agent.bvh" --wav="./samples/genea/trn_2023_v0_000_main-agent.wav"`

### Linux

Not tested.

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
