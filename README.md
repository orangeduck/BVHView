# BVHView

BVHView is a simple .bvh animation file format viewer built using [raylib](https://www.raylib.com/).

* [Download (Windows)](https://theorangeduck.com/media/uploads/BVHView/bvhview.zip)
* [Project Page](https://theorangeduck.com/page/bvhview)
* [Web Demo](https://theorangeduck.com/media/uploads/BVHView/bvhview.html)

TODO: Video

# Building

## Windows

Download [raylib](https://github.com/raysan5/raylib) into `C:/raylib/raylib`.

Download [raygui](https://github.com/raysan5/raygui) into `C:/raylib/raygui`.

Build raylib by going to `C:/raylib/raylib/src` and running `make`.

Download `BVHView` and run `make` in the main directory to build `bvhview.exe`.

To build a release version with optimizations enabled and no console window alongside run `make BUILD_MODE=RELEASE` in the main directory instead.

## Linux

Download [raylib](https://github.com/raysan5/raylib) into `~/raylib/raylib`.

Download [raygui](https://github.com/raysan5/raygui) into `~/raylib/raygui`.

Build raylib by going to `~/raylib/raylib/src` and running `make`.

Download `BVHView` and run `make` in the main directory to build `bvhview`.

To build a release version with optimizations enabled run `make BUILD_MODE=RELEASE` in the main directory instead.

