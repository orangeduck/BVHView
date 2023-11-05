# BVHView

BVHView is a simple .bvh animation file format viewer built using [raylib](https://www.raylib.com/).

TODO: Video

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
