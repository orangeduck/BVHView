# BVHView

BVHView 是一个基于 [raylib](https://www.raylib.com/) 的 `.bvh` 动画文件查看器。

[中文](README_ZH.md)｜[English](README.md)

![BVHView 演示](https://github.com/orangeduck/BVHView/assets/177299/9f976284-02c1-4a14-bca4-b8b3d7c9a774)

## 获取与运行

已发布版本请前往 [GitHub Releases](https://github.com/cnwxi/BVHView/releases) 下载与系统和 CPU 架构一致的 ZIP 文件，解压后运行：

```bash
# macOS 或 Linux
./bvhview path/to/animation.bvh

# Windows
bvhview.exe path\to\animation.bvh
```

也可以在程序中点击 `Open` 打开文件，或将 `.bvh` 文件直接拖入窗口。

## 基本操作

- `Ctrl` + 鼠标左键：旋转相机。
- `Ctrl` + 鼠标右键：平移相机。
- 鼠标滚轮：缩放。
- `H`：显示或隐藏界面。

## 命令行参数

参数格式为 `--名称=值`。位置参数会作为要加载的 `.bvh` 文件。

```bash
./bvhview walk.bvh --screenWidth=1920 --screenHeight=1080 \
  --playing=true --looping=true --playSpeed=1.25 --sampleMode=linear
```

常用参数如下：

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `screenWidth` | `1280` | 窗口宽度。 |
| `screenHeight` | `720` | 窗口高度。 |
| `playing` | `true` | 是否自动播放，取值为 `true` 或 `false`。 |
| `looping` | `false` | 是否循环播放。 |
| `playSpeed` | `1.0` | 播放速度倍率。 |
| `playTime` | `0.0` | 初始播放时间。 |
| `sampleMode` | `linear` | 采样方式：`nearest`、`linear` 或 `cubic`。 |
| `cameraFOV` | `45.0` | 相机视场角。 |
| `cameraAzimuth` | `0.0` | 相机水平角度。 |
| `cameraAltitude` | `0.4` | 相机垂直角度。 |
| `cameraDistance` | `4.0` | 相机距离。 |

## 从源码构建

项目需要 [raylib](https://github.com/raysan5/raylib) 与 [raygui](https://github.com/raysan5/raygui)。请先将它们克隆到各平台约定的位置，再构建 raylib。

### Windows

安装 [MinGW-w64](https://www.mingw-w64.org/) 或 [MSYS2](https://www.msys2.org/)。在 MSYS2 的 MinGW shell 中执行：

```bash
git clone https://github.com/raysan5/raylib.git C:/raylib/raylib
git clone https://github.com/raysan5/raygui.git C:/raylib/raygui
make -C C:/raylib/raylib/src PLATFORM=PLATFORM_DESKTOP
make BUILD_MODE=RELEASE
```

构建结果为 `bvhview.exe`。Windows 构建会从 `bvhview.rc` 生成与当前目标架构匹配的资源文件，因此需确保 `windres` 可用。

### macOS

先安装 Xcode Command Line Tools：

```bash
xcode-select --install
```

然后执行：

```bash
git clone https://github.com/raysan5/raylib.git ~/raylib/raylib
git clone https://github.com/raysan5/raygui.git ~/raylib/raygui
make -C ~/raylib/raylib/src PLATFORM=PLATFORM_DESKTOP
make BUILD_MODE=RELEASE PLATFORM_OS=Darwin
```

构建结果为 `bvhview`。

### Linux

以 Debian/Ubuntu 为例，先安装图形和音频开发依赖：

```bash
sudo apt-get update
sudo apt-get install --yes \
  libasound2-dev libgl1-mesa-dev libglu1-mesa-dev libwayland-dev \
  libx11-dev libxcursor-dev libxi-dev libxinerama-dev \
  libxkbcommon-dev libxrandr-dev
```

然后执行：

```bash
git clone https://github.com/raysan5/raylib.git ~/raylib/raylib
git clone https://github.com/raysan5/raygui.git ~/raylib/raygui
make -C ~/raylib/raylib/src PLATFORM=PLATFORM_DESKTOP
make BUILD_MODE=RELEASE
```

构建结果为 `bvhview`。

不传 `BUILD_MODE=RELEASE` 时会构建调试版本。执行 `make clean` 可删除当前平台的可执行文件。

## CI 与发布

GitHub Actions 会在 Pull Request 中构建以下平台组合：

- Windows x64 和 ARM64
- macOS x64 和 ARM64
- Linux x64 和 ARM64

推送 `v*` 格式的 tag 后，所有构建成功时会创建或更新对应的 GitHub Release，并上传六个平台的 ZIP 文件。

## 许可证

本项目采用 [MIT License](LICENSE)。

## 相关链接

- [项目主页](https://theorangeduck.com/page/bvhview)
- [在线演示](https://theorangeduck.com/media/uploads/BVHView/bvhview.html)
