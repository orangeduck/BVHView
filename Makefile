# PLATFORM_DESKTOP, PLATFORM_WEB
PLATFORM ?= PLATFORM_DESKTOP

# DEBUG, RELEASE
BUILD_MODE ?= DEBUG

# true, false
INCLUDE_CONSOLE ?= true
ifeq ($(INCLUDE_CONSOLE), true)
    MWINDOWS_FLAG =
else
    MWINDOWS_FLAG = -mwindows
endif

DEFINES = -D _DEFAULT_SOURCE -D RAYLIB_BUILD_MODE=$(BUILD_MODE) -D $(PLATFORM)
PLATFORM_OS ?= $(shell uname)

ifeq ($(PLATFORM), PLATFORM_DESKTOP)
    
    CC = gcc
    
    ifeq ($(findstring Linux,$(PLATFORM_OS)),Linux)
        EXT=
        RAYLIB_DIR = ~/raylib
        INCLUDE_DIR = -I ./ -I ./external -I $(RAYLIB_DIR)/raylib/src -I $(RAYLIB_DIR)/raygui/src
        LIBRARY_DIR = -L $(RAYLIB_DIR)/raylib/src
        ifeq ($(BUILD_MODE),RELEASE)
            CFLAGS ?= $(DEFINES) -Wall -Wno-format-truncation -D NDEBUG -O3 $(INCLUDE_DIR) $(LIBRARY_DIR)
        else
            CFLAGS ?= $(DEFINES) -Wall -Wno-format-truncation -pg $(INCLUDE_DIR) $(LIBRARY_DIR) -ggdb
        endif
        LIBS = -lraylib -lGL -lm
    endif
    
    ifeq ($(findstring MINGW,$(PLATFORM_OS)),MINGW)
        EXT = .exe
        RAYLIB_DIR = C:/raylib
        INCLUDE_DIR = -I ./ -I ./external -I $(RAYLIB_DIR)/raylib/src -I $(RAYLIB_DIR)/raygui/src
        LIBRARY_DIR = -L $(RAYLIB_DIR)/raylib/src
        ifeq ($(BUILD_MODE),RELEASE)
            CFLAGS ?= bvhview.res $(DEFINES) -Wall $(MWINDOWS_FLAG) -D NDEBUG -O3 $(INCLUDE_DIR) $(LIBRARY_DIR)
        else
            CFLAGS ?= bvhview.res $(DEFINES) -Wall -pg $(INCLUDE_DIR) $(LIBRARY_DIR) -ggdb
        endif
        LIBS = -lraylib -lopengl32 -lgdi32 -lwinmm
    endif

endif

ifeq ($(PLATFORM),PLATFORM_WEB)
    CC = emcc
    EXT = .html
    RAYLIB_DIR = C:/raylib
    INCLUDE_DIR = -I ./ -I ./external -I $(RAYLIB_DIR)/raylib/src -I $(RAYLIB_DIR)/raygui/src
    LIBRARY_DIR = -L $(RAYLIB_DIR)/raylib/src
    ifeq ($(BUILD_MODE),RELEASE)
        CFLAGS ?= $(DEFINES) $(RAYLIB_DIR)/raylib/src/libraylib.a -Os -s USE_GLFW=3 -s FORCE_FILESYSTEM=1 -s MAX_WEBGL_VERSION=2 -s ALLOW_MEMORY_GROWTH=1 --shell-file ./shell.html $(INCLUDE_DIR) $(LIBRARY_DIR)
    else
        CFLAGS ?= $(DEFINES) $(RAYLIB_DIR)/raylib/src/libraylib.a -Os -s ASSERTIONS -s USE_GLFW=3 -s FORCE_FILESYSTEM=1 -s MAX_WEBGL_VERSION=2 -s ALLOW_MEMORY_GROWTH=1 --shell-file ./shell.html $(INCLUDE_DIR) $(LIBRARY_DIR)
    endif
endif

.PHONY: all

all: bvhview

bvhview: bvhview.c additions.c external/cwalk.c
	mkdir -p build/BVHView
	cp -r assets build/BVHView
	cp -r samples build/BVHView
	cp example.ps1 build/BVHView
	cp example-with-audio.ps1 build/BVHView
	$(CC) -o build/BVHView/$@$(EXT) bvhview.c additions.c external/cwalk.c $(CFLAGS) $(LIBS)

clean:
	rm bvhview$(EXT)
