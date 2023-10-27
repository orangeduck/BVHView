PLATFORM ?= PLATFORM_DESKTOP
BUILD_MODE ?= DEBUG
RAYLIB_DIR = C:/raylib
INCLUDE_DIR = -I ./ -I $(RAYLIB_DIR)/raylib/src -I $(RAYLIB_DIR)/raygui/src
LIBRARY_DIR = -L $(RAYLIB_DIR)/raylib/src
DEFINES = -D _DEFAULT_SOURCE -D RAYLIB_BUILD_MODE=$(BUILD_MODE) -D $(PLATFORM)

ifeq ($(PLATFORM),PLATFORM_DESKTOP)
    CC = gcc
    EXT = .exe
    ifeq ($(BUILD_MODE),RELEASE)
        CFLAGS ?= bvhview.res $(DEFINES) -Wall -mwindows -D NDEBUG -O3 $(INCLUDE_DIR) $(LIBRARY_DIR) 
	else
        CFLAGS ?= bvhview.res $(DEFINES) -Wall -g $(INCLUDE_DIR) $(LIBRARY_DIR) 
	endif
    LIBS = -lraylib -lopengl32 -lgdi32 -lwinmm
endif

ifeq ($(PLATFORM),PLATFORM_WEB)
    CC = emcc
    EXT = .html
    ifeq ($(BUILD_MODE),RELEASE)
        CFLAGS ?= $(DEFINES) $(RAYLIB_DIR)/raylib/src/libraylib.a -Os -s USE_GLFW=3 -s FORCE_FILESYSTEM=1 -s MAX_WEBGL_VERSION=2 -s ALLOW_MEMORY_GROWTH=1 --shell-file ./shell.html $(INCLUDE_DIR) $(LIBRARY_DIR)
    else
        CFLAGS ?= $(DEFINES) $(RAYLIB_DIR)/raylib/src/libraylib.a -Os -s ASSERTIONS -s USE_GLFW=3 -s FORCE_FILESYSTEM=1 -s MAX_WEBGL_VERSION=2 -s ALLOW_MEMORY_GROWTH=1 --shell-file ./shell.html $(INCLUDE_DIR) $(LIBRARY_DIR)
    endif
endif

SOURCE = $(wildcard *.c)
HEADER = $(wildcard *.h)

.PHONY: all

all: bvhview

bvhview: $(SOURCE) $(HEADER)
	$(CC) -o $@$(EXT) $(SOURCE) $(CFLAGS) $(LIBS) 

clean:
	rm bvhview$(EXT)
