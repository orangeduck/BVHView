PLATFORM ?= PLATFORM_DESKTOP
BUILD_MODE ?= DEBUG
DEFINES = -D _DEFAULT_SOURCE -D RAYLIB_BUILD_MODE=$(BUILD_MODE) -D $(PLATFORM)
PLATFORM_OS ?= $(shell uname)

ifeq ($(PLATFORM),PLATFORM_DESKTOP)
    
    CC = gcc

    ifeq ($(PLATFORM_OS),Darwin) # macOS settings
        EXT=
        RAYLIB_DIR = ~/raylib
        INCLUDE_DIR = -I ./ -I $(RAYLIB_DIR)/raylib/src -I $(RAYLIB_DIR)/raygui/src
        LIBRARY_DIR = -L $(RAYLIB_DIR)/raylib/src
        ifeq ($(BUILD_MODE),RELEASE)
            CFLAGS = -D_DEFAULT_SOURCE -DPLATFORM_DESKTOP -Wall -D NDEBUG -O3 $(INCLUDE_DIR) $(LIBRARY_DIR)
        else
            CFLAGS = -D_DEFAULT_SOURCE -DPLATFORM_DESKTOP -Wall -g $(INCLUDE_DIR) $(LIBRARY_DIR)
        endif
        LIBS = -lraylib -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo -framework QuartzCore
    endif

    ifeq ($(findstring Linux,$(PLATFORM_OS)),Linux)
        EXT=
        RAYLIB_DIR = ~/raylib
        INCLUDE_DIR = -I ./ -I $(RAYLIB_DIR)/raylib/src -I $(RAYLIB_DIR)/raygui/src
        LIBRARY_DIR = -L $(RAYLIB_DIR)/raylib/src
        ifeq ($(BUILD_MODE),RELEASE)
            CFLAGS ?= $(DEFINES) -Wall -Wno-format-truncation -D NDEBUG -O3 $(INCLUDE_DIR) $(LIBRARY_DIR)
        else
            CFLAGS ?= $(DEFINES) -Wall -Wno-format-truncation -g $(INCLUDE_DIR) $(LIBRARY_DIR)
        endif
        LIBS = -lraylib -lGL -lm -lpthread -ldl -lrt -lX11
    endif
    
    ifneq ($(filter MINGW% MSYS%,$(PLATFORM_OS)),)
        EXT = .exe
        RESOURCE = bvhview.res
        WINDRES ?= windres
        RAYLIB_DIR = C:/raylib
        INCLUDE_DIR = -I ./ -I $(RAYLIB_DIR)/raylib/src -I $(RAYLIB_DIR)/raygui/src
        LIBRARY_DIR = -L $(RAYLIB_DIR)/raylib/src
        ifeq ($(BUILD_MODE),RELEASE)
            CFLAGS ?= $(DEFINES) -Wall -mwindows -D NDEBUG -O3 $(INCLUDE_DIR) $(LIBRARY_DIR)
        else
            CFLAGS ?= $(DEFINES) -Wall -g $(INCLUDE_DIR) $(LIBRARY_DIR)
        endif
        LIBS = -lraylib -lopengl32 -lgdi32 -lwinmm
    endif

endif

ifeq ($(PLATFORM),PLATFORM_WEB)
    CC = emcc
    EXT = .html
    RAYLIB_DIR = C:/raylib
    INCLUDE_DIR = -I ./ -I $(RAYLIB_DIR)/raylib/src -I $(RAYLIB_DIR)/raygui/src
    LIBRARY_DIR = -L $(RAYLIB_DIR)/raylib/src
    ifeq ($(BUILD_MODE),RELEASE)
        CFLAGS ?= $(DEFINES) $(RAYLIB_DIR)/raylib/src/libraylib.web.a -Os -s USE_GLFW=3 -s FORCE_FILESYSTEM=1 -s MAX_WEBGL_VERSION=2 -s ALLOW_MEMORY_GROWTH=1 --shell-file ./shell.html $(INCLUDE_DIR) $(LIBRARY_DIR)
    else
        CFLAGS ?= $(DEFINES) $(RAYLIB_DIR)/raylib/src/libraylib.web.a -Os -s ASSERTIONS -s USE_GLFW=3 -s FORCE_FILESYSTEM=1 -s MAX_WEBGL_VERSION=2 -s ALLOW_MEMORY_GROWTH=1 --shell-file ./shell.html $(INCLUDE_DIR) $(LIBRARY_DIR)
    endif
endif

.PHONY: all $(RESOURCE)

all: bvhview

bvhview: bvhview.c $(RESOURCE)
	$(CC) -o $@$(EXT) bvhview.c $(CFLAGS) $(LIBS)

bvhview.res: bvhview.rc bvhview.ico
	$(WINDRES) -O coff -o $@ $<

clean:
	rm bvhview$(EXT)

