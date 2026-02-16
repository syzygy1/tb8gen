FLAGS = -DMAGIC

# Use libzstd instead of LZ4 for compression of temporary files, and set
# compression level.
FLAGS += -DUSE_ZSTD -DZSTD_LEVEL=6
LDFLAGS += -lzstd -lm

# Specify number of threads used for (de)compression of temporary files.
FLAGS += -DCOMPRESSION_THREADS=12

debug = no

# Specify the compiler
# regular linux
CC = gcc -pthread
# cross compiling from Linux to Windows
#CC = x86_64-w64-mingw32-gcc

# Compilation flags
CFLAGS = -march=native -pipe -Wall -Wpedantic -std=c23
#CFLAGS = -march=native -pipe -Wall -Wpedantic -std=c23 -fsanitize=undefined

ifeq ($(debug),no)
	CFLAGS += -O3 -flto=auto -DNDEBUG
else
	CFLAGS += -g3
endif

export FLAGS CC CFLAGS LDFLAGS

tb8gen clean:
	@$(MAKE) -f Makefile.regular $@

all: tb8gen

.PHONY: tb8gen clean
