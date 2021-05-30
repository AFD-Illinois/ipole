# Problem to compile
MODEL = iharm

# Top directory of HDF5, or blank if using h5pcc
HDF5_DIR =
# Top directory of GSL, or blank if installed to system
GSL_DIR =
# System /lib equivalent (can be /usr/lib, /lib64, /usr/lib64)
# Can leave this blank if it's included automatically by GCC
SYSTEM_LIBDIR = /lib64

# Try pointing this to h5pcc or h5cc on your machine, before hunting down libraries
CC=h5cc
# Example CFLAGS for going fast with GCC
CFLAGS = -std=gnu99 -O3 -march=native -mtune=native -flto -fopenmp -funroll-loops
MATH_LIB = -lm
# ICC does not like -lm and uses different flags
#CFLAGS = -xCORE-AVX2 -Ofast -fstrict-aliasing -Wall -Werror -ipo -qopenmp
#MATH_LIB =

# Name of the executable
EXE = ipole

# Executables which apparently aren't standard
MD5=md5sum
ECHO=echo -e

# Overrides of the above for macOS
ifneq (,$(findstring Darwin,$(shell uname)))
	export HDF5_CC = /usr/local/opt/llvm/bin/clang
	export HDF5_CLINKER = /usr/local/opt/llvm/bin/clang

	GSL_DIR=/usr/local
	SYSTEM_LIBDIR=

	MD5=md5
	ECHO=echo
endif

# Override these defaults if we know the machine we're working with
# Once you know what compiles, add it as a machine def here
MAKEFILE_PATH := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
HOST := $(shell hostname)
ifneq (,$(findstring stampede2,$(HOST)))
	-include $(MAKEFILE_PATH)/machines/stampede2.make
endif
ifneq (,$(findstring frontera,$(HOST)))
        -include $(MAKEFILE_PATH)/machines/frontera.make
endif
# Hack to check only whether host begins with bh*
ifneq (,$(findstring beginsbh,begins$(HOST)))
        -include $(MAKEFILE_PATH)/machines/bh-cluster.make
endif
-include $(MAKEFILE_PATH)/machines/$(HOST).make

# Allow overrides of which cflags to add
CFLAGS += $(CFLAGS_CUSTOM)

# Everything below this should be static

## VERSION PRESERVATION ##
GIT_VERSION := $(shell cd $(MAKEFILE_PATH); git describe --dirty --always --tags)

## LINKING PARAMETERS ##
LINK = $(CC)
LDFLAGS = $(CFLAGS)

HDF5_LIB = -lhdf5_hl -lhdf5
GSL_LIB = -lgsl -lgslcblas

## LOGIC FOR PATHS ##
CORE_DIR := $(MAKEFILE_PATH)/src/
EMIS_DIR := $(MAKEFILE_PATH)/src/symphony/
MODEL_DIR := $(MAKEFILE_PATH)/model/$(MODEL)/
VPATH = $(CORE_DIR):$(EMIS_DIR):$(MODEL_DIR)

#ARC_DIR := $(MAKEFILE_PATH)/model/$(MODEL)/build_archive/
# TODO this is I think gmake-specific
ARC_DIR := $(CURDIR)/build_archive/

SRC := $(wildcard $(CORE_DIR)/*.c) $(wildcard $(EMIS_DIR)/*.c) $(wildcard $(MODEL_DIR)/*.c)
HEAD := $(wildcard $(CORE_DIR)/*.h) $(wildcard $(EMIS_DIR)/*.h) $(wildcard $(MODEL_DIR)/*.h)

HEAD_ARC := $(addprefix $(ARC_DIR)/, $(notdir $(HEAD)))
OBJ := $(addprefix $(ARC_DIR)/, $(notdir $(SRC:%.c=%.o)))

INC = -I$(ARC_DIR)
LIBDIR =
LIB = $(MATH_LIB) $(GSL_LIB)

# Add HDF and MPI directories only if compiler doesn't
ifneq ($(strip $(HDF5_DIR)),)
	INC += -I$(HDF5_DIR)/include/
	LIBDIR += -L$(HDF5_DIR)/lib/
	LIB += $(HDF5_LIB)
endif
ifneq ($(strip $(GSL_DIR)),)
	INC += -I$(GSL_DIR)/include/
	LIBDIR += -L$(GSL_DIR)/lib/
endif
ifneq ($(strip $(SYSTEM_LIBDIR)),)
	# Prefer user libraries (above) to system
	LIBDIR += -L$(SYSTEM_LIBDIR)
endif

## TARGETS ##
.PRECIOUS: $(ARC_DIR)/$(EXE) $(ARC_DIR)/%

default: build

build: $(EXE)
	@$(ECHO) "Completed build with model: $(MODEL)"
	@$(ECHO) "CFLAGS: $(CFLAGS)"
	@$(ECHO) "MD5: $(shell $(MD5) $(EXE))"

debug: CFLAGS += -g -Wall -Werror -Wno-unused-variable -Wno-unused-but-set-variable
debug: CFLAGS += -DDEBUG=1
debug: build

profile: CFLAGS += -g -pg
profile: build

vtune: CFLAGS += -g -Wall -Werror
vtune: CFLAGS += -debug inline-debug-info -shared-intel
vtune: build


clean:
	@$(ECHO) "Cleaning build files..."
	@rm -rf $(EXE) $(OBJ) $(ARC_DIR)

$(EXE): $(ARC_DIR)/$(EXE)
	@cp $(ARC_DIR)/$(EXE) .

$(ARC_DIR)/$(EXE): $(OBJ)
	@$(ECHO) "\tLinking $(EXE)"
	@$(LINK) $(LDFLAGS) $(OBJ) $(LIBDIR) $(LIB) -o $(ARC_DIR)/$(EXE)
	@rm $(OBJ) # This ensures full recompile

$(ARC_DIR)/%.o: $(ARC_DIR)/%.c $(HEAD_ARC)
	@$(ECHO) "\tCompiling $(notdir $<)"
	@$(CC) $(CFLAGS) $(INC) -DVERSION=$(GIT_VERSION) -DNOTES=$(NOTES) -DMODEL=$(MODEL) -c $< -o $@

$(ARC_DIR)/%: % | $(ARC_DIR)
	@cp $< $(ARC_DIR)

$(ARC_DIR):
	@mkdir $(ARC_DIR)
