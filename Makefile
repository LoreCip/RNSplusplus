# /*************************************************************************
# * MAKEFILE FOR RNS and RNSRegrid                   
# *************************************************************************/

MAIN = RNS_Diff_Jc
SIZE = -DMDIV=401 -DSDIV=801

CC = gcc
SRC_DIR = ./src
BUILD_DIR = $(SRC_DIR)/build
MAIN_DIR = $(SRC_DIR)/mains
# DEBUG = -DDDEBUG

# --- HDF5 DETECTION LOGIC (Priority: Manual > pkg-config > h5cc) ---

ifdef HDF5_DIR
    # 1. Manual Path provided by the user
    HDF5_CFLAGS  := -I$(HDF5_DIR)/include
    HDF5_LDFLAGS := -L$(HDF5_DIR)/lib -lhdf5
else
    # 2. Try using pkg-config
    HDF5_INSTALLED := $(shell pkg-config --exists hdf5 && echo yes)
    ifeq ($(HDF5_INSTALLED),yes)
        HDF5_CFLAGS  := $(shell pkg-config --cflags hdf5)
        HDF5_LDFLAGS := $(shell pkg-config --libs hdf5)
    else
        # 3. Fallback: Try to find h5cc in the PATH
        H5CC_PATH := $(shell which h5cc 2>/dev/null)
        ifneq ($(H5CC_PATH),)
            HDF5_ROOT    := $(shell h5cc -showconfig | grep 'Installation point' | awk '{print $$3}')
            HDF5_CFLAGS  := -I$(HDF5_ROOT)/include
            HDF5_LDFLAGS := -L$(HDF5_ROOT)/lib -lhdf5
        else
            # Error section if all methods fail
            $(info )
            $(info ************************************************************)
            $(info [ERROR] HDF5 library not found.)
            $(info ************************************************************)
            $(info Possible solutions:)
            $(info 1. Install HDF5 dev files)
            $(info 2. Pass the path manually: make HDF5_DIR=/your/path/to/hdf5)
            $(info ************************************************************)
            $(info )
            $(error Compilation aborted)
        endif
    endif
endif

# Append standard system dependencies for HDF5
HDF5_LDFLAGS += -lm -lz -ldl

# --- COMPILATION FLAGS ---

CFLAGS = -I$(SRC_DIR)/include $(HDF5_CFLAGS) -O3 -march=native -flto=auto -ffast-math
LDFLAGS = $(HDF5_LDFLAGS)

# Source files for RNS and RNSRegrid
CMN = nrutil.c parfile.c equil_util.c equil.c output.c main_util.c
SRCS_RNS = mains/$(MAIN).c $(CMN)
SRCS_RNSRegrid = RNSRegrid.c rnsregrid_util.c $(CMN)

# Object files
OBJS_RNS = $(addprefix $(BUILD_DIR)/, $(notdir $(SRCS_RNS:.c=.o)))
OBJS_RNSRegrid = $(addprefix $(BUILD_DIR)/, $(SRCS_RNSRegrid:.c=.o))

# Executables
EXEC_RNS = $(MAIN)
EXEC_RNSRegrid = RNSRegrid

# Default rule
all: $(EXEC_RNS) $(EXEC_RNSRegrid)

# Link RNS
$(EXEC_RNS): $(OBJS_RNS)
	$(CC) $(CFLAGS) $(SIZE) -o $(EXEC_RNS) $(OBJS_RNS) $(LDFLAGS)

# Link RNSRegrid
$(EXEC_RNSRegrid): $(OBJS_RNSRegrid)
	$(CC) $(CFLAGS) $(SIZE) -o $(EXEC_RNSRegrid) $(OBJS_RNSRegrid) $(LDFLAGS)

# Compile .c files from src directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(SIZE) -c $< -o $@

# Compile .c files from mains directory
$(BUILD_DIR)/%.o: $(MAIN_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(SIZE) -c $< -o $@

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean up
clean:
	rm -rf $(BUILD_DIR) $(EXEC_RNS) $(EXEC_RNSRegrid)

.PHONY: all clean debug