# MAKEFILE to NML v1.0.0
#
# Author: Yoel Monsalve | yymonsalve@gmail.com
# Date  : May, 2022. 

# flags
CC := gcc
CSTD := gnu99
CCFLAGS := -Wall -c -lm --std=$(CSTD)
CCFLAGS_EXAMPLES := -Wall -lm --std=$(CSTD)
CCFLAGS_TESTS    := -Wall -lm --std=$(CSTD)
AR := ar
ARFLAGS := crs
SRC_DIR := src
SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJ_DIR := obj
OBJ_NAMES := $(SRCS:$(SRC_DIR)/%.c=%.o)
OBJS := $(addprefix $(OBJ_DIR)/,$(OBJ_NAMES))
INC_DIR := include

# NOTE: if 'include/'' contains subdirectories, something more complex will be needed,
#       for example
#       
#       INCS := $(shell find $(INC_DIR) -type f -name '*.h')
INCS := $(shell find $(INC_DIR) -type f -name '*.h')
INC_FLAGS := $(addprefix -I,$(INC_DIR))

LIB_DIR  := ./lib#                                 directory for library
# NOTE: if 'lib/' contains subdirectories, something more complex will be needed,
#       for example
#       
#       LIBS := $(shell find $(LIB_DIR) -type f -name '*.a')
LIB_NAME := libnml.a#                              filename for the library
LIB_NAME_SIMPLE := nml#                            simple name of the library
LIBS     := $(addprefix $(LIB_DIR)/,$(LIB_NAME))#  full name of the library, e.g., ./libnml.a

DIST_DIR := ./dist
DISTS    := $(addprefix $(DIST_DIR)/,$(LIBS))
DISTS    += $(addprefix $(DIST_DIR)/,$(INCS))

EXAMPLE_DIR := ./examples
EXAMPLE_LIB := $(EXAMPLE_DIR)/lib
EXAMPLE_INC := $(EXAMPLE_DIR)/include
TEST_DIR    := ./tests
TEST_LIB    := $(TEST_DIR)/lib
TEST_INC    := $(TEST_DIR)/include

# color constants
RED:=\033[0;31m
GREEN:=\033[0;32m
YELLOW:=\033[0;33m
BLUE:=\033[0;34m
NC:=\033[0m

.PHONY: all clean print_info objects dist lib examples tests
#.SILENT: $(LIBS)

# PRIMARY TARGET
all: lib dist examples tests

print_info:
	@echo "${YELLOW}Files are:${NC}"
	@echo "${GREEN}SRCS${NC}: ${SRCS}"
	@echo "${GREEN}OBJS${NC}: ${OBJS}"
	@echo "${GREEN}INCS${NC}: ${INCS}"

# DISTRIBUTABLE FILES
#   libraries, headers
# ======================================================
dist: $(DISTS)

# build directory for dist, if it does not exist
dist_dir:
	@[ -d $(DIST_DIR) ] || { echo "${YELLOW}Creating '$(DIST_DIR)/'${NC}" \
		&& mkdir $(DIST_DIR); }

$(DISTS): $(LIBS) $(INCS)
	@[ -d $(DIST_DIR) ] || mkdir -p $(DIST_DIR)
	cp -rv ${LIB_DIR} ${DIST_DIR}/
	cp -rv ${INC_DIR} ${DIST_DIR}/
	
# NML LIBRARY
# =============================================================================
# create directory for objects
$(OBJ_DIR):
	@echo "Creating directory '${OBJ_DIR}/'"
	$(shell [ -d $(OBJ_DIR) ] || mkdir $(OBJ_DIR))

# compile objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(INCS)
	$(CC) $(CCFLAGS) -o $@ $< $(INC_FLAGS)

# build library
$(LIBS): $(OBJS)
	@printf "${YELLOW}Building Library${NC} -> ${LIB_NAME}\nFrom objects:\n"    # printf, instead of echo, for a more consistent behavior respect to \
	                                                                            # escaped characters, like '\n'
	@for object in ${OBJS}; \
	do \
		printf "\t$$object\n"; \
	done
	${AR} ${ARFLAGS} ${LIB_DIR}/${LIB_NAME} ${OBJS}
	@echo "${YELLOW}File created -> ${LIB_DIR}/${LIB_NAME}${NC}"

lib: lib_dir $(LIBS) 
	@echo "${YELLOW}Library done.${NC}"

# build directory for library, if it does not exist
lib_dir: $(OBJ_DIR)
	@[ -d $(LIB_DIR) ] || { echo "${YELLOW}Creating '$(LIB_DIR)/'${NC}" \
		&& mkdir $(LIB_DIR); }

objects: $(OBJ_DIR) $(OBJS)
	@echo "${YELLOW}Compiling objects done.${NC}"


# EXAMPLES 
# ======================================================
EXAMPLE_SRCS  := $(wildcard $(EXAMPLE_DIR)/src/*.c)
EXAMPLE_EXECS := $(notdir $(EXAMPLE_SRCS))
EXAMPLE_EXECS := $(addprefix $(EXAMPLE_DIR)/,$(EXAMPLE_EXECS:%.c=%))
EXAMPLE_DISTS := $(DISTS:$(DIST_DIR)/%=$(EXAMPLE_DIR)/%)

examples: example_dir $(EXAMPLE_EXECS)

# compiling executables
$(EXAMPLE_DIR)/%: $(EXAMPLE_DIR)/src/%.c $(EXAMPLE_DISTS)
	@echo "\t$(<) -> ${GREEN}$(@)${NC}"
	${CC} ${CCFLAGS} $< -L${EXAMPLE_LIB} -l${LIB_NAME_SIMPLE} -I${EXAMPLE_INC} -o $@

# prepare directory for examples
example_dir:
	@echo "${YELLOW}Preparing '$(EXAMPLE_DIR)/' folder with the latest version:${NC}"
	@[ -d $(EXAMPLE_DIR) ] || { echo "${YELLOW}Creating '$(EXAMPLE_DIR)/'${NC}" \
		&& mkdir $(EXAMPLE_DIR); }

# prepare library for examples
$(EXAMPLE_DISTS): $(DISTS)
	@echo "\tCopying ${DIST_DIR}/* to ${EXAMPLE_DIR}/*"
	cp -rv ${DIST_DIR}/* "${EXAMPLE_DIR}/"


# TESTS
# ======================================================
TEST_SRCS  := $(wildcard $(TEST_DIR)/src/*.c)
TEST_EXECS := $(notdir $(TEST_SRCS))
TEST_EXECS := $(addprefix $(TEST_DIR)/,$(TEST_EXECS:%.c=%))
TEST_DISTS := $(DISTS:$(DIST_DIR)/%=$(TEST_DIR)/%)

tests: test_dir $(TEST_EXECS)

# compiling executables
$(TEST_DIR)/%: $(TEST_DIR)/src/%.c $(TEST_DISTS)
	@echo "\t$(<) -> ${GREEN}$(@)${NC}"
	${CC} ${CCFLAGS_TESTS} $< -L${TEST_LIB} -l${LIB_NAME_SIMPLE} -I${TEST_INC} -o $@

# prepare directory for examples
test_dir:
	@echo "${YELLOW}Preparing '$(TEST_DIR)/' folder with the latest version:${NC}"
	@[ -d $(TEST_DIR) ] || { echo "${YELLOW}Creating '$(TEST_DIR)/'${NC}" \
		&& mkdir $(TEST_DIR); }

# prepare library for examples
$(TEST_DISTS): $(DISTS)
	@echo "\tCopying ${DIST_DIR}/* to ${TEST_DIR}/*"
	cp -rv ${DIST_DIR}/* "${TEST_DIR}/"


# OTHERS
# ===================================================================
check:
	@echo "${YELLOW}Sanity Checks:${NC}"
	# TO - DO

greet:
	@echo
	@echo "${GREEN}N${NC}eat ${GREEN}M${NC}atrix ${GREEN}L${NC}ibrary (libnml)"
	@echo

clean:
	# objects
	@echo "${YELLOW}Deleting:${NC}"
	@echo "\tObject files (*.o) and directory '$(OBJ_DIR)'"
	rm -rf $(OBJ_DIR)

	# lib
	[ -d $(LIB_DIR) ] && rm -r $(LIB_DIR) || true

	# examples
	@echo "\tClean '${EXAMPLE_DIR}'"
	@for file in `ls $(EXAMPLE_DIR)/*`; do \
		{ [ -e $${file} ] && rm $${file}; } || true; \
	done
	[ -d $(EXAMPLE_DIR)/include ] && rm -r $(EXAMPLE_DIR)/include || true
	[ -d $(EXAMPLE_DIR)/lib ] && rm -r $(EXAMPLE_DIR)/lib || true

	# tests
	@echo "\tClean '${TEST_DIR}'"
	@for file in `ls $(TEST_DIR)/*`; do \
		{ [ -e $${file} ] && rm $${file}; } || true; \
	done
	[ -d $(TEST_DIR)/include ] && rm -r $(TEST_DIR)/include || true
	[ -d $(TEST_DIR)/lib ] && rm -r $(TEST_DIR)/lib || true

	@echo "${YELLOW}Clean done${NC}"
