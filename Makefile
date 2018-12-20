

# Things needed from the submodules directory (relative to src)
# =====================================================================
THIRD_PARTY += ../third_party/patches/patches.cpp
THIRD_PARTY += ../third_party/patches/patches.hpp
THIRD_PARTY += ../third_party/ndarray/include/ndarray.hpp
THIRD_PARTY += ../third_party/visit_struct/include/visit_struct/visit_struct.hpp


# Build rules
# =====================================================================
-include Makefile.in
.PHONY: all clean

CXXFLAGS += -std=c++14
CXXFLAGS += -MMD -MP
CXXFLAGS += -Wall -Wextra -Wno-missing-braces
CXXFLAGS += -O3

SRC := $(wildcard src/*.cpp)
OBJ := $(SRC:%.cpp=%.o)
DEP := $(SRC:%.cpp=%.d)
EXE := atmo


# Build rules
# =====================================================================
all: post-build

pre-build:
	@cd src; ln -sf $(THIRD_PARTY) .

main-build: pre-build
	@$(MAKE) $(EXE)

post-build: main-build
	@find src -type l -delete

$(EXE): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJ) $(DEP) $(EXE)

-include $(DEP)
