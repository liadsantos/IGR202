# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build

# Include any dependencies generated for this target.
include CMakeFiles/tpSph.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/tpSph.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/tpSph.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tpSph.dir/flags.make

CMakeFiles/tpSph.dir/main.cpp.o: CMakeFiles/tpSph.dir/flags.make
CMakeFiles/tpSph.dir/main.cpp.o: /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/main.cpp
CMakeFiles/tpSph.dir/main.cpp.o: CMakeFiles/tpSph.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tpSph.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tpSph.dir/main.cpp.o -MF CMakeFiles/tpSph.dir/main.cpp.o.d -o CMakeFiles/tpSph.dir/main.cpp.o -c /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/main.cpp

CMakeFiles/tpSph.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tpSph.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/main.cpp > CMakeFiles/tpSph.dir/main.cpp.i

CMakeFiles/tpSph.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tpSph.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/main.cpp -o CMakeFiles/tpSph.dir/main.cpp.s

CMakeFiles/tpSph.dir/dep/glad/src/glad.c.o: CMakeFiles/tpSph.dir/flags.make
CMakeFiles/tpSph.dir/dep/glad/src/glad.c.o: /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glad/src/glad.c
CMakeFiles/tpSph.dir/dep/glad/src/glad.c.o: CMakeFiles/tpSph.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/tpSph.dir/dep/glad/src/glad.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/tpSph.dir/dep/glad/src/glad.c.o -MF CMakeFiles/tpSph.dir/dep/glad/src/glad.c.o.d -o CMakeFiles/tpSph.dir/dep/glad/src/glad.c.o -c /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glad/src/glad.c

CMakeFiles/tpSph.dir/dep/glad/src/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tpSph.dir/dep/glad/src/glad.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glad/src/glad.c > CMakeFiles/tpSph.dir/dep/glad/src/glad.c.i

CMakeFiles/tpSph.dir/dep/glad/src/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tpSph.dir/dep/glad/src/glad.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glad/src/glad.c -o CMakeFiles/tpSph.dir/dep/glad/src/glad.c.s

# Object files for target tpSph
tpSph_OBJECTS = \
"CMakeFiles/tpSph.dir/main.cpp.o" \
"CMakeFiles/tpSph.dir/dep/glad/src/glad.c.o"

# External object files for target tpSph
tpSph_EXTERNAL_OBJECTS =

tpSph: CMakeFiles/tpSph.dir/main.cpp.o
tpSph: CMakeFiles/tpSph.dir/dep/glad/src/glad.c.o
tpSph: CMakeFiles/tpSph.dir/build.make
tpSph: dep/glfw/src/libglfw3.a
tpSph: /usr/lib/x86_64-linux-gnu/librt.so
tpSph: /usr/lib/x86_64-linux-gnu/libm.so
tpSph: /usr/lib/x86_64-linux-gnu/libX11.so
tpSph: CMakeFiles/tpSph.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable tpSph"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tpSph.dir/link.txt --verbose=$(VERBOSE)
	/usr/bin/cmake -E copy /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/tpSph /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src

# Rule to build all files generated by this target.
CMakeFiles/tpSph.dir/build: tpSph
.PHONY : CMakeFiles/tpSph.dir/build

CMakeFiles/tpSph.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tpSph.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tpSph.dir/clean

CMakeFiles/tpSph.dir/depend:
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles/tpSph.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tpSph.dir/depend

