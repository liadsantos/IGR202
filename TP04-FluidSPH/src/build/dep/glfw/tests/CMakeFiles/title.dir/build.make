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
include dep/glfw/tests/CMakeFiles/title.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/glfw/tests/CMakeFiles/title.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/glfw/tests/CMakeFiles/title.dir/progress.make

# Include the compile flags for this target's objects.
include dep/glfw/tests/CMakeFiles/title.dir/flags.make

dep/glfw/tests/CMakeFiles/title.dir/title.c.o: dep/glfw/tests/CMakeFiles/title.dir/flags.make
dep/glfw/tests/CMakeFiles/title.dir/title.c.o: /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests/title.c
dep/glfw/tests/CMakeFiles/title.dir/title.c.o: dep/glfw/tests/CMakeFiles/title.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/glfw/tests/CMakeFiles/title.dir/title.c.o"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/title.dir/title.c.o -MF CMakeFiles/title.dir/title.c.o.d -o CMakeFiles/title.dir/title.c.o -c /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests/title.c

dep/glfw/tests/CMakeFiles/title.dir/title.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/title.dir/title.c.i"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests/title.c > CMakeFiles/title.dir/title.c.i

dep/glfw/tests/CMakeFiles/title.dir/title.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/title.dir/title.c.s"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests/title.c -o CMakeFiles/title.dir/title.c.s

dep/glfw/tests/CMakeFiles/title.dir/__/deps/glad_gl.c.o: dep/glfw/tests/CMakeFiles/title.dir/flags.make
dep/glfw/tests/CMakeFiles/title.dir/__/deps/glad_gl.c.o: /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/glad_gl.c
dep/glfw/tests/CMakeFiles/title.dir/__/deps/glad_gl.c.o: dep/glfw/tests/CMakeFiles/title.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object dep/glfw/tests/CMakeFiles/title.dir/__/deps/glad_gl.c.o"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/title.dir/__/deps/glad_gl.c.o -MF CMakeFiles/title.dir/__/deps/glad_gl.c.o.d -o CMakeFiles/title.dir/__/deps/glad_gl.c.o -c /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/glad_gl.c

dep/glfw/tests/CMakeFiles/title.dir/__/deps/glad_gl.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/title.dir/__/deps/glad_gl.c.i"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/glad_gl.c > CMakeFiles/title.dir/__/deps/glad_gl.c.i

dep/glfw/tests/CMakeFiles/title.dir/__/deps/glad_gl.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/title.dir/__/deps/glad_gl.c.s"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/glad_gl.c -o CMakeFiles/title.dir/__/deps/glad_gl.c.s

# Object files for target title
title_OBJECTS = \
"CMakeFiles/title.dir/title.c.o" \
"CMakeFiles/title.dir/__/deps/glad_gl.c.o"

# External object files for target title
title_EXTERNAL_OBJECTS =

dep/glfw/tests/title: dep/glfw/tests/CMakeFiles/title.dir/title.c.o
dep/glfw/tests/title: dep/glfw/tests/CMakeFiles/title.dir/__/deps/glad_gl.c.o
dep/glfw/tests/title: dep/glfw/tests/CMakeFiles/title.dir/build.make
dep/glfw/tests/title: dep/glfw/src/libglfw3.a
dep/glfw/tests/title: /usr/lib/x86_64-linux-gnu/libm.so
dep/glfw/tests/title: /usr/lib/x86_64-linux-gnu/librt.so
dep/glfw/tests/title: /usr/lib/x86_64-linux-gnu/libm.so
dep/glfw/tests/title: /usr/lib/x86_64-linux-gnu/libX11.so
dep/glfw/tests/title: dep/glfw/tests/CMakeFiles/title.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable title"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/title.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/glfw/tests/CMakeFiles/title.dir/build: dep/glfw/tests/title
.PHONY : dep/glfw/tests/CMakeFiles/title.dir/build

dep/glfw/tests/CMakeFiles/title.dir/clean:
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && $(CMAKE_COMMAND) -P CMakeFiles/title.dir/cmake_clean.cmake
.PHONY : dep/glfw/tests/CMakeFiles/title.dir/clean

dep/glfw/tests/CMakeFiles/title.dir/depend:
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests/CMakeFiles/title.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : dep/glfw/tests/CMakeFiles/title.dir/depend
