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
include dep/glfw/tests/CMakeFiles/monitors.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/glfw/tests/CMakeFiles/monitors.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/glfw/tests/CMakeFiles/monitors.dir/progress.make

# Include the compile flags for this target's objects.
include dep/glfw/tests/CMakeFiles/monitors.dir/flags.make

dep/glfw/tests/CMakeFiles/monitors.dir/monitors.c.o: dep/glfw/tests/CMakeFiles/monitors.dir/flags.make
dep/glfw/tests/CMakeFiles/monitors.dir/monitors.c.o: /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests/monitors.c
dep/glfw/tests/CMakeFiles/monitors.dir/monitors.c.o: dep/glfw/tests/CMakeFiles/monitors.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/glfw/tests/CMakeFiles/monitors.dir/monitors.c.o"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/monitors.dir/monitors.c.o -MF CMakeFiles/monitors.dir/monitors.c.o.d -o CMakeFiles/monitors.dir/monitors.c.o -c /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests/monitors.c

dep/glfw/tests/CMakeFiles/monitors.dir/monitors.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/monitors.dir/monitors.c.i"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests/monitors.c > CMakeFiles/monitors.dir/monitors.c.i

dep/glfw/tests/CMakeFiles/monitors.dir/monitors.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/monitors.dir/monitors.c.s"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests/monitors.c -o CMakeFiles/monitors.dir/monitors.c.s

dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/getopt.c.o: dep/glfw/tests/CMakeFiles/monitors.dir/flags.make
dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/getopt.c.o: /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/getopt.c
dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/getopt.c.o: dep/glfw/tests/CMakeFiles/monitors.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/getopt.c.o"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/getopt.c.o -MF CMakeFiles/monitors.dir/__/deps/getopt.c.o.d -o CMakeFiles/monitors.dir/__/deps/getopt.c.o -c /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/getopt.c

dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/getopt.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/monitors.dir/__/deps/getopt.c.i"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/getopt.c > CMakeFiles/monitors.dir/__/deps/getopt.c.i

dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/getopt.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/monitors.dir/__/deps/getopt.c.s"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/getopt.c -o CMakeFiles/monitors.dir/__/deps/getopt.c.s

dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/glad_gl.c.o: dep/glfw/tests/CMakeFiles/monitors.dir/flags.make
dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/glad_gl.c.o: /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/glad_gl.c
dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/glad_gl.c.o: dep/glfw/tests/CMakeFiles/monitors.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/glad_gl.c.o"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/glad_gl.c.o -MF CMakeFiles/monitors.dir/__/deps/glad_gl.c.o.d -o CMakeFiles/monitors.dir/__/deps/glad_gl.c.o -c /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/glad_gl.c

dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/glad_gl.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/monitors.dir/__/deps/glad_gl.c.i"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/glad_gl.c > CMakeFiles/monitors.dir/__/deps/glad_gl.c.i

dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/glad_gl.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/monitors.dir/__/deps/glad_gl.c.s"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/deps/glad_gl.c -o CMakeFiles/monitors.dir/__/deps/glad_gl.c.s

# Object files for target monitors
monitors_OBJECTS = \
"CMakeFiles/monitors.dir/monitors.c.o" \
"CMakeFiles/monitors.dir/__/deps/getopt.c.o" \
"CMakeFiles/monitors.dir/__/deps/glad_gl.c.o"

# External object files for target monitors
monitors_EXTERNAL_OBJECTS =

dep/glfw/tests/monitors: dep/glfw/tests/CMakeFiles/monitors.dir/monitors.c.o
dep/glfw/tests/monitors: dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/getopt.c.o
dep/glfw/tests/monitors: dep/glfw/tests/CMakeFiles/monitors.dir/__/deps/glad_gl.c.o
dep/glfw/tests/monitors: dep/glfw/tests/CMakeFiles/monitors.dir/build.make
dep/glfw/tests/monitors: dep/glfw/src/libglfw3.a
dep/glfw/tests/monitors: /usr/lib/x86_64-linux-gnu/libm.so
dep/glfw/tests/monitors: /usr/lib/x86_64-linux-gnu/librt.so
dep/glfw/tests/monitors: /usr/lib/x86_64-linux-gnu/libm.so
dep/glfw/tests/monitors: /usr/lib/x86_64-linux-gnu/libX11.so
dep/glfw/tests/monitors: dep/glfw/tests/CMakeFiles/monitors.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking C executable monitors"
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/monitors.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/glfw/tests/CMakeFiles/monitors.dir/build: dep/glfw/tests/monitors
.PHONY : dep/glfw/tests/CMakeFiles/monitors.dir/build

dep/glfw/tests/CMakeFiles/monitors.dir/clean:
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests && $(CMAKE_COMMAND) -P CMakeFiles/monitors.dir/cmake_clean.cmake
.PHONY : dep/glfw/tests/CMakeFiles/monitors.dir/clean

dep/glfw/tests/CMakeFiles/monitors.dir/depend:
	cd /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/dep/glfw/tests /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests /cal/exterieurs/dossantos-22/IGR202/TP04-FluidSPH/src/build/dep/glfw/tests/CMakeFiles/monitors.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : dep/glfw/tests/CMakeFiles/monitors.dir/depend

