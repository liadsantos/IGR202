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
CMAKE_SOURCE_DIR = /cal/exterieurs/dossantos-22/IGR202/Project/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cal/exterieurs/dossantos-22/IGR202/Project/src/build

# Include any dependencies generated for this target.
include dep/glfw/tests/CMakeFiles/inputlag.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/glfw/tests/CMakeFiles/inputlag.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/glfw/tests/CMakeFiles/inputlag.dir/progress.make

# Include the compile flags for this target's objects.
include dep/glfw/tests/CMakeFiles/inputlag.dir/flags.make

dep/glfw/tests/CMakeFiles/inputlag.dir/inputlag.c.o: dep/glfw/tests/CMakeFiles/inputlag.dir/flags.make
dep/glfw/tests/CMakeFiles/inputlag.dir/inputlag.c.o: /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/tests/inputlag.c
dep/glfw/tests/CMakeFiles/inputlag.dir/inputlag.c.o: dep/glfw/tests/CMakeFiles/inputlag.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/Project/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/glfw/tests/CMakeFiles/inputlag.dir/inputlag.c.o"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/inputlag.dir/inputlag.c.o -MF CMakeFiles/inputlag.dir/inputlag.c.o.d -o CMakeFiles/inputlag.dir/inputlag.c.o -c /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/tests/inputlag.c

dep/glfw/tests/CMakeFiles/inputlag.dir/inputlag.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/inputlag.dir/inputlag.c.i"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/tests/inputlag.c > CMakeFiles/inputlag.dir/inputlag.c.i

dep/glfw/tests/CMakeFiles/inputlag.dir/inputlag.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/inputlag.dir/inputlag.c.s"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/tests/inputlag.c -o CMakeFiles/inputlag.dir/inputlag.c.s

dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/getopt.c.o: dep/glfw/tests/CMakeFiles/inputlag.dir/flags.make
dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/getopt.c.o: /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/deps/getopt.c
dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/getopt.c.o: dep/glfw/tests/CMakeFiles/inputlag.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/Project/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/getopt.c.o"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/getopt.c.o -MF CMakeFiles/inputlag.dir/__/deps/getopt.c.o.d -o CMakeFiles/inputlag.dir/__/deps/getopt.c.o -c /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/deps/getopt.c

dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/getopt.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/inputlag.dir/__/deps/getopt.c.i"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/deps/getopt.c > CMakeFiles/inputlag.dir/__/deps/getopt.c.i

dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/getopt.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/inputlag.dir/__/deps/getopt.c.s"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/deps/getopt.c -o CMakeFiles/inputlag.dir/__/deps/getopt.c.s

dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/glad_gl.c.o: dep/glfw/tests/CMakeFiles/inputlag.dir/flags.make
dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/glad_gl.c.o: /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/deps/glad_gl.c
dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/glad_gl.c.o: dep/glfw/tests/CMakeFiles/inputlag.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cal/exterieurs/dossantos-22/IGR202/Project/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/glad_gl.c.o"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/glad_gl.c.o -MF CMakeFiles/inputlag.dir/__/deps/glad_gl.c.o.d -o CMakeFiles/inputlag.dir/__/deps/glad_gl.c.o -c /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/deps/glad_gl.c

dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/glad_gl.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/inputlag.dir/__/deps/glad_gl.c.i"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/deps/glad_gl.c > CMakeFiles/inputlag.dir/__/deps/glad_gl.c.i

dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/glad_gl.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/inputlag.dir/__/deps/glad_gl.c.s"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/deps/glad_gl.c -o CMakeFiles/inputlag.dir/__/deps/glad_gl.c.s

# Object files for target inputlag
inputlag_OBJECTS = \
"CMakeFiles/inputlag.dir/inputlag.c.o" \
"CMakeFiles/inputlag.dir/__/deps/getopt.c.o" \
"CMakeFiles/inputlag.dir/__/deps/glad_gl.c.o"

# External object files for target inputlag
inputlag_EXTERNAL_OBJECTS =

dep/glfw/tests/inputlag: dep/glfw/tests/CMakeFiles/inputlag.dir/inputlag.c.o
dep/glfw/tests/inputlag: dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/getopt.c.o
dep/glfw/tests/inputlag: dep/glfw/tests/CMakeFiles/inputlag.dir/__/deps/glad_gl.c.o
dep/glfw/tests/inputlag: dep/glfw/tests/CMakeFiles/inputlag.dir/build.make
dep/glfw/tests/inputlag: dep/glfw/src/libglfw3.a
dep/glfw/tests/inputlag: /usr/lib/x86_64-linux-gnu/libm.so
dep/glfw/tests/inputlag: /usr/lib/x86_64-linux-gnu/librt.so
dep/glfw/tests/inputlag: /usr/lib/x86_64-linux-gnu/libm.so
dep/glfw/tests/inputlag: /usr/lib/x86_64-linux-gnu/libX11.so
dep/glfw/tests/inputlag: dep/glfw/tests/CMakeFiles/inputlag.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cal/exterieurs/dossantos-22/IGR202/Project/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking C executable inputlag"
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/inputlag.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/glfw/tests/CMakeFiles/inputlag.dir/build: dep/glfw/tests/inputlag
.PHONY : dep/glfw/tests/CMakeFiles/inputlag.dir/build

dep/glfw/tests/CMakeFiles/inputlag.dir/clean:
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests && $(CMAKE_COMMAND) -P CMakeFiles/inputlag.dir/cmake_clean.cmake
.PHONY : dep/glfw/tests/CMakeFiles/inputlag.dir/clean

dep/glfw/tests/CMakeFiles/inputlag.dir/depend:
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cal/exterieurs/dossantos-22/IGR202/Project/src /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw/tests /cal/exterieurs/dossantos-22/IGR202/Project/src/build /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/tests/CMakeFiles/inputlag.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : dep/glfw/tests/CMakeFiles/inputlag.dir/depend

