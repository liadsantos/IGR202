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

# Utility rule file for uninstall.

# Include any custom commands dependencies for this target.
include dep/glfw/CMakeFiles/uninstall.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/glfw/CMakeFiles/uninstall.dir/progress.make

dep/glfw/CMakeFiles/uninstall:
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw && /usr/bin/cmake -P /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/cmake_uninstall.cmake

uninstall: dep/glfw/CMakeFiles/uninstall
uninstall: dep/glfw/CMakeFiles/uninstall.dir/build.make
.PHONY : uninstall

# Rule to build all files generated by this target.
dep/glfw/CMakeFiles/uninstall.dir/build: uninstall
.PHONY : dep/glfw/CMakeFiles/uninstall.dir/build

dep/glfw/CMakeFiles/uninstall.dir/clean:
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw && $(CMAKE_COMMAND) -P CMakeFiles/uninstall.dir/cmake_clean.cmake
.PHONY : dep/glfw/CMakeFiles/uninstall.dir/clean

dep/glfw/CMakeFiles/uninstall.dir/depend:
	cd /cal/exterieurs/dossantos-22/IGR202/Project/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cal/exterieurs/dossantos-22/IGR202/Project/src /cal/exterieurs/dossantos-22/IGR202/Project/src/dep/glfw /cal/exterieurs/dossantos-22/IGR202/Project/src/build /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw /cal/exterieurs/dossantos-22/IGR202/Project/src/build/dep/glfw/CMakeFiles/uninstall.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : dep/glfw/CMakeFiles/uninstall.dir/depend
