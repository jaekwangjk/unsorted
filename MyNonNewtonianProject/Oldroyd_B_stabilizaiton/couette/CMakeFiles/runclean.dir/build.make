# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.15.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.15.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/jaekwangkim/codes/MyAdjProject/Oldroyd_B_stabilizaiton/couette

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jaekwangkim/codes/MyAdjProject/Oldroyd_B_stabilizaiton/couette

# Utility rule file for runclean.

# Include the progress variables for this target.
include CMakeFiles/runclean.dir/progress.make

CMakeFiles/runclean:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/jaekwangkim/codes/MyAdjProject/Oldroyd_B_stabilizaiton/couette/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "runclean invoked"
	/usr/local/Cellar/cmake/3.15.2/bin/cmake -E remove *.log *.gmv *.gnuplot *.gpl *.eps *.pov *.vtk *.ucd *.d2

runclean: CMakeFiles/runclean
runclean: CMakeFiles/runclean.dir/build.make

.PHONY : runclean

# Rule to build all files generated by this target.
CMakeFiles/runclean.dir/build: runclean

.PHONY : CMakeFiles/runclean.dir/build

CMakeFiles/runclean.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/runclean.dir/cmake_clean.cmake
.PHONY : CMakeFiles/runclean.dir/clean

CMakeFiles/runclean.dir/depend:
	cd /Users/jaekwangkim/codes/MyAdjProject/Oldroyd_B_stabilizaiton/couette && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jaekwangkim/codes/MyAdjProject/Oldroyd_B_stabilizaiton/couette /Users/jaekwangkim/codes/MyAdjProject/Oldroyd_B_stabilizaiton/couette /Users/jaekwangkim/codes/MyAdjProject/Oldroyd_B_stabilizaiton/couette /Users/jaekwangkim/codes/MyAdjProject/Oldroyd_B_stabilizaiton/couette /Users/jaekwangkim/codes/MyAdjProject/Oldroyd_B_stabilizaiton/couette/CMakeFiles/runclean.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/runclean.dir/depend

