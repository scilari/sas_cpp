# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.12

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2018.2.5\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2018.2.5\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\ilari\CLionProjects\sas_cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug

# Include any dependencies generated for this target.
include test/CMakeFiles/allTests.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/allTests.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/allTests.dir/flags.make

test/CMakeFiles/allTests.dir/sas_test.cpp.obj: test/CMakeFiles/allTests.dir/flags.make
test/CMakeFiles/allTests.dir/sas_test.cpp.obj: test/CMakeFiles/allTests.dir/includes_CXX.rsp
test/CMakeFiles/allTests.dir/sas_test.cpp.obj: ../test/sas_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/allTests.dir/sas_test.cpp.obj"
	cd /d C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug\test && C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\allTests.dir\sas_test.cpp.obj -c C:\Users\ilari\CLionProjects\sas_cpp\test\sas_test.cpp

test/CMakeFiles/allTests.dir/sas_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/allTests.dir/sas_test.cpp.i"
	cd /d C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug\test && C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\ilari\CLionProjects\sas_cpp\test\sas_test.cpp > CMakeFiles\allTests.dir\sas_test.cpp.i

test/CMakeFiles/allTests.dir/sas_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/allTests.dir/sas_test.cpp.s"
	cd /d C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug\test && C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\ilari\CLionProjects\sas_cpp\test\sas_test.cpp -o CMakeFiles\allTests.dir\sas_test.cpp.s

# Object files for target allTests
allTests_OBJECTS = \
"CMakeFiles/allTests.dir/sas_test.cpp.obj"

# External object files for target allTests
allTests_EXTERNAL_OBJECTS =

test/allTests.exe: test/CMakeFiles/allTests.dir/sas_test.cpp.obj
test/allTests.exe: test/CMakeFiles/allTests.dir/build.make
test/allTests.exe: src/libsas.a
test/allTests.exe: lib/libgtest_maind.a
test/allTests.exe: lib/libgtestd.a
test/allTests.exe: test/CMakeFiles/allTests.dir/linklibs.rsp
test/allTests.exe: test/CMakeFiles/allTests.dir/objects1.rsp
test/allTests.exe: test/CMakeFiles/allTests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable allTests.exe"
	cd /d C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug\test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\allTests.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/allTests.dir/build: test/allTests.exe

.PHONY : test/CMakeFiles/allTests.dir/build

test/CMakeFiles/allTests.dir/clean:
	cd /d C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug\test && $(CMAKE_COMMAND) -P CMakeFiles\allTests.dir\cmake_clean.cmake
.PHONY : test/CMakeFiles/allTests.dir/clean

test/CMakeFiles/allTests.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\ilari\CLionProjects\sas_cpp C:\Users\ilari\CLionProjects\sas_cpp\test C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug\test C:\Users\ilari\CLionProjects\sas_cpp\cmake-build-debug\test\CMakeFiles\allTests.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/allTests.dir/depend

