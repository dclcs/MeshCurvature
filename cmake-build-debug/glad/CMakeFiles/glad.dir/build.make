# CMAKE generated file: DO NOT EDIT!
# Generated by "NMake Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE
NULL=nul
!ENDIF
SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = M:\Paper\new\CLion-2020.2.win\bin\cmake\win\bin\cmake.exe

# The command to remove a file.
RM = M:\Paper\new\CLion-2020.2.win\bin\cmake\win\bin\cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = E:\LIBIGL\MeshCurvature

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = E:\LIBIGL\MeshCurvature\cmake-build-debug

# Include any dependencies generated for this target.
include glad\CMakeFiles\glad.dir\depend.make

# Include the progress variables for this target.
include glad\CMakeFiles\glad.dir\progress.make

# Include the compile flags for this target's objects.
include glad\CMakeFiles\glad.dir\flags.make

glad\CMakeFiles\glad.dir\src\glad.c.obj: glad\CMakeFiles\glad.dir\flags.make
glad\CMakeFiles\glad.dir\src\glad.c.obj: E:\LIBIGL\libigl\external\glad\src\glad.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\LIBIGL\MeshCurvature\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object glad/CMakeFiles/glad.dir/src/glad.c.obj"
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug\glad
	D:\vs\VC\Tools\MSVC\14.27.29110\bin\Hostx86\x86\cl.exe @<<
 /nologo $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) /FoCMakeFiles\glad.dir\src\glad.c.obj /FdCMakeFiles\glad.dir\glad.pdb /FS -c E:\LIBIGL\libigl\external\glad\src\glad.c
<<
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug

glad\CMakeFiles\glad.dir\src\glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glad.dir/src/glad.c.i"
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug\glad
	D:\vs\VC\Tools\MSVC\14.27.29110\bin\Hostx86\x86\cl.exe > CMakeFiles\glad.dir\src\glad.c.i @<<
 /nologo $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E E:\LIBIGL\libigl\external\glad\src\glad.c
<<
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug

glad\CMakeFiles\glad.dir\src\glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glad.dir/src/glad.c.s"
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug\glad
	D:\vs\VC\Tools\MSVC\14.27.29110\bin\Hostx86\x86\cl.exe @<<
 /nologo $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) /FoNUL /FAs /FaCMakeFiles\glad.dir\src\glad.c.s /c E:\LIBIGL\libigl\external\glad\src\glad.c
<<
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug

# Object files for target glad
glad_OBJECTS = \
"CMakeFiles\glad.dir\src\glad.c.obj"

# External object files for target glad
glad_EXTERNAL_OBJECTS =

glad\glad.lib: glad\CMakeFiles\glad.dir\src\glad.c.obj
glad\glad.lib: glad\CMakeFiles\glad.dir\build.make
glad\glad.lib: glad\CMakeFiles\glad.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=E:\LIBIGL\MeshCurvature\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library glad.lib"
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug\glad
	$(CMAKE_COMMAND) -P CMakeFiles\glad.dir\cmake_clean_target.cmake
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug\glad
	D:\vs\VC\Tools\MSVC\14.27.29110\bin\Hostx86\x86\link.exe /lib /nologo /machine:X86 /out:glad.lib @CMakeFiles\glad.dir\objects1.rsp 
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug

# Rule to build all files generated by this target.
glad\CMakeFiles\glad.dir\build: glad\glad.lib

.PHONY : glad\CMakeFiles\glad.dir\build

glad\CMakeFiles\glad.dir\clean:
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug\glad
	$(CMAKE_COMMAND) -P CMakeFiles\glad.dir\cmake_clean.cmake
	cd E:\LIBIGL\MeshCurvature\cmake-build-debug
.PHONY : glad\CMakeFiles\glad.dir\clean

glad\CMakeFiles\glad.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" E:\LIBIGL\MeshCurvature E:\LIBIGL\libigl\external\glad E:\LIBIGL\MeshCurvature\cmake-build-debug E:\LIBIGL\MeshCurvature\cmake-build-debug\glad E:\LIBIGL\MeshCurvature\cmake-build-debug\glad\CMakeFiles\glad.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : glad\CMakeFiles\glad.dir\depend

