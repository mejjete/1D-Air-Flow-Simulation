# One dimensional air flow simulator for an underground mine

## Compile and run

1. Download and install additional dependencies:
- OpenMPI 5.0.5
- OpenMP (either libomp for clang or libgomp for gcc)
- Boost 1.66.0

2. Execute following commands from within the root directory of the project:
Create build directory
    ```
    $  mkdir build
    $  cd build
    ```
    Configure and build a project

    ```
    $  cmake -DCMAKE_CXX_COMPILER=mpicxx ..
    $  cmake --build .
    ```

    Additionaly, every executable instance supports output dump to produce a visualization reports that can later be used by tools like `dot` to do simple plotting. To enable this option, just pass `-DENABLE_DUMP=On` during cmake configuration.
    
    Example:
    
    ```
    $  cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=mpicxx -DENABLE_DUMP=On ..
    ```

3. Execute
Upon successfull configuration and build, it should generate `bin/` directory inside `build/` folder. There should be a couple of executables ready to run either in single or parallel environment alogside with configurational json files.

## Visualization
## Configuration file structure
## Sequential single-edge simulator 
## Sequential network simulator

