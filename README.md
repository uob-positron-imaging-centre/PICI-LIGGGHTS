```
🇵​​​​​🇮​​​​​🇨​​​​​🇮​​​​​-🇱​​​​​🇮​​​​​🇬​​​​​🇬​​​​​🇬​​​​​🇭​​​​​🇹​​​​​🇸​​​​​
```

# UoB Positron Imaging Centre's Improved LIGGGHTS Distribution - PICI-LIGGGGHTS-3.8.1
## Internal code used through the improved Python interface

The LIGGGHTS® distribution includes the following files and directories:

- README          this file
- LICENSE         the GNU General Public License (GPL)
- doc             documentation
- examples        simple example simulation setups
- lib             libraries LIGGGHTS® can be linked with
- python          Python wrapper on LIGGGHTS® as a library
- src             source files


## Install LIGGGHTS

There are multiple ways to build the LIGGGHTS package, including an executable, a shared library and a static archive. The easiest method uses CMake, which flexibly handles the process of finding the external libraries that LIGGGHTS depends on. For this, you will need:

- A relatively modern C++ compiler (Clang and GCC are great).
- CMake (at least version 2.8, but the latest one is recommended).
- MPI.
- (Optionally) VTK with MPI support.

Once you have those programs installed, we can create an "out-of-source" build of LIGGGHTS - that is, we compile the many files of LIGGGHTS in a separate folder, so that our source files are not clutttered. Create a `build` folder at the root of the repository:

```bash
$LIGGGHTS_DIR> mkdir build && cd build
$LIGGGHTS_DIR/build>
```

Where `LIGGGHTS_DIR` is the location where you cloned the repository to. Inside the `build` folder, run `cmake` with the path to the LIGGGHTS source files:

```bash
$LIGGGHTS/build> cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -march=native -fPIC" ../src/
```

Short explanation for the `CMAKE_CXX_FLAGS` that are passed to the compiler: `-O3` allows expensive optimisation passes, `-march=native` tunes the generated code to the specific CPU on your machine (so the resulting binary may not work on other computers), while `-fPIC` allows the code to be loaded as a shared library.

The command above will generate the required compilation scripts (Makefiles) to create the LIGGGHTS files. It should work straight out of the box, but if your VTK is installed in a non-default location, you need to point CMake to it using the flag `-DVTK_DIR=<filepath>`, where `<filepath>` should be substituted with the absolute path to the VTK installation (which should also contain some `.cmake` files). The script would then become:

```bash
$LIGGGHTS/build> cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -march=native -fPIC" -DVTK_DIR=/usr/local/VTK ../src/
```

Now that the right Makefiles were generated, we can (finally) compile LIGGGHTS! We can use multiple cores for the compilation, using the `-j<number of cores>` flag:

```bash
$LIGGGHTS/build> make -j4
```

Three files should have been created:
- "liggghts": the LIGGGHTS executable. This can run simulation scripts.
- "libliggghts.so": the LIGGGHTS shared library. This can be used to run simulation from Python.
- "libliggghts.a": the LIGGGHTS static library. This can be used to statically link against LIGGGHTS from another C++ program. Don't.


## Install Python Interface

In order to use the Python interface, the system needs to know where to find the LIGGGHTS shared library (`libliggghts.so`); you can either:
- Copy it to a standard location like `/usr/local/lib` on Unix-based systems.
- Add the its directory to the `LD_LIBRARY_PATH` environment variable, e.g. `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:LIGGGHTS_DIR/build` (recommended).

Finally, Python needs to know where to find the `LIGGGHTS_DIR/python/liggghts.py` module; again, you can either:
- Copy it to the Python "site-packages".
- Add its directory to the `PYTHONPATH` environment variable, e.g. `export PYTHONPATH=$PYTHONPATH:LIGGGHTS_DIR/python` (recommended).

Where `LIGGGHTS_DIR` is the location where you cloned this repository to. To test whether everything is in place.



## JKR model Implementation

The PICI-LIGGGGHTS fork comes with the implementation of the JKR - model. The model was initially implemented in a LIGGGHTS fork from [Tobias Eidevåg](https://github.com/eidevag/LIGGGHTS-PUBLIC-JKR). If you use this model please be kind and cite his [paper](https://doi.org/10.1016/j.powtec.2019.10.085).

An example use and desciption can be found in the [examples folder](examples/PICI-LIGGGHTS/).

