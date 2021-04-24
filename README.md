```
🇵​​​​​🇮​​​​​🇨​​​​​🇮​​​​​-🇱​​​​​🇮​​​​​🇬​​​​​🇬​​​​​🇬​​​​​🇭​​​​​🇹​​​​​🇸​​​​​
```

# UoB Positron Imaging Centre's Improved LIGGGHTS Distribution
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

There are many ways to build the LIGGGHTS package (including an executable, a shared library and a static archive). The easiest method uses CMake, which flexibly handles the process of finding the external libraries that LIGGGHTS depends on. For this, you will need:

- A relatively modern C++ compiler (Clang and GCC are great).
- CMake (at least version 2.8, but the latest one is recommended).
- MPI
- (Recommended) VTK with MPI support

Once you have those programs installed, we can create an "out-of-source" build of LIGGGHTS - that is, we compile the many files of LIGGGHTS in a separate folder, so that our source files are not clutttered. Create a `build` folder at the root of the repository:

```
$LIGGGHTS> mkdir build && cd build
$LIGGGHTS/build>
```

While inside the `build` folder, run `cmake` with the path to the LIGGGHTS source files:

```
$LIGGGHTS/build> cmake ../src/ -DCMAKE_BUILD_TYPE=Release
```

This will generate the required compilation scripts (Makefiles) to create the LIGGGHTS files. It should work straight out of the box, but if your VTK is installed in a non-default location, you need to point CMake to it using the flag `-DVTK_DIR=<filepath>`, where `<filepath>` should be substituted with the absolute path to the VTK installation. The script would then become:

```
$LIGGGHTS/build> cmake ../src/ -DCMAKE_BUILD_TYPE=Release -DVTK_DIR=/usr/local/VTK
```

Now that the right Makefiles were generated, we can (finally) compile LIGGGHTS! We can use multiple cores for the compilation, using the `-j<number of cores>` flag:

```
$LIGGGHTS/build> make -j4
```

Three files should have been created:
- "liggghts": the LIGGGHTS executable. This can run simulation scripts.
- "libliggghts.so": the LIGGGHTS shared library. This can be used to run simulation from Python.
- "libliggghts.a": the LIGGGHTS static library. This can be used to statically link against LIGGGHTS from another C++ program. Don't.


## Install Python Interface

In order to use the Python interface, we need the LIGGGHTS shared library (`libliggghts.so`) to be in a standard location - on Unix systems, that is `/usr/local/lib` - and we need to install the Python package that will communicate with it. Thankfully, both of those are done by the Python script at `python/install.py` - run it with your Python 3 interpreter:

```
$LIGGGHTS> cd python
$LIGGGHTS/python> python3 install.py
```

You might need to give the script `sudo` permissions to copy the `libliggghts.so` file to the system library folder. Just add `sudo` before running the Python script:

```
$LIGGGHTS/python> sudo python3 install.py
```


## JKR model Implementation

The PICI-LIGGGGHTS fork comes with the implementation of the JKR - model. The model was initially implemented in a LIGGGHTS fork from [Tobias Eidevåg](https://github.com/eidevag/LIGGGHTS-PUBLIC-JKR). If you use this model please be kind and cite his [paper](https://doi.org/10.1016/j.powtec.2019.10.085).

An example use and desciption can be found in the [examples folder](examples/PICI-LIGGGHTS/).

