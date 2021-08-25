![Build Status](https://github.com/mparno/sdot2d/actions/workflows/build-test.yml/badge.svg)

# Installation:

0. Install Dependencies:

```bash
conda install -c conda-forge cgal cmake
```

1. Clone sdot2d and dependent repositories

```bash
git clone --recurse-submodules git@github.com:mparno/sdot2d.git
```

2. Run CMake to configure sdot2d

```bash
cd sdot2d
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<Some/Install/Dir> ..
```
where you should replace `<Some/Install/Dir>` with the location where you want to install the SDOT2d library, for example `~/Installations/SDOT_INSTALL`.

You can also tell CMake where to look for GTest using the following command instead:
```bash
cmake -DCMAKE_INSTALL_PREFIX=<Some/Install/Dir> -DSDOT_GTEST_DIR=<gtest/install/dir> ..
```
The `gtest/install/dir` should be the folder that contains the `include` and
`lib` folders where GTest has been installed.

3. Build SDOT (including optional tests)

```bash
make -j4 install
```

4. (Optional) Run the tests
```bash
make test
```
or, to run just the unit tests
```bash
./RunUnitTests
```

5. Add the SDOT paths to your environment variables:

#### OSX:
```bash
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:<Some/Install/Dir>/lib
export PYTHONPATH=$PYTHONPATH:<Some/Install/Dir>/lib:<Some/Install/Dir>/lib/python/
```
#### Linux:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<Some/Install/Dir>/lib
export PYTHONPATH=$PYTHONPATH:<Some/Install/Dir>/lib:<Some/Install/Dir>/lib/python/
```

# Examples:
Try running the random points in a rectangle example

```bash
cd sdot2d/examples/python
python RandomRectangle.py
```
