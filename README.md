# Installation:

0. Install Dependencies
   a. CGAL:
```bash
conda install -c conda-forge cgal
```
    b. CMake
```bash
conda install -c anaconda cmake
```

1. Clone SDOT and dependenct repositories:
```bash
git clone --recurse-submodules git@public.git.erdc.dren.mil:sirc/sdot.git
```

2. Run CMake to configure SDOT:
```bash
cd <SDOT/source/dir>
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<Some/Install/Dir> ..
```

3. Build SDOT:
```bash
make -j4 install
```

4. Add the SDOT paths to your environment variables 
```bash
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH=<Some/Install/Dir>/lib
export PYTHONPATH=$PYTHONPATH:<Some/Install/Dir>/lib:<Some/Install/Dir>/lib/python/
```
    
# Examples:
Try running the random points in a rectangle example:
```bash
cd SDOT/examples/python
python RandomRectangle.py
```