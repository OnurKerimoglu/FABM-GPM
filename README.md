GPM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation (https://www.gnu.org/licenses/gpl.html), either version 3 of the License, or (at your option) any later version.
It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
A copy of the license is provided in COPYING.

Copyright 2018 Onur Kerimoglu (kerimoglu.o@gmail.com)

# A few words on GPM, the General Plankton Model

'Generality' of the model is about a generic module (plankton.F90), that describes both autotrophic and heterotrophic processes, and that can optionally resolve C, N, P and Si cycles. This module can be coupled to a BGC model multiple times on run time (i.e., without requiring re-compilation) through a .yaml configuration file (see the example yaml files in testcases folder), and at each instance, the autothrophy/heterotrophy fraction of the plankton can be prescribed (using the parameter 'fracaut'). Besides the plankton module, 2 additional modules to describe the recycling of C, N, P and Si in the pelagic (abio_pel.F90) and sediment (abio_sed.F90) are provided.

An earlier version of the model based on C- and P- only, and its 1-D application to Lake Bourget (France) was described by Kerimoglu et al., 2017 Ecological Modelling 359, 415–433

# Obtaining the code and building 

First obtain the FABM source code:

    git clone git://git.code.sf.net/p/fabm/code <FABMDIR>

(Replace `<FABMDIR>` with the directory with the FABM code, e.g., ~/fabm-git.)

Obtain the GPM code for FABM:

    git clone git@bitbucket.org:OnurKerimoglu/gpm.git <GPMDIR>

(Replace `<GPMDIR>` with the directory with the GPM code, e.g., ~/gpm-git.)

For letting FABM know about GPM, edit $FABMDIR/src/CMakeLists.txt and append 'gpm          # general plankton model' at the end of the DEFAULT_INSTITUTES list. 

FABM and GPM use object-oriented Fortran and therefore require a recent Fortran compiler, such as Intel Fortran 12.1 or higher and gfortran 4.7 or higher.

FABM and GPM use a platform-independent build system based on [cmake](http://www.cmake.org). You'll need version 2.8.8 or higher. First check whether you have that installed: execute `cmake --version` on the command line.

## GOTM + FABM + GPM

First obtain the latest (developers') version of the GOTM code from its git repository:

    git clone https://github.com/gotm-model/code.git gotm-git

To build GOTM with FABM support, create a build directory, call cmake to generate makefiles, and make to compile and install. For instance:

    mkdir -p ~/build/gotm && cd ~/build/gotm
    cmake <GOTMDIR>/src -DFABM_BASE=<FABMDIR> -DFABM_GPM_BASE=<GPMDIR>
    make install

In the above, replace `<GOTMDIR>` with the directory with the GOTM source code, e.g., ~/gotm-git if you executed `git clone` in you home directory. Also, replace `<FABMDIR>` with the directory with the FABM code, e.g., ~/fabm-git and `<GPMDIR>` with the directory with the GPM code, e.g., ~/gpm-git.

Now you should have a GOTM executable with FABM and GPM support at `~/local/gotm/bin/gotm`.

It is good practice to keep up to date with the latest code from the GPM, FABM and GOTM repositories by regularly executing `git pull` in a directory of each repository.

If either the GPM, FABM or GOTM source codes change (e.g., because changes you made to the code yourself, or after `git pull`), you will need to recompile. This does NOT require rerunning cmake. Instead, you need to return to the build directory and rerun `make install`. For instance `cd ~/build/gotm && make install`.

## FABM-GPM 0d

The 0d driver allows you to run FABM models in a "well-mixed box", under arbitrary (time-varying) environmental forcing.

To build the 0d driver, you need to create a directory to build the code in, call `cmake` to generate makefiles, and call `make` to compile and install the FABM library. Usually, the following suffices for this:

    mkdir -p ~/build/fabm-0d && cd ~/build/fabm-0d
    cmake <FABMDIR>/src/drivers/0d -DGOTM_BASE=<GOTMDIR> -DFABM_GPM_BASE=<GPMDIR>
    make install

In the above, replace `<GPMDIR>` with the path to directory with the GPM source code (e.g., ~/ersem-git), `<FABMDIR>` with the path to directory with the FABM source code (e.g., ~/fabm-git), and `<GOTMDIR>` with the path to directory with the GOTM source code (e.g., ~/gotm-git). The latter is needed because the 0d driver uses GOTM routines for input, output, time integration, etc. If you experience issues related to NetCDF, see [tips and tricks/troubleshooting](#tips-and-tricks-troubleshooting).

This will give you an executable at `~/local/fabm/0d/bin/fabm0d`.
