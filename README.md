GPM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation (https://www.gnu.org/licenses/gpl.html), either version 3 of the License, or (at your option) any later version.
It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
A copy of the license is provided in COPYING.

Copyright 2018 Onur Kerimoglu (kerimoglu.o@gmail.com)

# A few words on GPM, the General Plankton Model

'Generality' of the model is about a generic module (plankton.F90), that describes both autotrophic and heterotrophic processes, and that can optionally resolve C, N, P and Si cycles. This module can be coupled to a BGC model multiple times on run time (i.e., without requiring re-compilation) through a .yaml configuration file (see the example yaml file), and at each instance, the autothrophy/heterotrophy fraction of the plankton can be prescribed (using the parameter 'fracaut'). Besides the plankton module, 2 additional modules to describe the recycling of C, N, P and Si in the pelagic (abio_pel.F90) and sediment (abio_sed.F90) are provided.

A full description of the model and its 3-D application to the southern North Sea is provided by Kerimoglu et al., 2020, Biogoesciences (doi:https://doi.org/10.5194/bg-2020-1). The configuration (yaml) file used for this study is provided in the testcases/ directory.

# Obtaining the code and building 

First obtain the FABM source code:

    git clone git://git.code.sf.net/p/fabm/code <FABMDIR>

Obtain the GPM code for FABM:

    git clone https://github.com/OnurKerimoglu/FABM-GPM.git <GPMDIR>


For letting FABM know about GPM, edit $FABMDIR/src/CMakeLists.txt and append 'gpm          # general plankton model' at the end of the DEFAULT_INSTITUTES list. 

FABM and GPM use object-oriented Fortran and therefore require a recent Fortran compiler and a platform-independent build system based on cmake.  For the requirements and installation instructions, please check the FABM-wiki at https://github.com/fabm-model/fabm/wiki/Building-and-installing. 
