open-quantum-systems
====================

(C) Dr Agnieszka M. Werpachowska 2011-2017

Licensed under GNU Public License v3 (see https://www.gnu.org/licenses/gpl.html and LICENSE.txt). Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* Neither the name of Averisera Ltd nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

This is a C++ code for modelling open quantum systems, accompanying my paper 

A. M. Werpachowska, "Reduced operator approximation for modelling open quantum systems", Open Syst. Inf. Dyn., 22, 1550008 (2015)  http://www.worldscientific.com/doi/abs/10.1142/S1230161215500080 (http://arxiv.org/abs/1212.1753)

It includes the Reduced Operator Approximation, Quantum State Diffusion and Pseudomode methods, algorithms for modelling quantum systems with static noise, as well as some general numerical methods.

Requirements:
- C++ 11 compiler (I recommend g++ 4.8.x)
- Eigen 3.x library
- GNU Scientific Library library
- FFTW 3.x library
- SCons build system
- google-test library (for unit tests)
