# Spectral correlation
this repository is dedicated to provide c++ implementation of spectral correlation and spectral coherence for cyclostationary signal analysis.
Using visual studio 2017.
The project is based on c++ Win32 console application. 
the application calculates spectral correlation/coherence using ACP, fast ACP, FAM, Antoni et al Fast SC and 
Borghesani and Antoni fast dirichlet based SC.
Please cite the Fast ACP code as:
Jaafar K. Alsalaet, Fast Averaged Cyclic Periodogram method to compute spectral correlation and coherence,ISA Transactions,
2022, ISSN 0019-0578, https://doi.org/10.1016/j.isatra.2022.01.029.

# Project Dependency
- FFTW library is required to compile and run the project. Depending on your buil configuration CPU type, 
The header file, dll and intermediate library are already contained in the project.
Copy the file libfftw3-3.dll from projectdir/SigProcess/fftw32 for 32-bit building or projectdir/SigProcess/fftw64
for 64-bit building and put it in your application .exe folder

Current build configuration uses x86.
Once, you extract the downloaded zip file, open the project and select x86 build setting.
The latest c++ toolset is required to build the project.


# License
Please read https://github.com/jaafarkh/spectcorr/blob/add-license-1/LICENSE
