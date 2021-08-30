# Spectral correlation
this repository is dedicated to provide c++ implement spectral correlation and spectral coherence for cyclostationary signal analysis.
Using visual studio 2017
The project is based on c++ Win32 console application 
Project Dependency
- FFTW library is required to compile and run the project. Depending on your buil configuration CPU type, 
rename either libfftw3-3-x86.dll or libfftw3-3-x64.dll to libfftw3-3.dll and put it in your application .exe folder

Current build configuration uses x86.
The latest c++ toolset is required to build the project
the application calculates spectral correlation/coherence using ACP, fast ACP, FAM, Antoni etal Fast SC and 
Borghesani and antoni fast dirichlet based SC

# License
Please read https://github.com/jaafarkh/spectcorr/blob/add-license-1/LICENSE
