# LUDY

**Contents**

The LUDY folder consists of the following core files
- integrand.cpp: backend f64 evaluation of the integrand
- integrand_f128.cpp: backend f128 evaluation of the integrand
- integrand_wrapper.py: wrapper for python evaluation of C++ functions contained in integrand.cpp and integrand_f128.cpp
- NLO_integrands.py: module containing the class NLO_integrands, whose methods allow for the safe evaluation of the C++ function either in momentum space or in x space, and provide the basic infrastructure for binning

  The entire code can be compiled by running the Makefile. Other than the files given above and the Makefile, the folder also includes:
  - example.py: self-explanatory example of code
  - my_code.py: my own sandbox
  - dual: folder containing the implementation of dual numbers

**Requirements**

Other than standard libraries, the C++ code requires the following
- tgmath.h
- complex.h
- quadmath.h
- stdlib.h

**Methods**

Please check the self-sufficient example for a comprehensive look at the methods of the NLO_integrands class and their workings.

