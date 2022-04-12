mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))

ifeq ($(GMP_PATH),)
GMP_PATH=/Users/vjhirsch/HEP_programs/gmp-6.2.1_build
endif

ifeq ($(MPP_PATH),)
MPP_PATH=/Users/vjhirsch/MG5/3.0.2.py3/PLUGIN/tderivative_alphaloop/libraries/form/mppp-0.26
endif

ifeq ($(MPFR_PATH),)
MPFR_PATH=/Users/vjhirsch/HEP_programs/mpfr-4.1.0_build
endif

ifeq ($(GCC),)
GCC=g++
endif

ifeq ($(OPTLEVEL),)
OPTLEVEL=3
endif

all: new

integrand.o:
	$(GCC) -c -fPIC -O$(OPTLEVEL) integrand.cpp

integrand.so: integrand.o
	$(GCC) -shared -o integrand.so integrand.o

integrand_f128.o:
	$(GCC) -c -fPIC -O$(OPTLEVEL) -I $(GMP_PATH)/include -I $(MPP_PATH)/include -I $(MPFR_PATH)/include integrand_f128.cpp

integrand_f128.so: integrand_f128.o
	$(GCC) -shared -O$(OPTLEVEL) -L$(MPP_PATH) -L$(GMP_PATH)/lib -L$(MPFR_PATH)/lib -o integrand_f128.so integrand_f128.o -lmp++

new: integrand.so integrand_f128.so

clean:
	rm -f integrand.o integrand.so integrand_f128.o integrand_f128.so
