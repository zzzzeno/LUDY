new:
	g++ -c -fPIC -O3 integrand.cpp 
	g++ -shared -o integrand.so integrand.o
	g++ -c -fPIC -I/scratch/zcapatti/drell_yan_2 integrand_f128.cpp 
	g++ -shared -L/scratch/zcapatti/drell_yan_2/mppp-0.26 -o integrand_f128.so integrand_f128.o -lmp++
