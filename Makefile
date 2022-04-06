mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))

$(info $$current_dir is [${current_dir}])

new:
	g++ -c -fPIC -O3 integrand.cpp 
	g++ -shared -o integrand.so integrand.o
	g++ -c -fPIC -I ${current_dir} integrand_f128.cpp 
	g++ -shared -L ${current_dir}/mppp-0.26 -o integrand_f128.so integrand_f128.o -lmp++

clean:
	rm integrand.o integrand.so integrand_f128.o integrand_f128.so