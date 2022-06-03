import math
from random import random
import numpy as np
from scipy.spatial.transform import Rotation as Roto
import vegas
import yaml
from NLO_integrands import NLO_integrands, LO_integrands
import time
import matplotlib.pyplot as plt
from integrand_wrapper import cut_res, test_function_py
from integrand_wrapper import set_r, set_kin, set_sig, set_defo, set_MUV

if __name__ == "__main__":

    #loop integration basis
    basis=np.linalg.inv(np.array([[1,1,0],[0,1,0],[0,0,1]])) 
    #center of mass energy: used to bound Bjorken variables and normalise stuff
    eCM=2
    #a=x_1-x_2: fixing a fixes the difference between the Bjorken variables
    a=0.2
    #min_pt: minimum transverse momentum of final state jet
    min_pt=0.
    #max_pt: maximum transverse momentum of final state jet
    max_pt=6500.
    #n_bind: number of bins for the pt of the transverse momentum of final state jet
    n_bins=20
    #mZ: mass of final state boson
    mZ=1
    #res_c: collinear resolution of initial-state jets
    res_c=0.1
    #res_s: soft resolution of initial-state jets
    res_s=0.1
    #sigma: standard deviation of h(t) gaussian
    sigma=1
    #tag: supergraph of choice
    my_tag='st'
    #mode: determines the type of output. 
    #Default is "min_pt", should not be changed.
    my_mode='min_pt'
    #phase: allows to select real or imag part of the result
    my_phase='real'



    #INITIALISE NLO_integrands
    integrand=NLO_integrands(mZ, eCM, res_c, res_s, sigma, a, min_pt, max_pt, n_bins, basis, tag=my_tag, debug=False, mode=my_mode,phase=my_phase)


    ####################################################################
    #EVALUATE INTEGRAND
    #The output is a vector of structures. The single structure represents the evaluation of a single interference diagrams, 
    #and contains the following information
    #class cut_res:
    #   def __init__(self, res, jac, j1, spin1, j2, spin2, pg):
    #       self.res = res
    #       self.jac = jac
    #       self.j1 = j1
    #       self.spin1 = spin1
    #       self.j2 = j2
    #       self.spin2 = spin2
    #       self.pg = pg
    #
    # self.res is the full evaluation of the integrans, including jacobian 
    # self.jac is the jacobian, including that stemming from causal flow, Lorentz flow and various Dirac deltas
    # self.j1 momentum of the first jet
    # self.j2 momentum of the second jet
    # self.spin1 spin of the first jet
    # self.spin2 spin of the second jet
    # self.pg momentum of the final state jet
    ####################################################################


    ######################################################################
    #It is possible to evaluate the integrand directly in momentum space, 
    #by providing a set [[kx,ky,kz],[lx,ly,lz],[qx,qy,qz]] of vectors,
    #and obtain an output computed either with f64 or f128 precision
    ######################################################################

    # moms=[[1.0,1.2,0.17],[-0.56,0.33,0.17],[0.77,0.98,-1.54]]
    # f64_res=integrand.sample_f64_mom(moms)
    # f128_res=integrand.sample_f128_mom(moms)

    # for r in f64_res:
    #     print("f64 result")
    #     print(r.res)
    #     print(r.jac)
    #     print(r.j1)
    #     print(r.spin1)
    #     print(r.j2)
    #     print(r.spin2)
    #     print(r.pg)

    # for r in f128_res:
    #     print("f128 result")
    #     print(r.res)
    #     print(r.jac)
    #     print(r.j1)
    #     print(r.spin1)
    #     print(r.j2)
    #     print(r.spin2)
    #     print(r.pg)

    # ##########################################################################
    # #It is also possible to evaluate with input given in x space, and result being subject to a safety mechanism
    # ##########################################################################

    # xs=[0.1,0.2,0.1,0.17,0.44,0.11,0.17,0.6,0.21]

    # safe_res=integrand.safe_eval(xs)

    # for r in safe_res:
    #     print(r.res)
    #     print(r.jac)
    #     print(r.j1)
    #     print(r.spin1)
    #     print(r.j2)
    #     print(r.spin2)
    #     print(r.pg)

    ##########################################################################
    #Finally, there's a function that directly bins the distribution in the transverse momentum of the final-state jet.
    #The output is a vector of values attached to each bin
    ##########################################################################

    # xs=[0.1,0.2,0.5,0.17,0.14,0.11,0.17,0.6,0.21]

    # safe_res=integrand.safe_eval_bin(xs)

    # print(safe_res)






    ##########################################################################
    #Let's test this with some benchmarked number
    ##########################################################################

    #number of learning iterations for vegas
    n_iterations_learning=10 
    #size of learning iterations for vegas
    n_points_iteration_learning=100000
    #number of refining iterations for vegas
    n_iterations=10 
    #size of refining iterations for vegas
    n_points_iteration=1000000


    #################################################
    #We start with the "LO" cross-section
    #################################################

    # benchmark=1.3820603e-03

    # integrand0=LO_integrands(eCM=13000, m=91.1188, x1=0.595, n_bins=0)

    # m=9.118800e+01
    # s=13000**2

    # def f_LO(x):
    #     #integrand0.set_x1(x[6])
    #     return integrand0.LO_bin(x)
        
    # integ0 = vegas.Integrator([[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]])#, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])
    # integ0(f_LO, nitn=10, neval=100000)

    # result = integ0(f_LO, nitn=10, neval=1000000)
    # if integ0.mpi_rank == 0:
    #     print(result.summary())
    #     print('result = %s    Q = %.2f' % (result, result.Q))



    integrand0=LO_integrands(eCM=1, m=1, x1=0.595, n_bins=0)

    m=1#9.118800e+01

    set_kin(m,m)
    integrand0.set_a_jj(0)

    def f_LO(x):
        #integrand0.set_x1(x[6])
        return integrand0.LO_bin(x)
        
    # integ0 = vegas.Integrator([[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]])#, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])
    # integ0(f_LO, nitn=10, neval=100000)

    # result = integ0(f_LO, nitn=10, neval=1000000)
    # if integ0.mpi_rank == 0:
    #     print(result.summary())
    #     print('result = %s    Q = %.2f' % (result, result.Q))



    #################################################
    #and continue with the "NLO" cross-sections
    #################################################

    basis=np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))


    #Set kinematics and hyperparameters at the C++ level
    set_kin(91.188,13000)
    set_r(0.,0.)
    set_sig(2.)

    #Benchmark numbers to reproduce for given hyperparameters
    benchmarks={"st":-0.01199844,
                "tqg":0.13399,
                "s":0.0068855,
                "tqq":2*0.2875,
                "u":0.20195
        }

    my_tag="all_trees"

    my_integrand=NLO_integrands(91.188,13000,0.,0.,2,0.,10,2000,1,basis,tag=my_tag, obs_tag='JADE', obs_lambda=1, phase='real')
    my_integrand.set_digits(8,4)

    def f(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand.safe_eval_bin(xs)

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        return res

    integ = vegas.Integrator([[0.000001,0.999999], [0,1], [0,1], [0.000001,0.999999], [0,1], [0,1], [0.000001,0.999999], [0,1], [0,1],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

    #integ(f, nitn=n_iterations_learning, neval=n_points_iteration_learning, analyzer=vegas.reporter())
    integ(f, nitn=n_iterations_learning, neval=n_points_iteration_learning)
    result = integ(f, nitn=n_iterations, neval=n_points_iteration)
    if integ.mpi_rank == 0:
        print(result.summary())
    if integ.mpi_rank == 0:
        print('result = %s    Q = %.2f' % (result[0], result.Q))
        for biny in result:
            print(biny)

    print(my_integrand.n_unstable_f64)
    print(my_integrand.n_unstable_f128)


    #################################################
    #and continue with the "NLO" cross-sections
    #################################################

    # basis=np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))


    # #Set kinematics and hyperparameters at the C++ level
    # # set_kin(91.188,13000)
    # # set_r(0.,0.)
    # # set_sig(2.)


    # set_kin(91.188,13000)
    # set_r(0.1,0.1)
    # set_sig(2.)
    # set_MUV(1000)
    # set_defo(0.5,0.1,1)
    # basispp=np.linalg.inv(np.array([[1,-1,0],[0,1,0],[0,0,1]]))

    # #Benchmark numbers to reproduce for given hyperparameters
    # my_tag="dt"

    # my_integrand=NLO_integrands(91.188,13000,0.1,0.1,2,0.,0,10000,1,basispp,tag=my_tag, obs_tag='JADE', obs_lambda=1, phase='real')
    # my_integrand.set_digits(8,4)

    # def f(x):

    #     ap=10*math.tan(math.pi*(x[9]-1/2))
    #     my_integrand.set_a(ap)

    #     xs=[x[i] for i in range(0,9)]

    #     bins=my_integrand.safe_eval_bin(xs)

    #     res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

    #     #print("python res")
    #     #print(res)

    #     return res

    # integ = vegas.Integrator([[0.000001,0.999999], [0,1], [0,1], [0.000001,0.999999], [0,1], [0,1], [0.000001,0.999999], [0,1], [0,1],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

    # #integ(f, nitn=n_iterations_learning, neval=n_points_iteration_learning, analyzer=vegas.reporter())
    # integ(f, nitn=n_iterations_learning, neval=n_points_iteration_learning)
    # result = integ(f, nitn=n_iterations, neval=n_points_iteration)
    # if integ.mpi_rank == 0:
    #     print(result.summary())
    # if integ.mpi_rank == 0:
    #     print('result = %s    Q = %.2f' % (result[0], result.Q))
    #     for biny in result:
    #         print(biny)

    #print(my_integrand.n_unstable_f64)
    #print(my_integrand.n_unstable_f128)




    #basispp=np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))
    basispp=np.linalg.inv(np.array([[1,-1,0],[0,1,0],[0,0,1]]))
    ###IMAG###
    
    

    #set_defo(0.5,0.069,1)
    #eCM=13000
    #mZ=91.188
    #eCM=91.188
    #mZ=91.188
    mZ=1
    eCM=1
    lamb=1.0#eCM
    mij=eCM*0.069

    #TEST HYPERPARAMETER INDEPENDENCEEEEE
    lamb_bc=1.0

    set_MUV(mZ)

    set_defo(lamb,mij,lamb_bc)

    set_sig(2)
    set_r(0.0,0.0)
    set_kin(mZ,eCM)
    

    n_iterations_learning=10
    n_points_iteration_learning=100000
    n_iterations=10
    n_points_iteration=2000000


    for i in range(0,9):

        lambd=0.3/10.+i*0.3/10.

        my_integrand_real=NLO_integrands(mZ,eCM,0.05,0.05,2,0.,0,2000,1,basispp,tag="se",obs_tag='HEMI', obs_lambda=lambd,phase='real')
        my_integrand_real.set_digits(8,10)

    


        def freal(x):

            # ap=10*math.tan(math.pi*(x[9]-1/2))
            # my_integrand_real.set_a(ap)

            xs=[x[i] for i in range(0,9)]

            #WE WILL USE SAFE EVAL FOR NO CONSTRAINT ON x1 and x2
            #bins=my_integrand_real.safe_eval_bin(xs)

            bins=my_integrand_real.safe_eval(xs)
            
            #res=[biny for biny in bins]

            res=0
            for r in bins:
                res=res+r.res*r.jac

            return res#*10*math.pi/pow(math.sin(math.pi*x[9]),2)


        
        integ = vegas.Integrator([[0.0001,1.0],[0.0,1.0],[0.0,1.0],[0.0001,1.0],[0.0,1.0],[0.0,1.0],[0.0001,1.0],[0.0,1.0],[0.0,1.0],[0.0,1.0]])

        integ(freal, nitn=n_iterations_learning, neval=n_points_iteration_learning)
        result = integ(freal, nitn=n_iterations, neval=n_points_iteration)
        if integ.mpi_rank == 0:
            print("#######RUN WITH LAMBDA########")
            print(lambd)
            # print(result.summary())
        if integ.mpi_rank == 0:
            print('result = %s    Q = %.2f' % (result, result.Q))
            # print(my_integrand_real.n_unstable_f64)
            # print(my_integrand_real.n_unstable_f128)
            # print("-----")
            # for i in range(0,min(10,len(my_integrand_real.unstablef64))):
            #     print(my_integrand_real.unstablef64[i][0])
            #     print(my_integrand_real.unstablef64[i][1])
            # print("-----")
            # for i in range(0,min(10,len(my_integrand_real.unstablef128))):
            #     print(my_integrand_real.unstablef128[i][0])
                # print(my_integrand_real.unstablef128[i][1])

    # print("----max----")
    # print(my_integrand_real.max_val)
    # print(my_integrand_real.max_x)
    # print(my_integrand_real.max_k)


    # my_x=[0.86119181, 0.06631761, 0.40293758, 0.48049277, 0.37842339, 0.60848394, 0.73094879, 0.94802582, 0.42056274]
    my_x2=[0.9929777921706654, 0.9034888023489799, 0.26481297075893206, 0.8592341102346798, 0.2552656976761265, 0.22946627051380367, 0.8990129029914548, 0.2792449532812115, 0.2846998961319826,0.5]







    # # print(sum(r.res*r.jac for r in my_integrand_real.safe_eval(my_x)))

    # # listy=[sum(r.res*r.jac for r in my_integrand_real.safe_eval([0.86119181, 0.06631761+pow(10,-i), 0.40293758, 0.48049277+pow(10,-i), 0.37842339, 0.60848394+pow(10,-i), 0.73094879, 0.94802582, 0.42056274])) for i in range(0,10)]


    # # print(listy)

    # # print(my_integrand_real.x_param(my_x))

    # print("point number 1")
    # print(freal(my_x2))

    # print("point number 2")
    # print(freal(my_k3))

    #print([freal([1-pow(10,-i), 0.9034888023489799, 0.26481297075893206, 0.8592341102346798, 0.2552656976761265, 0.22946627051380367, 0.8990129029914548, 0.2792449532812115, 0.2846998961319826,0.5]) for i in range(1,6)])