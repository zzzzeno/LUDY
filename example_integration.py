import math
from random import random
import numpy as np
from scipy.spatial.transform import Rotation as Roto
from integrand_wrapper import cut_res, LO, NLO

import vegas
import yaml

import time
import numpy as np
from random import random

from integrand_wrapper import set_r, set_kin, set_sig, set_defo, set_MUV

import matplotlib.pyplot as plt


def rotator(ks):
    nx=random()
    ny=random()
    nz=random()

    theta=random()*2*math.pi

    axis=theta/(math.sqrt(nx*nx+ny*ny+nz*nz))*np.array([nx,ny,nz])

    r = Roto.from_rotvec(axis)

    rotks=[r.apply(ks[0]),r.apply(ks[1]),r.apply(ks[2])]

    return rotks



class NLO_integrands:

    def __init__(self, m, eCM, res_c, res_s, sigma, a, min_pt, max_pt, n_bins, basis, tag='st', debug=False, mode='min_pt',phase='real'):
        self.eCM = eCM
        self.mZ=m
        self.min_pt = min_pt
        self.max_pt = max_pt
        self.a = a
        self.n_bins = n_bins
        self.basis = basis
        self.debug = debug
        self.tag = tag
        self.mode = mode
        self.phase=phase
        set_kin(m,eCM)
        set_r(res_c,res_s)
        set_sig(sigma)
        self.n_digits_f128=0
        self.n_digits_f64=0
        self.n_unstable_f64=0
        self.n_unstable_f128=0
        self.unstablef64=[]
        self.unstablef128=[]

    def set_digits(self,n_digits_f64,n_digits_f128):
        self.n_digits_f128=n_digits_f128
        self.n_digits_f64=n_digits_f64

    def set_a(self, a):
        self.a = a


    def set_basis(self, basis):
        self.basis = basis


    def x_param(self, x):
        eCMn=1#1300#self.eCM/
        dr1=eCMn*x[0]/(1-x[0])
        dth1=x[1]*2*math.pi
        dph1=x[2]*math.pi
        dr2=eCMn*x[3]/(1-x[3])
        dth2=x[4]*2*math.pi
        dph2=x[5]*math.pi
        dr3=eCMn*x[6]/(1-x[6])
        dth3=x[7]*2*math.pi
        dph3=x[8]*math.pi


        k=np.array([dr1*math.cos(dth1)*math.sin(dph1),dr1*math.sin(dth1)*math.sin(dph1),dr1*math.cos(dph1)])
        l=np.array([dr2*math.cos(dth2)*math.sin(dph2),dr2*math.sin(dth2)*math.sin(dph2),dr2*math.cos(dph2)])
        q=np.array([dr3*math.cos(dth3)*math.sin(dph3),dr3*math.sin(dth3)*math.sin(dph3),dr3*math.cos(dph3)])

        ks=[]

        for b in self.basis:
            kb=b[0]*k+b[1]*l+b[2]*q
            ks.append(kb)

        jacques=pow(eCMn,3)*math.pow(dr1,2)*math.sin(dph1)*math.pow(dr2,2)*math.sin(dph2)*math.pow(dr3,2)*math.sin(dph3)*math.pow(2*math.pi,3)*math.pow(math.pi,3)/(math.pow(1-x[0],2)*math.pow(1-x[3],2)*math.pow(1-x[6],2))

        return [jacques,ks]

    def f64_eval_mom(self, moms):
        res=NLO(moms[0],moms[1],moms[2],self.a,self.tag,'f64',self.phase)
        totres=[]
        for r in res:
            totres.append(r.res*r.jac)#r.res*r.jac)#*r.jac)#r.res*r.jac

        return totres

    def f128_eval_mom(self, moms):

        res=NLO(moms[0],moms[1],moms[2],self.a,self.tag,'f128',self.phase)
        totres=[]
    
        for r in res:
            totres.append(r.res*r.jac)#r.res*r.jac)#r.res*r.jac

        return totres


    def f64_eval_bin(self, x, prec='f64'):

        res_inclusive=0
        res_bin=[0 for i in range(0, self.n_bins+1)]

        [jacques,ks]=self.x_param(x)

        res=NLO(ks[0],ks[1],ks[2],self.a,self.tag,prec,self.phase)

        for r in res:
            r.jac=r.jac*jacques

        if self.mode=='min_pt':
            for r in res:
                ptg=math.sqrt(math.pow(r.pg[1],2.)+math.pow(r.pg[2],2.))
                x1px2=2*(r.j1[0])/self.eCM
                x1mx2=2*(r.j1[3])/self.eCM
                x1=(x1px2+x1mx2)/2
                x2=(x1px2-x1mx2)/2

                if x1<1 and x2<1 and x1>0 and x2>0 and ptg<=self.max_pt and ptg>=self.min_pt:
                #c=1
                #if c==1:
                    res_bin[0]+=r.res*r.jac
                    bin_index=int(ptg/(self.max_pt)*self.n_bins)+1
                    res_bin[bin_index]+=r.res*r.jac

        return res_bin



    def safe_eval_bin(self, x):
        [jacques, ks]=self.x_param(x)
        
        res_inclusive=0
        res_bin=[0 for i in range(0, self.n_bins+1)]

        ksnew=rotator(ks)

        resf64=NLO(ks[0],ks[1],ks[2],self.a,self.tag,'f64',self.phase)
        resnewf64=NLO(ksnew[0],ksnew[1],ksnew[2],self.a,self.tag,'f64',self.phase)

        totresf64=0
        totresnewf64=0
        for (r,newr) in zip(resf64,resnewf64):
            totresf64+=r.res*jacques
            totresnewf64+=newr.res*jacques
            
        c=0
        # print("distance")
        # print(self.n_digits_f64)
        # print(x)
        # print(totresf64)
        # print(totresnewf64)
        # print(abs(totresf64-totresnewf64))
        # print(pow(10.,-self.n_digits_f64)*(abs(totresf64)+abs(totresnewf64)))
        if abs(totresf64-totresnewf64)>pow(10.,-self.n_digits_f64)*(abs(totresf64)+abs(totresnewf64)):
            
            self.n_unstable_f64+=1
            self.unstablef64.append([self.x_param(x),x])

            resf128=NLO(ks[0],ks[1],ks[2],self.a,self.tag,'f128',self.phase)
            resnewf128=NLO(ksnew[0],ksnew[1],ksnew[2],self.a,self.tag,'f128',self.phase)

            totresf128=0
            totresnewf128=0
            for (r,newr) in zip(resf128,resnewf128):
                totresf128+=r.res*jacques
                totresnewf128+=newr.res*jacques

            if abs(totresf128-totresnewf128)<pow(10.,-self.n_digits_f128)*(abs(totresf128)+abs(totresnewf128)):
                #print(x)
                for r in resf128:
                    r.jac=r.jac*jacques

                for r in resf128:
                    ptg=math.sqrt(math.pow(r.pg[1],2.)+math.pow(r.pg[2],2.))
                    x1px2=2*(r.j1[0])/self.eCM
                    x1mx2=2*(r.j1[3])/self.eCM
                    x1=(x1px2+x1mx2)/2
                    x2=(x1px2-x1mx2)/2

                    if x1<1 and x2<1 and x1>0 and x2>0 and ptg<=self.max_pt and ptg>=self.min_pt:
                        res_bin[0]+=r.res*r.jac
                        bin_index=int(ptg/(self.max_pt)*self.n_bins)+1
                        res_bin[bin_index]+=r.res*r.jac

            else:
                self.unstablef128.append(self.x_param(x))
                self.n_unstable_f128+=1
        else:

            for r in resf64:
                r.jac=r.jac*jacques

            for r in resf64:
                ptg=math.sqrt(math.pow(r.pg[1],2.)+math.pow(r.pg[2],2.))
                x1px2=2*(r.j1[0])/self.eCM
                x1mx2=2*(r.j1[3])/self.eCM
                x1=(x1px2+x1mx2)/2
                x2=(x1px2-x1mx2)/2

                a=1
                if x1<1 and x2<1 and x1>0 and x2>0 and ptg<=self.max_pt and ptg>=self.min_pt:
                #if a==1:
                    res_bin[0]+=r.res*r.jac
                    bin_index=int(ptg/(self.max_pt)*self.n_bins)+1
                    res_bin[bin_index]+=r.res*r.jac



        return res_bin










basis=np.linalg.inv(np.array([[1,1,0],[0,1,0],[0,0,1]]))

#my_integrand=NLO_integrands(91.188,13000,0.,0.,1,0.,20,2000,20,basis)
n_iterations_learning=10
n_points_iteration_learning=100000
n_iterations=10
n_points_iteration=100000



topo="none"
#topo="NLO_st"
#print(my_integrand.f64_eval_bin([0.1,0.2,0.3,0.17,0.08,0.4,0.21,0.33,0.13]))


if topo=='limit_checker':
    eCM=2
    a=0.2
    min_pt=0.
    max_pt=6500.
    n_bins=1
    mZ=1
    res_c=0.1
    res_s=0.1
    sigma=1
    basis=np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))


    integrand1=NLO_integrands(mZ, eCM, res_c, res_s, sigma, a, min_pt, max_pt, n_bins, basis, tag='st', debug=False, mode='min_pt',phase='real')

    def collinear_p(eps):
        k=[0.2,0.67,0.73]
        l=[0,eps,0.5]
        klq=[0,0,1]

        q=[]
        for i in range(0,3):
            q.append(klq[i]-k[i]-l[i])

        return [k,l,q]

    points1f64=[]
    points1f128=[]

    for i in range(0,20):
        points1f64.append(integrand1.f64_eval_mom(collinear_p(pow(10,-i/10))))
    for i in range(0,20):
        points1f128.append(integrand1.f128_eval_mom(collinear_p(pow(10,-i/10))))

    print("-------")
    for p in points1f64:
        print(p)
    print("-------")
    for p in points1f128:
        print(p)
    print("-------")

    sump1f64=[sum(p for p in ps) for ps in points1f64]
    sump1f128=[sum(p for p in ps) for ps in points1f128]

    for p in sump1f64:
        print(p)
    print("-------")
    for p in sump1f128:
        print(p)






    integrand2=NLO_integrands(mZ, eCM, res_c, res_s, sigma, a, min_pt, max_pt, n_bins, basis, tag='tqg', debug=False, mode='min_pt',phase='real')

    def collinear_p(eps):
        k=[0,0,0.73]
        l=[0,eps,0.5]
        q=[0.92,0.4,1]

        return [k,l,q]

    points1f64=[]
    points1f128=[]

    for i in range(0,20):
        points1f64.append(integrand2.f64_eval_mom(collinear_p(pow(10,-i/10))))
    for i in range(0,20):
        points1f128.append(integrand2.f128_eval_mom(collinear_p(pow(10,-i/10))))

    print("-------")
    for p in points1f64:
        print(p)
    print("-------")
    for p in points1f128:
        print(p)

    sump1f64=[sum(p for p in ps) for ps in points1f64]
    sump1f128=[sum(p for p in ps) for ps in points1f128]

    print("-------")
    for p in sump1f64:
        print(p)
    print("-------")
    for p in sump1f128:
        print(p)

    

    eCM=2
    a=0.2
    min_pt=0.
    max_pt=6500.
    n_bins=1
    mZ=1
    res_c=0.1
    res_s=0.1
    sigma=1
    basis=np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))


    integrand3=NLO_integrands(mZ, eCM, res_c, res_s, sigma, a, min_pt, max_pt, n_bins, basis, tag='u', debug=False, mode='min_pt',phase='real')

    def collinear_p1(eps):
        k=[0,0,0.73]
        l=[eps,eps,0.5]
        q=[0.92,0.4,1]

        return [k,l,q]

    def collinear_p2(eps):
        kq=[0,0,-0.5]
        l=[eps,eps,0.8]
        q=[0.92,0.4,1]

        k=[]
        for i in range(0,3):
            k.append(kq[i]-q[i])

        return [k,l,q]

    points1f64=[]
    points1f128=[]
    points2f64=[]
    points2f128=[]

    for i in range(0,20):
        points1f64.append(integrand3.f64_eval_mom(collinear_p1(pow(10,-i/10))))
    for i in range(0,20):
        points1f128.append(integrand3.f128_eval_mom(collinear_p1(pow(10,-i/10))))

    print("-------")
    for p in points1f64:
        print(p)
    print("-------")
    for p in points1f128:
        print(p)

    sump1f64=[sum(p for p in ps) for ps in points1f64]
    sump1f128=[sum(p for p in ps) for ps in points1f128]

    print("-------")
    for p in sump1f64:
        print(p)
    print("-------")
    for p in sump1f128:
        print(p)




    for i in range(0,20):
        points2f64.append(integrand3.f64_eval_mom(collinear_p2(pow(10,-i/10))))
    for i in range(0,20):
        points2f128.append(integrand3.f128_eval_mom(collinear_p2(pow(10,-i/10))))

    print("-------")
    for p in points2f64:
        print(p)
    print("-------")
    for p in points2f128:
        print(p)

    sump2f64=[sum(p for p in ps) for ps in points2f64]
    sump2f128=[sum(p for p in ps) for ps in points2f128]

    print("-------")
    for p in sump2f64:
        print(p)
    print("-------")
    for p in sump2f128:
        print(p)



    integrand4=NLO_integrands(mZ, eCM, res_c, res_s, sigma, a, min_pt, max_pt, n_bins, basis, tag='tqq', debug=False, mode='min_pt',phase='real')

    print("set")
    set_sig(1.)
    set_r(0.1,0.1)
    set_kin(91.188,130)

    def collinear_p(eps):
        k=[0,0,0.73]
        l=[0,eps,0.5]
        q=[0.92,0.4,1]

        return [k,l,q]

    pointsf64=[]
    pointsf128=[]


    for i in range(0,30):
        pointsf64.append(pow(10,-i/5)*np.array(integrand4.f64_eval_mom(collinear_p(pow(10,-i/5)))))

    for i in range(0,30):
        pointsf128.append(pow(10,-i/5)*np.array(integrand4.f128_eval_mom(collinear_p(pow(10,-i/5)))))


    print("here")

    print("-------")
    for p in pointsf64:
        print(p)

    print("-----")

    for p in pointsf128:
        print(p)

    sumppf64=[sum(p for p in ps) for ps in pointsf64]
    sumppf128=[sum(p for p in ps) for ps in pointsf128]

    print("-------")
    for p in sumppf64:
        print(p)
    print("-----")

    for p in sumppf128:
        print(p)


#topo="NLO_st"

if topo=='NLO_st':

    my_integrand=NLO_integrands(91.18,13000,0.,0.,2,0.,10,2000,1,basis,tag="st",phase='real')

    my_integrand.set_digits(6,2)

    def f(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        #bins=my_integrand.f64_eval_bin(xs,prec='f64')
        bins=my_integrand.safe_eval_bin(xs)

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        #print(xs)
        #print(res)

        return res

    #print(f([0.584955400549903, 0.007929236073818628, 1, 0.9, 0.438426391965813475, 0.5, 0.5, 0.2, 0, 0.5]))

    
    set_sig(1)

    integ = vegas.Integrator([[0.000,0.99], [0,1], [0,1], [0.000,0.99], [0,1], [0,1], [0.000,0.99], [0,1], [0,1],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

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

#topo="NLO_t_qg"

if topo=='NLO_t_qg':

    my_integrand=NLO_integrands(91.18,13000,0.,0.,2,0.,10,2000,1,basis,tag="tqg",phase='real')

    def f(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand.safe_eval_bin(xs)

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        #print(xs)
        #print(res)

        return res

    #print(f([0.584955400549903, 0.007929236073818628, 1, 0.9, 0.438426391965813475, 0.5, 0.5, 0.2, 0, 0.5]))

    
    integ = vegas.Integrator([[0.0001,0.9999], [0,1], [0,1], [0.0001,0.9999], [0,1], [0,1], [0.0001,0.9999], [0,1], [0,1],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

    integ(f, nitn=n_iterations_learning, neval=n_points_iteration_learning)
    result = integ(f, nitn=n_iterations, neval=n_points_iteration)
    if integ.mpi_rank == 0:
        print(result.summary())
    if integ.mpi_rank == 0:
        print('result = %s    Q = %.2f' % (result[0], result.Q))
        for biny in result:
            print(biny)

#topo="NLO_s"

if topo=='NLO_s':

    my_integrand=NLO_integrands(91.18,13000,0.,0.,2,0.,10,2000,1,basis,tag="s",phase='real')

    def f(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand.safe_eval_bin(xs)

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        #print(xs)
        #print(res)

        return res

    #print(f([0.584955400549903, 0.007929236073818628, 1, 0.9, 0.438426391965813475, 0.5, 0.5, 0.2, 0, 0.5]))

    
    integ = vegas.Integrator([[0.0001,0.9999], [0,1], [0,1], [0.0001,0.9999], [0,1], [0,1], [0.0001,0.9999], [0,1], [0,1],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

    integ(f, nitn=n_iterations_learning, neval=n_points_iteration_learning)
    result = integ(f, nitn=n_iterations, neval=n_points_iteration)
    if integ.mpi_rank == 0:
        print(result.summary())
    if integ.mpi_rank == 0:
        print('result = %s    Q = %.2f' % (result[0], result.Q))
        for biny in result:
            print(biny)

#topo="NLO_t_qq"

if topo=='NLO_t_qq':

    my_integrand=NLO_integrands(91.18,13000,0.,0.,2,0.,10,2000,1,basis,tag="tqq",phase='real')

    def f(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand.safe_eval_bin(xs)

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        #print(xs)
        #print(res)

        return res

    #print(f([0.584955400549903, 0.007929236073818628, 1, 0.9, 0.438426391965813475, 0.5, 0.5, 0.2, 0, 0.5]))

    
    integ = vegas.Integrator([[0.0001,0.9999], [0,1], [0,1], [0.0001,0.9999], [0,1], [0,1], [0.0001,0.9999], [0,1], [0,1],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

    integ(f, nitn=n_iterations_learning, neval=n_points_iteration_learning)
    result = integ(f, nitn=n_iterations, neval=n_points_iteration)
    if integ.mpi_rank == 0:
        print(result.summary())
    if integ.mpi_rank == 0:
        print('result = %s    Q = %.2f' % (result[0], result.Q))
        for biny in result:
            print(biny)


#topo="NLO_u"

if topo=='NLO_u':

    my_integrand=NLO_integrands(91.18,13000,0.,0.,2,0.,10,2000,1,basis,tag="u",phase='real')

    def f(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand.safe_eval_bin(xs)

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        #print(xs)
        #print(res)

        return res

    #print(f([0.584955400549903, 0.007929236073818628, 1, 0.9, 0.438426391965813475, 0.5, 0.5, 0.2, 0, 0.5]))

    
    integ = vegas.Integrator([[0.0001,0.9999], [0,1], [0,1], [0.0001,0.9999], [0,1], [0,1], [0.0001,0.9999], [0,1], [0,1],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

    integ(f, nitn=n_iterations_learning, neval=n_points_iteration_learning)
    result = integ(f, nitn=n_iterations, neval=n_points_iteration)
    if integ.mpi_rank == 0:
        print(result.summary())
    if integ.mpi_rank == 0:
        print('result = %s    Q = %.2f' % (result[0], result.Q))
        for biny in result:
            print(biny)

topo="none"
do="nothing"

if topo=='NLO_dt':

    basisp=np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))

    my_integrand=NLO_integrands(91.18,13000,0.1,0.1,2,0.,0,2000,1,basis,tag="dt",phase='real')

    my_integrand_test=NLO_integrands(1,2,0.1,0.1,2,0.,0,2000,1,basis,tag="dt",phase='real')

    set_defo(-1,0.069,0.99)

    print(my_integrand_test.f64_eval_bin([0.1,0.2,0.3,0.17,0.08,0.27,0.11,0.09,0.3]))
    print(my_integrand_test.f64_eval_bin([0.1,0.2,0.3,0.17,0.08,0.27,0.11,0.09,0.3],prec='f128'))

    if do=="test_ir_limits":
        def collinear_p1(eps):
            k=np.array([0,0,1])
            l=np.array([0,eps,0.5])
            q=np.array([1,0.1,1])

            return [k,l,q]

        def collinear_p2(eps):
            k=[0,0,0.5]
            l=[0,eps,1]
            q=[1,0.1,1]

            return [k,l,q]

        points1=[]
        points2=[]

        for i in range(0,5):
            points1.append(my_integrand_test.f64_eval_mom(collinear_p1(pow(10,-i))))

        for p in points1:
            print(p)

        sump1=[sum(p for p in ps) for ps in points1]

        for p in sump1:
            print(p)


        for i in range(0,5):
            points2.append(my_integrand_test.f64_eval_mom(collinear_p2(pow(10,-i))))

        for p in points2:
            print(p)

        sump2=[sum(p for p in ps) for ps in points2]

        for p in sump2:
            print(p)

    #set_MUV(1.e99)
    set_MUV(1)
    set_kin(1,10)
    set_r(0.1,0.1)


    if do=="test_uv_limits":
        def uv_p1(eps):
            k=[0,0,eps]
            l=[0.17,-0.3,0.5]
            q=[0.34,0.1,1]

            return [k,l,q]

        def uv_p2(eps):
            l=[0,0,eps]
            k=[1,-0.07,0.4]
            q=[0.34,0.1,1]

            return [k,l,q]

        points1=[]
        points1128=[]
        points2=[]
        points2128=[]

        for i in range(0,10):
            points1.append([pow(10,i), pow(10,3*i)*np.array(my_integrand_test.f64_eval_mom(uv_p2(pow(10,i))))])
        for i in range(0,10):
            points1128.append([pow(10,i), pow(10,3*i)*np.array(my_integrand_test.f128_eval_mom(uv_p2(pow(10,i))))])

        for p in points1:
            print("%e : %e %e %e %e"%(tuple([p[0],]+list(p[1]))))
        print("-------")
        for p in points1128:
            print("%e : %e %e %e %e"%(tuple([p[0],]+list(p[1]))))
        print("-------")
        sump1=[sum(p for p in ps[1]) for ps in points1]
        sump1128=[sum(p for p in ps[1]) for ps in points1128]
        print("-------")
        for p in sump1:
            print(p)
        print("-------")
        for p in sump1128:
            print(p)





        # for i in range(0,20):
        #     points2.append([pow(10,i/5), pow(10,3*i/5)*np.array(my_integrand_test.f64_eval_mom(uv_p2(pow(10,i/5))))])
        # for i in range(0,20):
        #     points2128.append([pow(10,i/5), pow(10,3*i/5)*np.array(my_integrand_test.f128_eval_mom(uv_p2(pow(10,i/5))))])

        # for p in points2:
        #     print("%e : %e %e %e %e"%(tuple([p[0],]+list(p[1]))))
        # print("-------")
        # for p in points2128:
        #     print("%e : %e %e %e %e"%(tuple([p[0],]+list(p[1]))))
        # print("-------")


        # sump2128=[sum(p for p in ps) for ps in points2128]
        # sump2=[sum(p for p in ps) for ps in points2]

        # for p in sump2:
        #     print("%e %e"%(p[0], p[1]))
        # print("-------")

        # for p in sump2128:
        #     print("%e %e"%(p[0], p[1]))
        # print("-------")


    if do=="test_soft":
        def soft_p(eps):
            k=[0.1,-0.3+eps,0.5+eps]
            l=[0.1,-0.3,0.5]
            q=[0.34,0.1,1]

            return [k,l,q]


        points1=[]
        points1128=[]
        points2=[]
        points2128=[]

        for i in range(0,10):
            points1.append([pow(10,-i), np.array(my_integrand_test.f64_eval_mom(soft_p(pow(10,-i))))])
        for i in range(0,10):
            points1128.append([pow(10,-i), np.array(my_integrand_test.f128_eval_mom(soft_p(pow(10,-i))))])

        for p in points1:
            print("%e : %e %e %e %e"%(tuple([p[0],]+list(p[1]))))
        print("-------")
        for p in points1128:
            print("%e : %e %e %e %e"%(tuple([p[0],]+list(p[1]))))
        print("-------")
        sump1=[sum(p for p in ps[1]) for ps in points1]
        sump1128=[sum(p for p in ps[1]) for ps in points1128]
        print("-------")
        for p in sump1:
            print(p)
        print("-------")
        for p in sump1128:
            print(p)


    basispp=np.linalg.inv(np.array([[1,-1,0],[0,1,0],[0,0,1]]))
    #basispp=np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))


    my_integrand_imag=NLO_integrands(1,10,0.05,0.05,2,0.,0,2000,1,basispp,tag="dt",phase='imag')
    my_integrand_real=NLO_integrands(1,10,0.05,0.05,2,0.,0,2000,1,basispp,tag="dt",phase='real')
    set_MUV(1)
    my_integrand_real.set_digits(8,2)
    my_integrand_imag.set_digits(8,2)

    set_defo(1,0.069,1)

    print("test")
    for i in range(0,10):
        print(np.array(my_integrand_real.safe_eval_bin([pow(10,-i)*0.1,0.2,0.3,0.1,0.25,0.5,0.24,0.17,0.07,0.5])))
    print("test")
    for i in range(0,10):
        print(np.array(my_integrand_real.f64_eval_bin([pow(10,-i)*0.1,0.2,0.3,0.1,0.25,0.5,0.24,0.17,0.07,0.5])))

    # print("test")
    # for i in range(0,10):
    #     print(my_integrand_real.safe_eval_bin([pow(10,-i)*0.1,0.2,0.3,0.17,0.11,0.67,0.24,0.17,0.07,0.5]))
    #     #print(my_integrand_real.x_param([pow(10,-i)*0.1,0.2,0.3,0.17,0.11,0.67,0.24,0.17,0.07,0.5]))
    # print("test")
    # for i in range(1,10):
    #     print(my_integrand_real.f64_eval_bin([pow(10,-i/10)*0.1,0.2,0.3,0.17,0.2+pow(10,-i/10),0.3+pow(10,-i/10),0.24,0.17,0.07,0.5]))
    

    





    set_defo(1,0.069,0.99)
    #CHANGE DEFO HYPERPARAMETERS MIJ
    #set_MUV(91.188)
    #set_MUV(10)
    #set_MUV(13000)
    set_sig(2)
    #set_kin(91.188,13000)
    set_r(0.1,0.1)
    my_integrand_real.set_digits(12,8)
    my_integrand_imag.set_digits(12,8)

    n_iterations_learning=10
    n_points_iteration_learning=100000
    n_iterations=10
    n_points_iteration=200000

    


    def freal(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand_real.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand_real.safe_eval_bin(xs)
        #bins=my_integrand_real.f64_eval_bin(xs,prec="f128")

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        #print(xs)
        #print(res)

        return res

    def fimag(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand_imag.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand_imag.safe_eval_bin(xs)

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        return res

    #print(f([0.584955400549903, 0.007929236073818628, 1, 0.9, 0.438426391965813475, 0.5, 0.5, 0.2, 0, 0.5]))

    if do=="integrate_real" or do=="integrate_both":    
        integ = vegas.Integrator([[0.001,0.9999], [0.0,1.0], [0.0,1.0], [0.001,0.9999], [0.0,1.0], [0.0,1.0], [0.001,0.9999], [0.0,1.0], [0.0,1.0],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

        integ(freal, nitn=n_iterations_learning, neval=n_points_iteration_learning)
        result = integ(freal, nitn=n_iterations, neval=n_points_iteration)
        if integ.mpi_rank == 0:
            print(result.summary())
        if integ.mpi_rank == 0:
            print('result = %s    Q = %.2f' % (result[0], result.Q))
            for biny in result:
                print(biny)
            print(my_integrand_real.n_unstable_f128)
            for i in range(0,1):
                print(my_integrand_real.unstablef128[i])


    if do=="integrate_imag" or do=="integrate_both":    
        integ = vegas.Integrator([[0.001,0.9999], [0.0,1.0], [0.0,1.0], [0.001,0.9999], [0.0,1.0], [0.0,1.0], [0.001,0.9999], [0.0,1.0], [0.0,1.0],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

        integ(fimag, nitn=n_iterations_learning, neval=n_points_iteration_learning)
        result = integ(fimag, nitn=n_iterations, neval=n_points_iteration)
        if integ.mpi_rank == 0:
            print(result.summary())
        if integ.mpi_rank == 0:
            print('result = %s    Q = %.2f' % (result[0], result.Q))
            for biny in result:
                print(biny)
            

if topo=='NLO_dt_real':

    basispp=np.linalg.inv(np.array([[1,-1,0],[0,1,0],[0,0,1]]))
    my_integrand_real=NLO_integrands(1,10,0.05,0.05,2,0.,0,2000,1,basispp,tag="dt",phase='real')

    set_MUV(1)

    #set_defo(0.5,0.069,1)
    set_defo(0.5,0.1,1)
    set_sig(2)
    set_r(0.1,0.1)
    my_integrand_real.set_digits(12,8)

    n_iterations_learning=10
    n_points_iteration_learning=100000
    n_iterations=10
    n_points_iteration=10000000



    def freal(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand_real.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand_real.safe_eval_bin(xs)
        #bins=my_integrand_real.f64_eval_bin(xs,prec="f128")

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        return res


 
    integ = vegas.Integrator([[0.001,0.9999], [0.0,1.0], [0.0,1.0], [0.001,0.9999], [0.0,1.0], [0.0,1.0], [0.001,0.9999], [0.0,1.0], [0.0,1.0],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

    integ(freal, nitn=n_iterations_learning, neval=n_points_iteration_learning)
    result = integ(freal, nitn=n_iterations, neval=n_points_iteration)
    if integ.mpi_rank == 0:
        print(result.summary())
    if integ.mpi_rank == 0:
        print('result = %s    Q = %.2f' % (result[0], result.Q))
        for biny in result:
            print(biny)
        print(my_integrand_real.n_unstable_f128)
        for i in range(0,1):
            print(my_integrand_real.unstablef128[i])
            
    


if topo=='NLO_se':

    basisp=np.linalg.inv(np.array([[1,0,0],[0,1,0],[0,0,1]]))

    #my_integrand=NLO_integrands(91.18,13000,0.1,0.1,2,0.,0,2000,1,basis,tag="dt",phase='real')

    my_integrand_test=NLO_integrands(1,2,0.1,0.1,2,0.,0,2000,1,basis,tag="se",phase='real')

    set_defo(-1,0.069,-0.99)

    #print(my_integrand.f64_eval_bin([0.1,0.2,0.3,0.17,0.08,0.27,0.11,0.09,0.3]))

    if do=="test_ir_limits":
        def collinear_p1(eps):
            k=[0,0,1]
            l=[0,eps,0.5]
            q=[1,0.1,1]

            return [k,l,q]

        points1=[]

        for i in range(0,30):
            points1.append(pow(10,-i/10)*np.array(my_integrand_test.f64_eval_mom(collinear_p1(pow(10,-i/10)))))

        for p in points1:
            print(p)

        sump1=[sum(p for p in ps) for ps in points1]

        for p in sump1:
            print(p)




    if do=="test_uv_limits":

        def uv_p2(eps):
            l=[0,0,eps]
            k=[0.17,-0.3,0.5]
            q=[0.34,0.1,1]

            return [k,l,q]

        points2=[]

        for i in range(0,40):
            points2.append(pow(10,3*i/10)*np.array(my_integrand_test.f64_eval_mom(uv_p2(pow(10,i/10)))))

        for p in points2:
            print(np.array(p))

        sump2=[sum(p for p in ps) for ps in points2]

        for p in sump2:
            print(p)

    #strp=['{:.16f}'.format(p) for p in points]

    #print("List1", '[%s]'%', '.join(map(str, strp)))


    basisp=np.linalg.inv(np.array([[1,-1,0],[0,1,0],[0,0,1]]))

    my_integrand_imag=NLO_integrands(91.18,13000,0.05,0.05,2,0.,0,2000,1,basisp,tag="se",phase='imag')
    my_integrand_real=NLO_integrands(91.18,13000,0.05,0.05,2,0.,0,2000,1,basisp,tag="se",phase='real')

    set_defo(1.,0.069,0.99)

    set_MUV(1)


    def freal(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand_real.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand_real.f64_eval_bin(xs)

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        #print(xs)
        #print(res)

        return res

    def fimag(x):

        ap=10*math.tan(math.pi*(x[9]-1/2))
        my_integrand_imag.set_a(ap)

        xs=[x[i] for i in range(0,9)]

        bins=my_integrand_imag.f64_eval_bin(xs)

        res=[biny*10*math.pi/pow(math.sin(math.pi*x[9]),2) for biny in bins]

        #print(xs)
        #print(res)

        return res

    #print(f([0.584955400549903, 0.007929236073818628, 1, 0.9, 0.438426391965813475, 0.5, 0.5, 0.2, 0, 0.5]))

    if do=="integrate_real" or do=="integrate_both":    
        integ = vegas.Integrator([[0.001,0.999], [0,1], [0,1], [0.001,0.999], [0,1], [0,1], [0.001,0.999], [0,1], [0,1],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

        integ(freal, nitn=n_iterations_learning, neval=n_points_iteration_learning)
        result = integ(freal, nitn=n_iterations, neval=n_points_iteration)
        if integ.mpi_rank == 0:
            print(result.summary())
        if integ.mpi_rank == 0:
            print('result = %s    Q = %.2f' % (result[0], result.Q))
            for biny in result:
                print(biny)

    if do=="integrate_imag" or do=="integrate_both":    
        integ = vegas.Integrator([[0.001,0.999], [0,1], [0,1], [0.001,0.999], [0,1], [0,1], [0.001,0.999], [0,1], [0,1],[0,1]]) #, [r3_min, r3_max], [th3_min, th3_max], [ph3_min, ph3_max]])

        integ(fimag, nitn=n_iterations_learning, neval=n_points_iteration_learning)
        result = integ(fimag, nitn=n_iterations, neval=n_points_iteration)
        if integ.mpi_rank == 0:
            print(result.summary())
        if integ.mpi_rank == 0:
            print('result = %s    Q = %.2f' % (result[0], result.Q))
            for biny in result:
                print(biny)