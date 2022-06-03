import math
from random import random
import numpy as np
from scipy.spatial.transform import Rotation as Roto
from integrand_wrapper import cut_res, LO, NLO
from integrand_wrapper import set_r, set_kin, set_sig, set_defo, set_MUV, set_observable, set_observable_lambda


def rotator(ks):
    nx=random()
    ny=random()
    nz=random()

    theta=random()*2*math.pi

    axis=theta/(math.sqrt(nx*nx+ny*ny+nz*nz))*np.array([nx,ny,nz])

    r = Roto.from_rotvec(axis)

    rotks=[r.apply(ks[0]),r.apply(ks[1]),r.apply(ks[2])]

    return rotks


class LO_integrands:

    def __init__(self, eCM, m, x1, n_bins, tag='dy'):
        self.eCM = eCM
        self.m = m
        self.x1 = x1
        self.x2 = m**2/(4*(eCM**2)*x1)
        self.jac = (self.x1+self.x2)*self.x2/(4*pow(self.x1*self.x2,3/2))
        self.a = (self.x1-self.x2)/(2*math.sqrt(self.x1*self.x2))
        self.n_bins = n_bins
        self.tag = tag
        set_kin(m,eCM)

    def set_x1(self, x1):
        self.x1 = x1
        self.x2 = self.m**2/(4*(self.eCM**2)*x1)
        self.jac = (self.x1+self.x2)*self.x2/(4*pow(self.x1*self.x2,3/2))
        self.a = (self.x1-self.x2)/(2*math.sqrt(self.x1*self.x2))

    def set_a_jj(self,a):
        self.a=a
        self.x1=0.5
        self.x2=0.5
        self.jac=1


    def LO_bin(self, x):
        res=0

        dr1=x[0]/(1-x[0])
        dth1=x[1]*2*math.pi
        dph1=x[2]*math.pi
        dr2=x[3]/(1-x[3])
        dth2=x[4]*2*math.pi
        dph2=x[5]*math.pi

        p=[dr1*math.cos(dth1)*math.sin(dph1),dr1*math.sin(dth1)*math.sin(dph1),dr1*math.cos(dph1)]
        q=[dr2*math.cos(dth2)*math.sin(dph2),dr2*math.sin(dth2)*math.sin(dph2),dr2*math.cos(dph2)]

        if self.x2<1 and self.x2>0:
            res=LO(p,q,self.a,'dy','f64')

        jacques=math.pow(dr1,2)*math.sin(dph1)*math.pow(dr2,2)*math.sin(dph2)*math.pow(2*math.pi,2)*math.pow(math.pi,2)/(math.pow(1-x[0],2)*math.pow(1-x[3],2))

        return self.jac*jacques*res


# class LO_integrands:

#     def __init__(self, m, eCM, sigma, a, debug=False):
#         self.eCM = eCM
#         self.mZ=m
#         self.a = a
#         set_kin(m,eCM)
#         set_sig(sigma)
#         set_r(0.,0.)

#     def set_a(self, ma):
#         self.a = ma

#     def eval(self, x):
#         eCMn=1.
#         dr1=eCMn*x[0]/(1-x[0])
#         dth1=x[1]*2*math.pi
#         dph1=x[2]*math.pi
#         dr2=eCMn*x[3]/(1-x[3])
#         dth2=x[4]*2*math.pi
#         dph2=x[5]*math.pi

#         k=np.array([dr1*math.cos(dth1)*math.sin(dph1),dr1*math.sin(dth1)*math.sin(dph1),dr1*math.cos(dph1)])
#         q=np.array([dr2*math.cos(dth2)*math.sin(dph2),dr2*math.sin(dth2)*math.sin(dph2),dr2*math.cos(dph2)])

#         jacques=pow(eCMn,2)*math.pow(dr1,2)*math.sin(dph1)*math.pow(dr2,2)*math.sin(dph2)*math.pow(2*math.pi,2)*math.pow(math.pi,2)/(math.pow(1-x[0],2)*math.pow(1-x[3],2))

#         j1=[0,0,0,0]
#         j2=[0,0,0,0]

#         res=cut_res(LO(k,q,self.a,'dy'),jacques,j1,1,j2,-1,[0,0,0,0])

#         return res


class NLO_integrands:

    def __init__(self, m, eCM, res_c, res_s, sigma, a, min_pt, max_pt, n_bins, basis, tag='st', obs_tag='JADE', obs_lambda=0,debug=False, mode='min_pt',phase='real'):
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
        self.max_val=0
        self.max_x=0
        self.max_k=0
        if obs_tag=='JADE':
            set_observable(0)
        if obs_tag=='HEMI':
            set_observable(1)
            set_observable_lambda(obs_lambda)
        if obs_tag=='JJ':
            set_observable(2)
        if obs_tag!='JADE' and obs_tag!='HEMI' and obs_tag!='JJ':
            print("not a valid clustering algorithm")
        


    def set_digits(self,n_digits_f64,n_digits_f128):
        self.n_digits_f128=n_digits_f128
        self.n_digits_f64=n_digits_f64

    def set_a(self, a):
        self.a = a

    def set_basis(self, basis):
        self.basis = basis

    def set_obs(obs_tag, lamb):
        if obs_tag=='JADE':
            set_observable(0)
        if obs_tag=='HEMI':
            set_observable(1)
            set_observable_lambda(obs_lambda)
        if obs_tag=='JJ':
            set_observable(2)
        if obs_tag!='JADE' and obs_tag!='HEMI' and obs_tag!='JJ':
            print("not a valid clustering algorithm")


    def x_param(self, x):
        eCMn=self.eCM/10.
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

        # print("jacques x")
        # print(jacques)
        # print(dr1)
        # print(dr2)
        # print(dr3)

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


    def sample_f64_mom(self, moms):
        res=NLO(moms[0],moms[1],moms[2],self.a,self.tag,'f64',self.phase)
        for r in res:
            r.res=r.res*r.jac#r.res*r.jac)#*r.jac)#r.res*r.jac
        #print(res)
        return res

    def sample_f128_mom(self, moms):

        res=NLO(moms[0],moms[1],moms[2],self.a,self.tag,'f128',self.phase)
    
        for r in res:
            r.res=r.res*r.jac

        return res


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
                x1px2=2*(r.j1[0]+r.j2[0])/self.eCM
                x1mx2=2*(r.j1[3]+r.j2[3])/self.eCM

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

        resf64=[]
        resnewf64=[]

        if self.tag=='all_trees':
            alltags=['st','tqq','tqg','s','u']
            for tags in alltags:
                for (v1,v2) in zip(NLO(ks[0],ks[1],ks[2],self.a,tags,'f64',self.phase),NLO(ksnew[0],ksnew[1],ksnew[2],self.a,tags,'f64',self.phase)):
                    resf64.append(v1)
                    resnewf64.append(v2)
        else:
            resf64=NLO(ks[0],ks[1],ks[2],self.a,self.tag,'f64',self.phase)
            resnewf64=NLO(ksnew[0],ksnew[1],ksnew[2],self.a,self.tag,'f64',self.phase)

        # resf64=NLO(ks[0],ks[1],ks[2],self.a,self.tag,'f64',self.phase)
        # resnewf64=NLO(ksnew[0],ksnew[1],ksnew[2],self.a,self.tag,'f64',self.phase)

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
            self.unstablef64.append([x, self.x_param(x)])

            resf128=[]
            resnewf128=[]

            if self.tag=='all_trees':
                alltags=['st','tqq','tqg','s','u']
                for tags in alltags:
                    for (v1,v2) in zip(NLO(ks[0],ks[1],ks[2],self.a,tags,'f128',self.phase),NLO(ksnew[0],ksnew[1],ksnew[2],self.a,tags,'f128',self.phase)):
                        resf128.append(v1)
                        resnewf128.append(v2)
            else:
                resf128=NLO(ks[0],ks[1],ks[2],self.a,self.tag,'f128',self.phase)
                resnewf128=NLO(ksnew[0],ksnew[1],ksnew[2],self.a,self.tag,'f128',self.phase)

            # resf128=NLO(ks[0],ks[1],ks[2],self.a,self.tag,'f128',self.phase)
            # resnewf128=NLO(ksnew[0],ksnew[1],ksnew[2],self.a,self.tag,'f128',self.phase)

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
                    x1px2=2*(r.j1[0]+r.j2[0])/self.eCM
                    x1mx2=2*(r.j1[3]+r.j2[3])/self.eCM

                    x1=(x1px2+x1mx2)/2
                    x2=(x1px2-x1mx2)/2

                    if x1<1 and x2<1 and x1>0 and x2>0 and ptg<=self.max_pt and ptg>=self.min_pt:
                        res_bin[0]+=r.res*r.jac
                        bin_index=int(ptg/(self.max_pt)*self.n_bins)+1
                        res_bin[bin_index]+=r.res*r.jac

            else:

                self.unstablef128.append([x,self.x_param(x)])
                self.n_unstable_f128+=1


        else:

            for r in resf64:
                r.jac=r.jac*jacques

            for r in resf64:
                ptg=math.sqrt(math.pow(r.pg[1],2.)+math.pow(r.pg[2],2.))
                x1px2=2*(r.j1[0]+r.j2[0])/self.eCM
                x1mx2=2*(r.j1[3]+r.j2[3])/self.eCM

                x1=(x1px2+x1mx2)/2
                x2=(x1px2-x1mx2)/2

                a=1
                if x1<1 and x2<1 and x1>0 and x2>0 and ptg<=self.max_pt and ptg>=self.min_pt:
                #if a==1:
                    res_bin[0]+=r.res*r.jac
                    bin_index=int(ptg/(self.max_pt)*self.n_bins)+1
                    res_bin[bin_index]+=r.res*r.jac



        return res_bin



    def safe_eval(self, x):
        [jacques, ks]=self.x_param(x)
        
        res_inclusive=0
        res_bin=[0 for i in range(0, self.n_bins+1)]

        ksnew=rotator(ks)

        resf64=[]
        resnewf64=[]

        if self.tag=='all_trees':
            alltags=['st','tqq','tqg','s','u']
            for tags in alltags:
                for (v1,v2) in zip(NLO(ks[0],ks[1],ks[2],self.a,tags,'f64',self.phase),NLO(ksnew[0],ksnew[1],ksnew[2],self.a,tags,'f64',self.phase)):
                    resf64.append(v1)
                    resnewf64.append(v2)
        else:
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

            resf128=[]
            resnewf128=[]

            if self.tag=='all_trees':
                alltags=['st','tqq','tqg','s','u']
                for tags in alltags:
                    for (v1,v2) in zip(NLO(ks[0],ks[1],ks[2],self.a,tags,'f128',self.phase),NLO(ksnew[0],ksnew[1],ksnew[2],self.a,tags,'f128',self.phase)):
                        resf128.append(v1)
                        resnewf128.append(v2)
            else:
                resf128=NLO(ks[0],ks[1],ks[2],self.a,self.tag,'f128',self.phase)
                resnewf128=NLO(ksnew[0],ksnew[1],ksnew[2],self.a,self.tag,'f128',self.phase)

            # resf128=NLO(ks[0],ks[1],ks[2],self.a,self.tag,'f128',self.phase)
            # resnewf128=NLO(ksnew[0],ksnew[1],ksnew[2],self.a,self.tag,'f128',self.phase)

            totresf128=0
            totresnewf128=0
            for (r,newr) in zip(resf128,resnewf128):
                totresf128+=r.res*jacques
                totresnewf128+=newr.res*jacques

            # print("****************")
            # print(totresf64)
            # print(totresnewf64)
            # print("----")
            # print(totresf128)
            # print(totresnewf128)

            if abs(totresf128-totresnewf128)<pow(10.,-self.n_digits_f128)*(abs(totresf128)+abs(totresnewf128)):
                #print(x)
                res=0
                for r in resf128:
                    r.jac=r.jac*jacques
                    res+=r.res*r.jac
                
                if abs(res)>abs(self.max_val):
                    self.max_val=res
                    self.max_x=x
                    self.max_k=ks

                return resf128

            else:
                for r in resf128:
                    r.res=0
                self.unstablef128.append(self.x_param(x))
                self.n_unstable_f128+=1
                return resf128
        else:

            res=0
            for r in resf64:
                r.jac=r.jac*jacques
                res+=r.res*r.jac
                
            if abs(res)>abs(self.max_val):
                self.max_val=res
                self.max_x=x
                self.max_k=ks
                

            return resf64

        res=0
        for r in resnewf64:
            r.res=0
            res+=r.res*r.jac
                
        if abs(res)>abs(self.max_val):
            self.max_val=res
            self.max_x=x
            self.max_k=ks
            
        return resnewf64
        

