from ctypes import *



###################SET PARAMETERS##################


#########F64########

set_kinematics = cdll.LoadLibrary('./integrand.so').set_kinematics
set_res = cdll.LoadLibrary('./integrand.so').set_res
set_sigma = cdll.LoadLibrary('./integrand.so').set_sigma
set_defo_parameters = cdll.LoadLibrary('./integrand.so').set_defo_parameters
set_uv_mass = cdll.LoadLibrary('./integrand.so').set_uv_mass

set_observable_tag = cdll.LoadLibrary('./integrand.so').set_observable_tag
set_hemisphere_lambda = cdll.LoadLibrary('./integrand.so').set_hemisphere_lambda


set_kinematics.argtypes = [c_double, c_double]
set_kinematics.restype = None

set_res.argtypes = [c_double, c_double]
set_res.restype = None

set_sigma.argtypes = [c_double]
set_sigma.restype = None

set_defo_parameters.argtypes = [c_double, c_double, c_double]
set_defo_parameters.restype = None

set_uv_mass.argtypes = [c_double]
set_uv_mass.restype = None


#########F128########

set_kinematics_f128 = cdll.LoadLibrary('./integrand_f128.so').set_kinematics
set_res_f128 = cdll.LoadLibrary('./integrand_f128.so').set_res
set_sigma_f128 = cdll.LoadLibrary('./integrand_f128.so').set_sigma
set_defo_parameters_f128 = cdll.LoadLibrary('./integrand_f128.so').set_defo_parameters
set_uv_mass_f128 = cdll.LoadLibrary('./integrand_f128.so').set_uv_mass

set_observable_tag_f128 = cdll.LoadLibrary('./integrand.so').set_observable_tag
set_hemisphere_lambda_f128 = cdll.LoadLibrary('./integrand.so').set_hemisphere_lambda


set_kinematics_f128.argtypes = [c_double, c_double]
set_kinematics_f128.restype = None

set_res_f128.argtypes = [c_double, c_double]
set_res_f128.restype = None

set_sigma_f128.argtypes = [c_double]
set_sigma_f128.restype = None

set_defo_parameters_f128.argtypes = [c_double, c_double, c_double]
set_defo_parameters_f128.restype = None

set_uv_mass_f128.argtypes = [c_double]
set_uv_mass_f128.restype = None



def set_kin(m,eCM):
    set_kinematics(m,eCM)
    set_kinematics_f128(m,eCM)

def set_r(res_c,res_s):
    set_res(res_c,res_s)
    set_res_f128(res_c,res_s)

def set_sig(sigma):
    set_sigma(sigma)
    set_sigma_f128(sigma)

def set_defo(lamb, mij, bc_lamb):
    set_defo_parameters(lamb, mij, bc_lamb)
    set_defo_parameters_f128(lamb, mij, bc_lamb)

def set_MUV(nMUV):
    set_uv_mass(nMUV)
    set_uv_mass_f128(nMUV)

def set_observable(tag):
    set_observable_tag(tag)
    set_observable_tag_f128(tag)

def set_observable_lambda(lamb):
    set_hemisphere_lambda(lamb)
    set_hemisphere_lambda_f128(lamb)




###################INTEGRANDS#####################


########F64########


LO_scalar_eval = cdll.LoadLibrary('./integrand.so').LO_scalar_eval
LO_DrellYan_eval = cdll.LoadLibrary('./integrand.so').LO_DrellYan_eval
NLO_st_channel_DrellYan_eval = cdll.LoadLibrary('./integrand.so').NLO_st_channel_DrellYan_eval
NLO_t_channel_qq_DrellYan_eval = cdll.LoadLibrary('./integrand.so').NLO_t_channel_qq_DrellYan_eval
NLO_t_channel_qg_DrellYan_eval = cdll.LoadLibrary('./integrand.so').NLO_t_channel_qg_DrellYan_eval
NLO_u_channel_DrellYan_eval = cdll.LoadLibrary('./integrand.so').NLO_u_channel_DrellYan_eval
NLO_s_channel_DrellYan_eval = cdll.LoadLibrary('./integrand.so').NLO_s_channel_DrellYan_eval
NLO_dt_DrellYan_eval = cdll.LoadLibrary('./integrand.so').NLO_dt_DrellYan_eval
NLO_se_DrellYan_eval = cdll.LoadLibrary('./integrand.so').NLO_se_DrellYan_eval
#NLO_dt_jj_eval = cdll.LoadLibrary('./integrand.so').NLO_dt_jj_eval
#NLO_se_jj_eval = cdll.LoadLibrary('./integrand.so').NLO_se_jj_eval
get_res_eval = cdll.LoadLibrary('./integrand.so').get_res_eval
get_res_complex_eval = cdll.LoadLibrary('./integrand.so').get_res_complex_eval



LO_scalar_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double]
LO_scalar_eval.restype = c_double

LO_DrellYan_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double]
LO_DrellYan_eval.restype = c_double

NLO_st_channel_DrellYan_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_st_channel_DrellYan_eval.restype = None

NLO_t_channel_qq_DrellYan_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_t_channel_qq_DrellYan_eval.restype = None

NLO_t_channel_qg_DrellYan_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_t_channel_qg_DrellYan_eval.restype = None

NLO_u_channel_DrellYan_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_u_channel_DrellYan_eval.restype = None

NLO_s_channel_DrellYan_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_s_channel_DrellYan_eval.restype = None

NLO_dt_DrellYan_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_dt_DrellYan_eval.restypes = None

NLO_se_DrellYan_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_se_DrellYan_eval.restypes = None

#NLO_dt_jj_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int, c_int]
#NLO_dt_jj_eval.restype = None

#NLO_se_jj_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int, c_int]
#NLO_se_jj_eval.restype = None

get_res_eval.argtypes = [c_int, c_int]
get_res_eval.restype = c_double

get_res_complex_eval.argtypes = [c_int, c_int, c_int]
get_res_complex_eval.restype = c_double




########F128########


LO_scalar_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').LO_scalar_eval
LO_DrellYan_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').LO_DrellYan_eval
NLO_st_channel_DrellYan_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').NLO_st_channel_DrellYan_eval
NLO_t_channel_qq_DrellYan_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').NLO_t_channel_qq_DrellYan_eval
NLO_t_channel_qg_DrellYan_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').NLO_t_channel_qg_DrellYan_eval
NLO_u_channel_DrellYan_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').NLO_u_channel_DrellYan_eval
NLO_s_channel_DrellYan_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').NLO_s_channel_DrellYan_eval
NLO_dt_DrellYan_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').NLO_dt_DrellYan_eval
NLO_se_DrellYan_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').NLO_se_DrellYan_eval
#NLO_dt_jj_eval = cdll.LoadLibrary('./integrand.so').NLO_dt_jj_eval
#NLO_se_jj_eval = cdll.LoadLibrary('./integrand.so').NLO_se_jj_eval
get_res_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').get_res_eval
get_res_complex_eval_f128 = cdll.LoadLibrary('./integrand_f128.so').get_res_complex_eval





LO_scalar_eval_f128.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double]
LO_scalar_eval_f128.restype = c_double

LO_DrellYan_eval_f128.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double]
LO_DrellYan_eval_f128.restype = c_double

NLO_st_channel_DrellYan_eval_f128.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_st_channel_DrellYan_eval_f128.restype = None

NLO_t_channel_qq_DrellYan_eval_f128.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_t_channel_qq_DrellYan_eval_f128.restype = None

NLO_t_channel_qg_DrellYan_eval_f128.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_t_channel_qg_DrellYan_eval_f128.restype = None

NLO_u_channel_DrellYan_eval_f128.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_u_channel_DrellYan_eval_f128.restype = None

NLO_s_channel_DrellYan_eval_f128.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_s_channel_DrellYan_eval_f128.restype = None

NLO_dt_DrellYan_eval_f128.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_dt_DrellYan_eval_f128.restypes = None

NLO_se_DrellYan_eval_f128.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
NLO_se_DrellYan_eval_f128.restypes = None

#NLO_dt_jj_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int, c_int]
#NLO_dt_jj_eval.restype = None

#NLO_se_jj_eval.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int, c_int]
#NLO_se_jj_eval.restype = None

get_res_eval_f128.argtypes = [c_int, c_int]
get_res_eval_f128.restype = c_double

get_res_complex_eval_f128.argtypes = [c_int, c_int, c_int]
get_res_complex_eval_f128.restype = c_double




class cut_res:
    def __init__(self, res, jac, j1, spin1, j2, spin2, pg):
        self.res = res
        self.jac = jac
        self.j1 = j1
        self.spin1 = spin1
        self.j2 = j2
        self.spin2 = spin2
        self.pg = pg

    @classmethod
    def set_equal(cls, cut_res_n):
        res = cut_res_n.res
        jac = cut_res_n.res
        j1 = cut_res_n.res
        spin1 = cut_res_n.res
        j2 = cut_res_n.res
        spin2 = cut_res_n.res
        pg = cut_res_n.res

        return cls(res,jac,j1,spin1,j2,spin2,pg)




def LO(k,l,a,tag='dy',prec='f64'):

    resf=0

    if tag=='scalar':
        if prec=='f64':
            resf=LO_scalar_eval(k[0],k[1],k[2],l[0],l[1],l[2],a)

    if tag=='dy':
        if prec=='f64':
            resf=LO_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],a)

    return resf


def NLO(k,l,q,a,tag='st',prec='f64',phase='real'):

    resf=[]

    if tag=='s':
        if prec=='f64':
            NLO_s_channel_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

        if prec=='f128':
            NLO_s_channel_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

    if tag=='st':
        if prec=='f64':
            NLO_st_channel_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

            NLO_st_channel_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

        if prec=='f128':
            NLO_st_channel_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

            NLO_st_channel_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

    if tag=='tqq':
        if prec=='f64':
            NLO_t_channel_qq_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

            NLO_t_channel_qq_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

        if prec=='f128':
            NLO_t_channel_qq_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

            NLO_t_channel_qq_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

    if tag=='tqg':
        if prec=='f64':
            NLO_t_channel_qg_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

            NLO_t_channel_qg_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

        if prec=='f128':
            NLO_t_channel_qg_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

            NLO_t_channel_qg_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

    if tag=='u':
        if prec=='f64':
            NLO_u_channel_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

            NLO_u_channel_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

            NLO_u_channel_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,3)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

            NLO_u_channel_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,4)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

            NLO_u_channel_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,5)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

            NLO_u_channel_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,6)
            j1=[get_res_eval(2,0),get_res_eval(2,1),get_res_eval(2,2),get_res_eval(2,3)]
            j2=[get_res_eval(4,0),get_res_eval(4,1),get_res_eval(4,2),get_res_eval(4,3)]
            ptg=[get_res_eval(6,0),get_res_eval(6,1),get_res_eval(6,2),get_res_eval(6,3)]
            resf.append(cut_res(get_res_eval(0,0),get_res_eval(1,0),j1,get_res_eval(3,0),j2,get_res_eval(5,0),ptg))

        if prec=='f128':
            NLO_u_channel_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

            NLO_u_channel_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

            NLO_u_channel_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,3)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

            NLO_u_channel_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,4)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

            NLO_u_channel_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,5)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

            NLO_u_channel_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,6)
            j1=[get_res_eval_f128(2,0),get_res_eval_f128(2,1),get_res_eval_f128(2,2),get_res_eval_f128(2,3)]
            j2=[get_res_eval_f128(4,0),get_res_eval_f128(4,1),get_res_eval_f128(4,2),get_res_eval_f128(4,3)]
            ptg=[get_res_eval_f128(6,0),get_res_eval_f128(6,1),get_res_eval_f128(6,2),get_res_eval_f128(6,3)]
            resf.append(cut_res(get_res_eval_f128(0,0),get_res_eval_f128(1,0),j1,get_res_eval_f128(3,0),j2,get_res_eval_f128(5,0),ptg))

    if tag=='dt':
        if prec=='f64':
            if phase=='real':
                NLO_dt_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
                j1=[get_res_complex_eval(2,0,0),get_res_complex_eval(2,1,0),get_res_complex_eval(2,2,0),get_res_complex_eval(2,3,0)]
                j2=[get_res_complex_eval(4,0,0),get_res_complex_eval(4,1,0),get_res_complex_eval(4,2,0),get_res_complex_eval(4,3,0)]
                ptg=[get_res_complex_eval(6,0,0),get_res_complex_eval(6,1,0),get_res_complex_eval(6,2,0),get_res_complex_eval(6,3,0)]
                resf.append(cut_res(get_res_complex_eval(0,0,0),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,0),j2,get_res_complex_eval(5,0,0),ptg))

                NLO_dt_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
                j1=[get_res_complex_eval(2,0,0),get_res_complex_eval(2,1,0),get_res_complex_eval(2,2,0),get_res_complex_eval(2,3,0)]
                j2=[get_res_complex_eval(4,0,0),get_res_complex_eval(4,1,0),get_res_complex_eval(4,2,0),get_res_complex_eval(4,3,0)]
                ptg=[get_res_complex_eval(6,0,0),get_res_complex_eval(6,1,0),get_res_complex_eval(6,2,0),get_res_complex_eval(6,3,0)]
                resf.append(cut_res(get_res_complex_eval(0,0,0),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,0),j2,get_res_complex_eval(5,0,0),ptg))

                NLO_dt_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,3)
                j1=[get_res_complex_eval(2,0,0),get_res_complex_eval(2,1,0),get_res_complex_eval(2,2,0),get_res_complex_eval(2,3,0)]
                j2=[get_res_complex_eval(4,0,0),get_res_complex_eval(4,1,0),get_res_complex_eval(4,2,0),get_res_complex_eval(4,3,0)]
                ptg=[get_res_complex_eval(6,0,0),get_res_complex_eval(6,1,0),get_res_complex_eval(6,2,0),get_res_complex_eval(6,3,0)]
                resf.append(cut_res(get_res_complex_eval(0,0,0),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,0),j2,get_res_complex_eval(5,0,0),ptg))

                NLO_dt_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,4)
                j1=[get_res_complex_eval(2,0,0),get_res_complex_eval(2,1,0),get_res_complex_eval(2,2,0),get_res_complex_eval(2,3,0)]
                j2=[get_res_complex_eval(4,0,0),get_res_complex_eval(4,1,0),get_res_complex_eval(4,2,0),get_res_complex_eval(4,3,0)]
                ptg=[get_res_complex_eval(6,0,0),get_res_complex_eval(6,1,0),get_res_complex_eval(6,2,0),get_res_complex_eval(6,3,0)]
                resf.append(cut_res(get_res_complex_eval(0,0,0),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,0),j2,get_res_complex_eval(5,0,0),ptg))

            if phase=='imag':
                NLO_dt_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
                j1=[get_res_complex_eval(2,0,1),get_res_complex_eval(2,1,1),get_res_complex_eval(2,2,1),get_res_complex_eval(2,3,1)]
                j2=[get_res_complex_eval(4,0,1),get_res_complex_eval(4,1,1),get_res_complex_eval(4,2,1),get_res_complex_eval(4,3,1)]
                ptg=[get_res_complex_eval(6,0,1),get_res_complex_eval(6,1,1),get_res_complex_eval(6,2,1),get_res_complex_eval(6,3,1)]
                resf.append(cut_res(get_res_complex_eval(0,0,1),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,1),j2,get_res_complex_eval(5,0,1),ptg))

                NLO_dt_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
                j1=[get_res_complex_eval(2,0,1),get_res_complex_eval(2,1,1),get_res_complex_eval(2,2,1),get_res_complex_eval(2,3,1)]
                j2=[get_res_complex_eval(4,0,1),get_res_complex_eval(4,1,1),get_res_complex_eval(4,2,1),get_res_complex_eval(4,3,1)]
                ptg=[get_res_complex_eval(6,0,1),get_res_complex_eval(6,1,1),get_res_complex_eval(6,2,1),get_res_complex_eval(6,3,1)]
                resf.append(cut_res(get_res_complex_eval(0,0,1),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,1),j2,get_res_complex_eval(5,0,1),ptg))

                NLO_dt_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,3)
                j1=[get_res_complex_eval(2,0,1),get_res_complex_eval(2,1,1),get_res_complex_eval(2,2,1),get_res_complex_eval(2,3,1)]
                j2=[get_res_complex_eval(4,0,1),get_res_complex_eval(4,1,1),get_res_complex_eval(4,2,1),get_res_complex_eval(4,3,1)]
                ptg=[get_res_complex_eval(6,0,1),get_res_complex_eval(6,1,1),get_res_complex_eval(6,2,1),get_res_complex_eval(6,3,1)]
                resf.append(cut_res(get_res_complex_eval(0,0,1),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,1),j2,get_res_complex_eval(5,0,1),ptg))

                NLO_dt_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,4)
                j1=[get_res_complex_eval(2,0,1),get_res_complex_eval(2,1,1),get_res_complex_eval(2,2,1),get_res_complex_eval(2,3,1)]
                j2=[get_res_complex_eval(4,0,1),get_res_complex_eval(4,1,1),get_res_complex_eval(4,2,1),get_res_complex_eval(4,3,1)]
                ptg=[get_res_complex_eval(6,0,1),get_res_complex_eval(6,1,1),get_res_complex_eval(6,2,1),get_res_complex_eval(6,3,1)]
                resf.append(cut_res(get_res_complex_eval(0,0,1),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,1),j2,get_res_complex_eval(5,0,1),ptg))

        if prec=='f128':
            if phase=='real':
                NLO_dt_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
                j1=[get_res_complex_eval_f128(2,0,0),get_res_complex_eval_f128(2,1,0),get_res_complex_eval_f128(2,2,0),get_res_complex_eval_f128(2,3,0)]
                j2=[get_res_complex_eval_f128(4,0,0),get_res_complex_eval_f128(4,1,0),get_res_complex_eval_f128(4,2,0),get_res_complex_eval_f128(4,3,0)]
                ptg=[get_res_complex_eval_f128(6,0,0),get_res_complex_eval_f128(6,1,0),get_res_complex_eval_f128(6,2,0),get_res_complex_eval_f128(6,3,0)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,0),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,0),j2,get_res_complex_eval_f128(5,0,0),ptg))

                NLO_dt_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
                j1=[get_res_complex_eval_f128(2,0,0),get_res_complex_eval_f128(2,1,0),get_res_complex_eval_f128(2,2,0),get_res_complex_eval_f128(2,3,0)]
                j2=[get_res_complex_eval_f128(4,0,0),get_res_complex_eval_f128(4,1,0),get_res_complex_eval_f128(4,2,0),get_res_complex_eval_f128(4,3,0)]
                ptg=[get_res_complex_eval_f128(6,0,0),get_res_complex_eval_f128(6,1,0),get_res_complex_eval_f128(6,2,0),get_res_complex_eval_f128(6,3,0)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,0),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,0),j2,get_res_complex_eval_f128(5,0,0),ptg))

                NLO_dt_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,3)
                j1=[get_res_complex_eval_f128(2,0,0),get_res_complex_eval_f128(2,1,0),get_res_complex_eval_f128(2,2,0),get_res_complex_eval_f128(2,3,0)]
                j2=[get_res_complex_eval_f128(4,0,0),get_res_complex_eval_f128(4,1,0),get_res_complex_eval_f128(4,2,0),get_res_complex_eval_f128(4,3,0)]
                ptg=[get_res_complex_eval_f128(6,0,0),get_res_complex_eval_f128(6,1,0),get_res_complex_eval_f128(6,2,0),get_res_complex_eval_f128(6,3,0)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,0),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,0),j2,get_res_complex_eval_f128(5,0,0),ptg))

                NLO_dt_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,4)
                j1=[get_res_complex_eval_f128(2,0,0),get_res_complex_eval_f128(2,1,0),get_res_complex_eval_f128(2,2,0),get_res_complex_eval_f128(2,3,0)]
                j2=[get_res_complex_eval_f128(4,0,0),get_res_complex_eval_f128(4,1,0),get_res_complex_eval_f128(4,2,0),get_res_complex_eval_f128(4,3,0)]
                ptg=[get_res_complex_eval_f128(6,0,0),get_res_complex_eval_f128(6,1,0),get_res_complex_eval_f128(6,2,0),get_res_complex_eval_f128(6,3,0)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,0),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,0),j2,get_res_complex_eval_f128(5,0,0),ptg))

            if phase=='imag':
                NLO_dt_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
                j1=[get_res_complex_eval_f128(2,0,1),get_res_complex_eval_f128(2,1,1),get_res_complex_eval_f128(2,2,1),get_res_complex_eval_f128(2,3,1)]
                j2=[get_res_complex_eval_f128(4,0,1),get_res_complex_eval_f128(4,1,1),get_res_complex_eval_f128(4,2,1),get_res_complex_eval_f128(4,3,1)]
                ptg=[get_res_complex_eval_f128(6,0,1),get_res_complex_eval_f128(6,1,1),get_res_complex_eval_f128(6,2,1),get_res_complex_eval_f128(6,3,1)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,1),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,1),j2,get_res_complex_eval_f128(5,0,1),ptg))

                NLO_dt_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
                j1=[get_res_complex_eval_f128(2,0,1),get_res_complex_eval_f128(2,1,1),get_res_complex_eval_f128(2,2,1),get_res_complex_eval_f128(2,3,1)]
                j2=[get_res_complex_eval_f128(4,0,1),get_res_complex_eval_f128(4,1,1),get_res_complex_eval_f128(4,2,1),get_res_complex_eval_f128(4,3,1)]
                ptg=[get_res_complex_eval_f128(6,0,1),get_res_complex_eval_f128(6,1,1),get_res_complex_eval_f128(6,2,1),get_res_complex_eval_f128(6,3,1)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,1),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,1),j2,get_res_complex_eval_f128(5,0,1),ptg))

                NLO_dt_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,3)
                j1=[get_res_complex_eval_f128(2,0,1),get_res_complex_eval_f128(2,1,1),get_res_complex_eval_f128(2,2,1),get_res_complex_eval_f128(2,3,1)]
                j2=[get_res_complex_eval_f128(4,0,1),get_res_complex_eval_f128(4,1,1),get_res_complex_eval_f128(4,2,1),get_res_complex_eval_f128(4,3,1)]
                ptg=[get_res_complex_eval_f128(6,0,1),get_res_complex_eval_f128(6,1,1),get_res_complex_eval_f128(6,2,1),get_res_complex_eval_f128(6,3,1)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,1),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,1),j2,get_res_complex_eval_f128(5,0,1),ptg))

                NLO_dt_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,4)
                j1=[get_res_complex_eval_f128(2,0,1),get_res_complex_eval_f128(2,1,1),get_res_complex_eval_f128(2,2,1),get_res_complex_eval_f128(2,3,1)]
                j2=[get_res_complex_eval_f128(4,0,1),get_res_complex_eval_f128(4,1,1),get_res_complex_eval_f128(4,2,1),get_res_complex_eval_f128(4,3,1)]
                ptg=[get_res_complex_eval_f128(6,0,1),get_res_complex_eval_f128(6,1,1),get_res_complex_eval_f128(6,2,1),get_res_complex_eval_f128(6,3,1)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,1),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,1),j2,get_res_complex_eval_f128(5,0,1),ptg))


    if tag=='se':
        if prec=='f64':
            if phase=='real':

                NLO_se_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
                j1=[get_res_complex_eval(2,0,0),get_res_complex_eval(2,1,0),get_res_complex_eval(2,2,0),get_res_complex_eval(2,3,0)]
                j2=[get_res_complex_eval(4,0,0),get_res_complex_eval(4,1,0),get_res_complex_eval(4,2,0),get_res_complex_eval(4,3,0)]
                ptg=[get_res_complex_eval(6,0,0),get_res_complex_eval(6,1,0),get_res_complex_eval(6,2,0),get_res_complex_eval(6,3,0)]
                resf.append(cut_res(get_res_complex_eval(0,0,0),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,0),j2,get_res_complex_eval(5,0,0),ptg))

                NLO_se_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
                j1=[get_res_complex_eval(2,0,0),get_res_complex_eval(2,1,0),get_res_complex_eval(2,2,0),get_res_complex_eval(2,3,0)]
                j2=[get_res_complex_eval(4,0,0),get_res_complex_eval(4,1,0),get_res_complex_eval(4,2,0),get_res_complex_eval(4,3,0)]
                ptg=[get_res_complex_eval(6,0,0),get_res_complex_eval(6,1,0),get_res_complex_eval(6,2,0),get_res_complex_eval(6,3,0)]
                resf.append(cut_res(get_res_complex_eval(0,0,0),get_res_complex_eval(1,0,0),j1,get_res_complex_eval(3,0,0),j2,get_res_complex_eval(5,0,0),ptg))

            if phase=='imag':
                NLO_se_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
                j1=[get_res_complex_eval(2,0,1),get_res_complex_eval(2,1,1),get_res_complex_eval(2,2,1),get_res_complex_eval(2,3,1)]
                j2=[get_res_complex_eval(4,0,1),get_res_complex_eval(4,1,1),get_res_complex_eval(4,2,1),get_res_complex_eval(4,3,1)]
                ptg=[get_res_complex_eval(6,0,1),get_res_complex_eval(6,1,1),get_res_complex_eval(6,2,1),get_res_complex_eval(6,3,1)]
                resf.append(cut_res(get_res_complex_eval(0,0,1),get_res_complex_eval(1,0,1),j1,get_res_complex_eval(3,0,1),j2,get_res_complex_eval(5,0,1),ptg))

                NLO_se_DrellYan_eval(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
                j1=[get_res_complex_eval(2,0,1),get_res_complex_eval(2,1,1),get_res_complex_eval(2,2,1),get_res_complex_eval(2,3,1)]
                j2=[get_res_complex_eval(4,0,1),get_res_complex_eval(4,1,1),get_res_complex_eval(4,2,1),get_res_complex_eval(4,3,1)]
                ptg=[get_res_complex_eval(6,0,1),get_res_complex_eval(6,1,1),get_res_complex_eval(6,2,1),get_res_complex_eval(6,3,1)]
                resf.append(cut_res(get_res_complex_eval(0,0,1),get_res_complex_eval(1,0,1),j1,get_res_complex_eval(3,0,1),j2,get_res_complex_eval(5,0,1),ptg))

        if prec=='f128':
            if phase=='real':

                NLO_se_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
                j1=[get_res_complex_eval_f128(2,0,0),get_res_complex_eval_f128(2,1,0),get_res_complex_eval_f128(2,2,0),get_res_complex_eval_f128(2,3,0)]
                j2=[get_res_complex_eval_f128(4,0,0),get_res_complex_eval_f128(4,1,0),get_res_complex_eval_f128(4,2,0),get_res_complex_eval_f128(4,3,0)]
                ptg=[get_res_complex_eval_f128(6,0,0),get_res_complex_eval_f128(6,1,0),get_res_complex_eval_f128(6,2,0),get_res_complex_eval_f128(6,3,0)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,0),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,0),j2,get_res_complex_eval_f128(5,0,0),ptg))

                NLO_se_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
                j1=[get_res_complex_eval_f128(2,0,0),get_res_complex_eval_f128(2,1,0),get_res_complex_eval_f128(2,2,0),get_res_complex_eval_f128(2,3,0)]
                j2=[get_res_complex_eval_f128(4,0,0),get_res_complex_eval_f128(4,1,0),get_res_complex_eval_f128(4,2,0),get_res_complex_eval_f128(4,3,0)]
                ptg=[get_res_complex_eval_f128(6,0,0),get_res_complex_eval_f128(6,1,0),get_res_complex_eval_f128(6,2,0),get_res_complex_eval_f128(6,3,0)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,0),get_res_complex_eval_f128(1,0,0),j1,get_res_complex_eval_f128(3,0,0),j2,get_res_complex_eval_f128(5,0,0),ptg))

            if phase=='imag':
                NLO_se_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,1)
                j1=[get_res_complex_eval_f128(2,0,1),get_res_complex_eval_f128(2,1,1),get_res_complex_eval_f128(2,2,1),get_res_complex_eval_f128(2,3,1)]
                j2=[get_res_complex_eval_f128(4,0,1),get_res_complex_eval_f128(4,1,1),get_res_complex_eval_f128(4,2,1),get_res_complex_eval_f128(4,3,1)]
                ptg=[get_res_complex_eval_f128(6,0,1),get_res_complex_eval_f128(6,1,1),get_res_complex_eval_f128(6,2,1),get_res_complex_eval_f128(6,3,1)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,1),get_res_complex_eval_f128(1,0,1),j1,get_res_complex_eval_f128(3,0,1),j2,get_res_complex_eval_f128(5,0,1),ptg))

                NLO_se_DrellYan_eval_f128(k[0],k[1],k[2],l[0],l[1],l[2],q[0],q[1],q[2],a,2)
                j1=[get_res_complex_eval_f128(2,0,1),get_res_complex_eval_f128(2,1,1),get_res_complex_eval_f128(2,2,1),get_res_complex_eval_f128(2,3,1)]
                j2=[get_res_complex_eval_f128(4,0,1),get_res_complex_eval_f128(4,1,1),get_res_complex_eval_f128(4,2,1),get_res_complex_eval_f128(4,3,1)]
                ptg=[get_res_complex_eval_f128(6,0,1),get_res_complex_eval_f128(6,1,1),get_res_complex_eval_f128(6,2,1),get_res_complex_eval_f128(6,3,1)]
                resf.append(cut_res(get_res_complex_eval_f128(0,0,1),get_res_complex_eval_f128(1,0,1),j1,get_res_complex_eval_f128(3,0,1),j2,get_res_complex_eval_f128(5,0,1),ptg))

    return resf


myres=NLO([1,2,3],[1.7,0.8,-1.5],[-0.93,-2.4,-1],0.2,'dt',phase='real')
myres=NLO([1,2,3],[1.7,0.8,-1.5],[-0.93,-2.4,-1],0.2,'dt',phase='imag')
myres=NLO([1,2,3],[1.7,0.8,-1.5],[-0.93,-2.4,-1],0.2,'st',prec='f128',phase='real')
