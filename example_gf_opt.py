""" Example of using the Green's Function approach to optimize parameters

Observation data as used in Ungermann et al. (2017, JGR Oceans)
"""
from __future__ import print_function

import numpy as np
from scipy.io import netcdf_file
import gf_opt as gf

class myData: 
    """ container for model and observational data to use in GF approach
    " data - list of arrays with data
    " datanames - dictionary {"dataname" : index}
    " sigma - measurement uncertainty
    " sigmanames - dictionary, differing from datanames for velocities
    " mask_conc - sea ice mask for observations
    """ 
    def __init__(self,kind=''):
        #
        if kind=='Obs':
            self.data=[]
            self.datanames=dict()
            self.sigma=[]
            self.sigmanames=dict()
            self.mask_conc=None
        elif kind=='Model':
            self.data=[]
            self.datanames=dict()
        else: 
            print(" Unknown kind {0}! ".format(kind))
            kind=None
        
class Dataflags: 
    """ collection of flags to distinguish between runs and used data
    """
    def __init__(self):
        # choice of used datasets 
        self.use_osiconc=1
        self.use_icesat=1
        self.use_osidrift=1
        self.use_kimura=1 
        # choice of used simulations  
        self.flag_strgf_Pmean=0 # Hibler Strength
        self.flag_strgf_Pdist=0 # Full Rothrock Strength

def fill_info(dataflags): 
    """ sets necessary run informations
    provides paths for baseline run, perturbed runs and possibly optimised run
    provides perturbations, original values and names for parameters
    flags decide which runs to use in optimsiation
    """ 
    toolpath='/home/ollie/mungerma/costfunction_NAOSIM/toolkit/'
    pert_path=[]
    pert=[]
    orig_parms=[]
    parm_names=[]
    optim_path=None
    basel_path=None
    if dataflags.flag_strgf_Pdist: 
        if basel_path is None:
            basel_path='/work/ollie/mungerma/arctic_output/str_gf/Pdist_bl/results/'
        else:
            print('WARNING: Using perturbation experiments with different' +
                  'baselines in the same optimization!')
        optim_path='/work/ollie/mungerma/arctic_output/str_gf/Pdist_opt1/results/'
        pert_path+=[
               '/work/ollie/mungerma/arctic_output/str_gf/Pdist_cf/results/',
               '/work/ollie/mungerma/arctic_output/str_gf/Pdist_G/results/',
               '/work/ollie/mungerma/arctic_output/str_gf/Pdist_H/results/',
                ]
        pert += [
                (0.16 - 0.12),      # G* (perturbed - baseline)
                (40. - 25.),        # H* (perturbed - baseline)
                (12. - 14.),        # C_f (perturbed - baseline)
                ]
        orig_parms +=[
                0.12,
                25.0,
                14.,
                ]
        parm_names += [
                'G*',
                'H*',
                'C_f',
                ] 
        
    if dataflags.flag_strgf_Pmean: 
        if basel_path is None:
            basel_path='/work/ollie/mungerma/arctic_output/str_gf/Pmean_bl/results/'
        else:
            print('WARNING: Using perturbation experiments with different' +
                  'baselines in the same optimization!')
        optim_path='/work/ollie/mungerma/arctic_output/str_gf/Pmean_opt1/results/'
        pert_path+=[
               '/work/ollie/mungerma/arctic_output/str_gf/Pmean_C/results/',
               '/work/ollie/mungerma/arctic_output/str_gf/Pmean_P/results/',
                ]
        pert += [
                (0.16 - 0.12),       # C~ (perturbed - baseline)
                (18250. - 14250.),   # P* (perturbed - baseline)
                ]
        orig_parms +=[
                0.12,
                14250.0,
                ]
        parm_names += [
                'C~',
                'P*',
                ] 

    if not optim_path:
        optim_path='/work/ollie/mungerma/arctic_output/str_gf/Pdist_bl/results/'
    pert=np.array(pert)
    orig_parms=np.array(orig_parms)
    parm_names=np.array(parm_names)

    return basel_path, pert_path, optim_path, toolpath, orig_parms, pert, parm_names


def readObs(basel_path, toolpath, dataflags):
    """
    read observational data ''obs''
    includes lists 
        obs.data
        obs.datanames
        obs.sigma
        obs.sigmanames
        obs.mask_conc
    """
    obs=myData(kind='Obs')
    i=0
    j=0
    if dataflags.use_osiconc:
        osisaf_file=netcdf_file(basel_path+'osisaf-miss.nc','r')
        osi_conc=osisaf_file.variables['ice_conc'][:]
        osi_err=osisaf_file.variables['standard_error'][:]
        osisaf_file.close()
        mask_file=netcdf_file(basel_path+'costfct-mask-2d.nc')
        mask_conc=mask_file.variables['SIarea'][:]
        mask_file.close()
        obs.data.append(osi_conc)
        obs.sigma.append(osi_err)
        obs.mask_conc=mask_conc
        obs.datanames['osiconc']=i
        i+=1
        obs.sigmanames['osiconc']=j
        j+=1
    # read Thickness
    if dataflags.use_icesat:
        icesat_file=netcdf_file(basel_path+'icesat-mask.nc','r')
        icesat_thick=icesat_file.variables['var1'][:]
        icesat_file.close()
        icesaterr_file=netcdf_file(toolpath+'icesat-mask-newerror.nc','r')
        icesat_err=icesaterr_file.variables['new_err'][:].squeeze()
        icesaterr_file.close()
        obs.data.append(icesat_thick)
        obs.sigma.append(icesat_err)
        obs.datanames['icesat']=i
        i+=1
        obs.sigmanames['icesat']=j
        j+=1
    # read OSISAF drift
    if dataflags.use_osidrift:
        osidrift_file=netcdf_file(toolpath+'ice_drift_NAOSIM-rot.nc','r')
        osidrift_u=osidrift_file.variables['urot'][:]
        osidrift_v=osidrift_file.variables['vrot'][:]
        osidrift_file.close()
        osidrift_errfile=netcdf_file(toolpath+'uncertainty_NAOSIM_0.25.nc')
        osidrift_err=osidrift_errfile.variables['uncertainty'][:]
        osidrift_errfile.close()
        obs.data.append(osidrift_u)
        obs.datanames['osidrift_u']=i
        i+=1
        obs.data.append(osidrift_v)
        obs.datanames['osidrift_v']=i
        i+=1
        obs.sigma.append(osidrift_err)
        obs.sigmanames['osidrift']=j
        j+=1
    # read Kimura drift
    if dataflags.use_kimura:
        kimura_file=netcdf_file(toolpath+'kimura-summer-rot.nc','r')
        kimura_u=kimura_file.variables['urot'][:]*0.01
        kimura_v=kimura_file.variables['vrot'][:]*0.01
        kimura_file.close()
        kimura_errfile=netcdf_file(toolpath+'kimura-summer-naosim.0.25.nc')
        kimura_err=kimura_errfile.variables['uncertainty'][:]*0.01
        kimura_errfile.close()
        obs.data.append(kimura_u)
        obs.datanames['kimura_u']=i
        i+=1
        obs.data.append(kimura_v)
        obs.datanames['kimura_v']=i
        i+=1
        obs.sigma.append(kimura_err)
        obs.sigmanames['kimura']=j
        j+=1
    return obs

def readModel(model_path, dataflags):
    """ Read model data ''model''
    includes 
        base.data       (list)
        base.datanames  (dict)
    """
    model=myData(kind='Model')
    i=0
    if dataflags.use_osiconc:
        conc_file=netcdf_file(model_path+'ice_concn-miss.nc')
        base_conc=conc_file.variables['SIarea'][:].squeeze()
        conc_file.close()
        model.data.append(base_conc)
        model.datanames['osiconc']=i
        i+=1
    if dataflags.use_icesat:
        thick_file=netcdf_file(model_path+'ice_thick_icesat-cut7.nc','r')
        base_thick=thick_file.variables['SIheff'][:].squeeze()
        thick_file.close()
        model.data.append(base_thick)
        model.datanames['icesat']=i
        i+=1
    if dataflags.use_osidrift:
        uvel_file=netcdf_file(model_path+'u_icevelocity-osisaf-rot.nc','r')
        base_uvel=uvel_file.variables['urot'][:].squeeze()
        uvel_file.close()
        vvel_file=netcdf_file(model_path+'v_icevelocity-osisaf-rot.nc','r')
        base_vvel=vvel_file.variables['vrot'][:].squeeze()
        vvel_file.close()
        model.data.append(base_uvel)
        model.datanames['osidrift_u']=i
        i+=1
        model.data.append(base_vvel)
        model.datanames['osidrift_v']=i
        i+=1
    if dataflags.use_kimura:
        uvel_file=netcdf_file(model_path+'u_icevelocity-kimura-rot.nc','r')
        base_ukim=uvel_file.variables['urot'][:].squeeze()
        uvel_file.close()
        vvel_file=netcdf_file(model_path+'v_icevelocity-kimura-rot.nc','r')
        base_vkim=vvel_file.variables['vrot'][:].squeeze()
        vvel_file.close()
        model.data.append(base_ukim)
        model.datanames['kimura_u']=i
        i+=1
        model.data.append(base_vkim)
        model.datanames['kimura_v']=i
        i+=1
    return model

def getData(dataflags):
    """ read data specified by used dataflags
    " returns:
    " x - array of baseline model values
    " y - array of obsevations
    " x_pert - list of arrays of perturbed model values
    " pertvals - magnitude of perturbation of parameters
    " nr_sigma - startin indices for each subset of used data
    """ 
    basel_path,pert_path,optim_path,toolpath,orig_parms,pertvals,parm_names = fill_info(dataflags)
    # read data
    obs=readObs(basel_path,toolpath,dataflags)
    base=readModel(basel_path,dataflags)
    pert=[]
    npert=0
    for ppath in pert_path:
        pert.append(readModel(ppath,dataflags))
        npert+=1
    # write common masks for comparison
    if dataflags.use_osiconc:
        areamask=np.logical_and(obs.mask_conc==1.,base.data[base.datanames['osiconc']]<1.e10)
        areamask=np.logical_and(areamask,obs.sigma[obs.sigmanames['osiconc']]<1.e10)
        for i in np.arange(npert):
            areamask=np.logical_and(areamask,pert[i].data[pert[i].datanames['osiconc']]<1.e10)
    if dataflags.use_icesat:
        thickmask=np.logical_and(obs.data[obs.datanames['icesat']]>-1.,
            base.data[base.datanames['icesat']]<1.e10)
        thickmask=np.logical_and(thickmask,obs.sigma[obs.sigmanames['icesat']]<10.)
        for i in np.arange(npert):
            thickmask=np.logical_and(thickmask,pert[i].data[pert[i].datanames['icesat']]<1.e10)
    if dataflags.use_osidrift:
        osidrift_mask= np.logical_and(base.data[base.datanames['osidrift_u']]>-100.,
            base.data[base.datanames['osidrift_v']]>-100.)
        for i in np.arange(npert):
            osidrift_mask=np.logical_and(osidrift_mask,
                    np.logical_and(pert[i].data[pert[i].datanames['osidrift_u']]>-100.,
                        pert[i].data[pert[i].datanames['osidrift_v']]>-100))
        osidrift_mask=np.logical_and(osidrift_mask,
                obs.sigma[obs.sigmanames['osidrift']]>-1.e2)
    if dataflags.use_kimura:
        kimura_mask= np.logical_and(base.data[base.datanames['kimura_u']]>-100.,
            base.data[base.datanames['kimura_v']]>-100.)
        for i in np.arange(npert):
            kimura_mask=np.logical_and(kimura_mask,
                    np.logical_and(pert[i].data[pert[i].datanames['kimura_u']]>-100.,
                        pert[i].data[pert[i].datanames['kimura_v']]>-100.))
        kimura_mask=np.logical_and(kimura_mask,
                obs.sigma[obs.sigmanames['kimura']]>-1.e2)
        kimura_mask=np.logical_and(kimura_mask,
                np.logical_and(obs.data[obs.datanames['kimura_u']]>-1.e5,
                    obs.data[obs.datanames['kimura_v']]>-1.e5))
        kimura_mask=np.logical_and(kimura_mask,
                np.logical_and(obs.data[obs.datanames['kimura_u']]!=0.,
                    obs.data[obs.datanames['kimura_v']]!=0.))
    # write vectors y,sigma for used observations
    y=np.array([],dtype='>f4')
    sigma=np.array([],dtype='>f4')
    nr_sigma=[]
    if dataflags.use_osiconc:
        y_aux=obs.data[obs.datanames['osiconc']][areamask]
        y=np.concatenate([y,y_aux])
        sigma_aux=obs.sigma[obs.sigmanames['osiconc']][areamask]*np.sqrt(len(y_aux))
        sigma=np.concatenate([sigma,sigma_aux])
        nr_sigma.append(np.size(sigma_aux))
    if dataflags.use_icesat:
        y_aux=obs.data[obs.datanames['icesat']][thickmask]
        y=np.concatenate([y,y_aux])
        sigma_aux=obs.sigma[obs.sigmanames['icesat']][thickmask]*np.sqrt(len(y_aux))
        sigma=np.concatenate([sigma,sigma_aux])
        nr_sigma.append(np.size(sigma_aux))
    if dataflags.use_osidrift:
        y_aux=np.sqrt(np.square(obs.data[obs.datanames['osidrift_u']][osidrift_mask])+
            np.square(obs.data[obs.datanames['osidrift_v']][osidrift_mask]))
        y=np.concatenate([y,y_aux])
        sigma_aux=obs.sigma[obs.sigmanames['osidrift']][osidrift_mask]*np.sqrt(len(y_aux))
        sigma=np.concatenate([sigma,sigma_aux])
        nr_sigma.append(np.size(sigma_aux))
    if dataflags.use_kimura:
        y_aux=np.sqrt(np.square(obs.data[obs.datanames['kimura_u']][kimura_mask])+
            np.square(obs.data[obs.datanames['kimura_v']][kimura_mask]))
        y=np.concatenate([y,y_aux])
        sigma_aux=obs.sigma[obs.sigmanames['kimura']][kimura_mask]*np.sqrt(len(y_aux))
        sigma=np.concatenate([sigma,sigma_aux])
        nr_sigma.append(np.size(sigma_aux))
    nr_sigma=np.concatenate(([0],np.cumsum(nr_sigma)))
    # write vector x for baseline model data
    x=np.array([],dtype='>f4')
    if dataflags.use_osiconc:
        x_aux=base.data[base.datanames['osiconc']][areamask]
        x=np.concatenate([x,x_aux])
    if dataflags.use_icesat:
        x_aux=base.data[base.datanames['icesat']][thickmask]
        x=np.concatenate([x,x_aux])
    if dataflags.use_osidrift:
        x_aux=np.sqrt(np.square(base.data[base.datanames['osidrift_u']][osidrift_mask])+
            np.square(base.data[base.datanames['osidrift_v']][osidrift_mask]))
        x=np.concatenate([x,x_aux])
    if dataflags.use_kimura:
        x_aux=np.sqrt(np.square(base.data[base.datanames['kimura_u']][kimura_mask])+
            np.square(base.data[base.datanames['kimura_v']][kimura_mask]))
        x=np.concatenate([x,x_aux])
    # write list of vectors x_pert for perturbed model data
    x_pert=[]
    for i in range(npert):
        x_pert_aux=np.array([],dtype='>f4')
        if dataflags.use_osiconc:
            x_aux=pert[i].data[pert[i].datanames['osiconc']][areamask]
            x_pert_aux=np.concatenate([x_pert_aux,x_aux])
        if dataflags.use_icesat:
            x_aux=pert[i].data[pert[i].datanames['icesat']][thickmask]
            x_pert_aux=np.concatenate([x_pert_aux,x_aux])
        if dataflags.use_osidrift:
            x_aux=np.sqrt(np.square(pert[i].data[pert[i].datanames['osidrift_u']][osidrift_mask])+
                np.square(pert[i].data[pert[i].datanames['osidrift_v']][osidrift_mask]))
            x_pert_aux=np.concatenate([x_pert_aux,x_aux])
        if dataflags.use_kimura:
            x_aux=np.sqrt(np.square(pert[i].data[pert[i].datanames['kimura_u']][kimura_mask])+
                np.square(pert[i].data[pert[i].datanames['kimura_v']][kimura_mask]))
            x_pert_aux=np.concatenate([x_pert_aux,x_aux])
        x_pert.append(x_pert_aux)
    return x,y,x_pert,sigma,pertvals,nr_sigma

if __name__ == '__main__':
    ## set up choice of validation data and simulations
    dataflags=Dataflags()
    dataflags.flag_strgf_Pmean = 1
    dataflags.flag_strgf_Pdist = 0
    dataflags.use_osiconc=0 # produces MemoryError (19/03/12)
    dataflags.use_icesat=1
    dataflags.use_osidrift=1
    dataflags.use_kimura=1

    # read data
    basel_path, pert_path, optim_path, toolpath, orig_parms, pert, parm_names = fill_info(dataflags)
    x, y, x_pert, sigma, pertvals, nr_sigma = getData(dataflags)

    # caculate optimal parameters
    G = gf.calcDataKernel(x, x_pert, pertvals)
    eta_new = gf.optimiseGF(x, y, sigma, G)
    new_parms = orig_parms + eta_new

    # output
    gf.printParms(parm_names,orig_parms,pert,eta_new)
