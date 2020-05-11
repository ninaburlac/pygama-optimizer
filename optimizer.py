#!/usr/bin/env python3
import os, json
import argparse
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('clint.mpl')
from pprint import pprint
import scipy.signal as signal
import itertools

from pygama import DataSet
import pygama.utils as pu
import pygama.analysis.histograms as ph
import pygama.analysis.peak_fitting as pf
import h5py, sys
import pygama.io.lh5 as lh5

from pygama.dsp.ProcessingChain import ProcessingChain
from pygama.dsp.processors import *
from pygama.dsp.units import *


def main():
    """
    To get the best energy resolution, we want to explore the possible values
    of our DSP processor list, especially trap filter and RC decay constants.
    Inclusion of dedicated functions for the optimization of the ZAC filter.
    
    Modified by:
    V. D'Andrea
    """
    
    filters = ['trapE','trapEftp','zacE','cuspE','trapEcorr','trapEftpcorr','zacEcorr']    
    
    par = argparse.ArgumentParser(description="pygama dsp optimizer")
    arg, st, sf = par.add_argument, "store_true", "store_false"
    arg("-d", "--dir", nargs=1, action="store", help="analysis directory")
    arg("-ds", nargs='*', action="store", help="load runs for a DS")
    arg("-r", "--run", nargs=1, help="load a single run")
    arg("-g", "--grid", action=st, help="set DSP parameters to be varied")
    arg("-w", "--window", action=st, help="generate a small waveform file")
    arg("-p", "--process", action=st, help="run DSP processing")
    arg("-f", "--fit", action=st, help="fit outputs to peakshape function")
    arg("-t", "--plot", action=st, help="find optimal parameters & make plots")
    arg("-c", "--case", nargs=1, help="energy filter: 0=trapE, 1=trapEftp, 2=zacE, 3=cuspE, 4=trapEcorr, 5=trapEftpcorr,6=zacEcorr")
    arg("-cf", "--compare", nargs='*', help="compare energy filter")
    arg("-v", "--verbose", action=st, help="set verbose mode")
    args = vars(par.parse_args())
    
    local_dir = "."
    if args["dir"]: local_dir = args["dir"][0]
    run = args["run"][0]
    
    case = 0
    if args["case"]: case = int(args["case"][0])
    if case < len(filters):
        print(f"Run {filters[case]} optimization")
    else:
        print("Filter Case not valid")
        return

    #ds = pu.get_dataset_from_cmdline(args, f"{local_dir}/meta/runDB.json", f"{local_dir}/meta/calDB.json")
    with open(f"{local_dir}/meta/runDB.json") as frun:
        config = json.load(frun)
    #pprint(ds.paths)
    
    d_out = f"{local_dir}/run{run:0>4}"
    try: os.mkdir(d_out)
    except: pass
    
    f_tier1 = f"{d_out}/cage_optimizer_raw.h5"
    f_grid = f"{d_out}/{filters[case]}_optimizer_grid.h5"
    f_opt = f"{d_out}/{filters[case]}_optimizer_dsp.h5"
    
    # -- run routines --
    if args["grid"]: set_grid(f_grid,filters[case])
    
    # generate a small single-peak file w/ uncalibrated energy to reanalyze
    #if args["window"]: window_ds(ds, f_tier1)
    if args["window"]: window_ds(config, f_tier1)
    
    # create a file with DataFrames for each set of parameters
    if args["process"]: process_ds(f_grid, f_opt, f_tier1, d_out, filters[case])
    
    # fit all outputs to the peakshape function and find the best resolution
    if args["fit"]: get_fwhm(f_grid, f_opt, filters[case],verbose=args["verbose"])
    
    # show results>
    if args["plot"]: plot_fwhm(f_grid,f_opt,d_out,filters[case],verbose=args["verbose"])

    # compare fwhm results
    if args["compare"]:
        filter1, filter2 = filters[int(args["compare"][0])], filters[int(args["compare"][1])]
        try:
            filter3 = filters[int(args["compare"][2])]
        except:
            filter3 = 0;
        compare_fwhm(d_out,filter1,filter2,filter3,verbose=args["verbose"])
    
    
def set_grid(f_grid,efilter):
    """
    create grid with set of parameters
    """
    if 'corr' in efilter:
        print("Creation of grid for charge trapping correction")
        ct_consts = np.linspace(-0.05, 0., 6, dtype='float')
        lists = [ct_consts]
        prod = list(itertools.product(*lists))
        df = pd.DataFrame(prod, columns=['ct_const']) 
    elif 'trapE' in efilter:
        print(f"Creation of grid for {efilter} optimization")
        rises = np.linspace(1, 10, 10, dtype='float')
        flats = np.linspace(0.5, 4.5, 5, dtype='float')
        rcs = np.linspace(500, 500, 1, dtype='float')
        lists = [rises, flats, rcs]
        prod = list(itertools.product(*lists))
        df = pd.DataFrame(prod, columns=['rise','flat','rc']) 
    elif 'zacE' in efilter or 'cuspE' in efilter:
        print(f"Creation of grid for {efilter} optimization")
        sigmas = np.linspace(1, 50, 10, dtype='float')
        flats =  np.linspace(0.5, 4.5, 5, dtype='float')
        decays =  np.linspace(500, 500, 1, dtype='float')
        lists = [sigmas, flats, decays]
        prod = list(itertools.product(*lists))
        df = pd.DataFrame(prod, columns=['sigma', 'flat','decay'])     
        
    print(df)
    df.to_hdf(f_grid, key="pygama_optimization")
    print("Wrote grid file:", f_grid)
    

#def window_ds(ds, f_tier1):
def window_ds(config, f_tier1):
    """
    Take a DataSet and window it so that the output file only contains 
    events near the calibration peak at 2614.5 keV.
    """
    print("Creating windowed raw file:",f_tier1)
    f_win = h5py.File(f_tier1, 'w')

    #raw_dir = ds.config["raw_dir"]
    #geds = ds.config["daq_to_raw"]["ch_groups"]["g{ch:0>3d}"]["ch_range"]
    raw_dir = config["raw_dir"]
    geds = config["daq_to_raw"]["ch_groups"]["g{ch:0>3d}"]["ch_range"]
    cols = ['energy','baseline','ievt','numtraces','timestamp','wf_max','wf_std','waveform/values','waveform/dt']
    
    for ged in range(geds[0],geds[1]+1):
        ged = f"g{ged:0>3d}"
        count = 0
        for p, d, files in os.walk(raw_dir):
            for f in files:
                if (f.endswith(".lh5")) & ("calib" in f):
                    print("Opening raw file:",f)
                    f_raw = h5py.File(f"{raw_dir}/{f}",'r')
                    if count == 0:
                        #cdate, ctime = f.split('run')[-1].split('-')[1], f.split('run')[-1].split('-')[2]
                        dsets = [ f_raw[ged]['raw'][col][()]  for col in cols ]
                    else:
                        for i, col in enumerate(cols): dsets[i] = np.append(dsets[i],f_raw[ged]['raw'][col][()],axis=0)
                    count += 1
        
        # search for 2.6 MeV peak
        energies = dsets[0]
        maxe = np.amax(energies)
        h, b, v = ph.get_hist(energies, bins=3500, range=(maxe/4,maxe))
        #xp = b[np.where(h > h.max()*0.1)][-1]
        #h, b = h[np.where(b < xp-200)], b[np.where(b < xp-200)]
        bin_max = b[np.where(h == h.max())][0]
        min_ene = int(bin_max*0.95)
        max_ene = int(bin_max*1.05)
        hist, bins, var = ph.get_hist(energies, bins=500, range=(min_ene, max_ene))
        print(ged,"Raw energy max",maxe,"histogram max",h.max(),"at",bin_max )
        
        # windowing
        for i, col in enumerate(cols):
            dsets[i] = dsets[i][(energies>min_ene) & (energies<max_ene)]
            d_dt = f_win.create_dataset(ged+"/raw/"+col,dtype='f',data=dsets[i])
            d_dt.attrs['units'] = 'ns'
        
        f_win.attrs['datatype'] = 'table{cols}'
        print("Created datasets",ged+"/raw")
        
    f_win.close()        
    print("wrote file:", f_tier1)
    
    
def process_ds(f_grid, f_opt, f_tier1, d_out, efilter):
    """
    process the windowed raw file 'f_tier1' and create the DSP file 'f_opt'
    
    """
    print("Grid file:",f_grid)
    df_grid = pd.read_hdf(f_grid)
    
    if os.path.exists(f_opt):
        os.remove(f_opt)

    if 'corr' in efilter:
        bfilter = efilter.split('corr')[0]
        try:
            df_res = pd.read_hdf(f'{d_out}/{bfilter}_results.h5',key='results')
            print("Extraction of best parameters for", bfilter)
        except:
            print(bfilter,"not optimized")
            return
    
    # open raw file
    lh5_in = lh5.Store()
    #groups = lh5_in.ls(f_tier1, '*/raw')
    f = h5py.File(f_tier1,'r')
    #print("File info: ",f.keys())
    
    t_start = time.time()
    #for group in groups:
    for idx, ged in enumerate(f.keys()):
        if idx == 4:
            diff = time.time() - t_start
            tot = diff/5 * len(df_grid) / 60
            tot -= diff / 60
            print(f"Estimated remaining time: {tot:.2f} mins")
        
        print("Detector:",ged)
        #data = lh5_in.read_object(group, f_tier1)
        data =  f[ged]['raw']
        
        #wf_in = data['waveform']['values'].nda
        #dt = data['waveform']['dt'].nda[0] * unit_parser.parse_unit(data['waveform']['dt'].attrs['units'])
        wf_in = data['waveform']['values'][()]
        dt = data['waveform']['dt'][0] * unit_parser.parse_unit(data['waveform']['dt'].attrs['units'])
        bl_in = data['baseline'][()] #flashcam baseline values
        
        # Set up DSP processing chain -- very minimal
        block = 8 #waveforms to process simultaneously
        proc = ProcessingChain(block_width=block, clock_unit=dt, verbosity=False)
        proc.add_input_buffer("wf", wf_in, dtype='float32')
        proc.add_input_buffer("bl", bl_in, dtype='float32')
        
        wsize = wf_in.shape[1]
        dt0 = data['waveform']['dt'][0]*0.001
        
        #proc.add_processor(mean_stdev, "wf[0:1000]", "bl", "bl_sig")
        proc.add_processor(np.subtract, "wf", "bl", "wf_blsub")
        for i, row in df_grid.iterrows():
            if 'corr' in efilter: ct_const = row
            if 'trapE' in efilter:
                if 'corr' in efilter: rise, flat, rc = float(df_res['rise'][idx]), float(df_res['flat'][idx]), float(df_res['rc'][idx])
                else: rise, flat, rc = row
                proc.add_processor(pole_zero, "wf_blsub", rc*us, "wf_pz")
                proc.add_processor(trap_norm, "wf_pz", rise*us, flat*us, f"wf_trap_{i}")
                proc.add_processor(asymTrapFilter, "wf_pz", 0.05*us, 4*us, 4*us, "wf_atrap")
                proc.add_processor(time_point_thresh, "wf_pz", 0, "tp_0")
                proc.add_processor(np.amax, f"wf_trap_{i}", 1, f"trapE_{i}", signature='(n),()->()', types=['fi->f'])
                proc.add_processor(fixed_time_pickoff, f"wf_trap_{i}", f"tp_0+({rise*us}+{flat*us})", f"trapEftp_{i}")
            if 'zacE' in efilter:
                if 'corr' in efilter: sigma, flat, decay = float(df_res['sigma'][idx]), float(df_res['flat'][idx]), float(df_res['decay'][idx])
                else: sigma, flat, decay = row
                proc.add_processor(zac_filter(wsize, sigma/dt0, flat/dt0, decay/dt0),"wf", f"wf_zac_{i}(101, f)")
                proc.add_processor(np.amax, f"wf_zac_{i}", 1, f"zacE_{i}", signature='(n),()->()', types=['fi->f'])
            if 'cuspE' in efilter:
                if 'corr' in efilter: sigma, flat, decay = float(df_res['sigma'][idx]), float(df_res['flat'][idx]), float(df_res['decay'][idx])
                else: sigma, flat, decay = row
                proc.add_processor(cusp_filter(wsize, sigma/dt0, flat/dt0, decay/dt0),"wf_blsub", f"wf_cusp_{i}(101, f)")
                proc.add_processor(np.amax, f"wf_cusp_{i}", 1, f"cuspE_{i}", signature='(n),()->()', types=['fi->f'])
            if 'corr' in efilter:
                proc.add_processor(trap_pickoff, "wf_pz", 1.5*us, 0, "tp_0", "ct_corr")
                #proc.add_processor(trap_pickoff, "wf_pz", rise*us, flat*us, "tp_0", "ct_corr")
                proc.add_processor(np.multiply, ct_const, "ct_corr", f"ct_corr_cal_{i}")
                proc.add_processor(np.add, f"ct_corr_cal_{i}", f"{bfilter}_{i}", f"{efilter}_{i}")
                
        # Set up the LH5 output
        lh5_out = lh5.Table(size=proc._buffer_len)
        for i, row in df_grid.iterrows():
            lh5_out.add_field(f"{efilter}_{i}", lh5.Array(proc.get_output_buffer(f"{efilter}_{i}"), attrs={"units":"ADC"}))
        
        print("Processing:\n",proc)
        proc.execute()
        
        #groupname = group[:group.rfind('/')+1]+"data"
        #groupname = df_key+"/"+group+"/data"
        groupname = ged+"/data"
        print("Writing to: " + f_opt + "/" + groupname)
        lh5_in.write_object(lh5_out, groupname, f_opt)
        print("")
    
    #list the datasets of the output file
    data_opt = lh5_in.ls(f_opt)
    #data_opt_0 = lh5_in.ls(f_opt,'opt_0/*')
    data_opt_0 = lh5_in.ls(f_opt,'g024/data/*')
    diff = time.time() - t_start
    print(f"Time to process: {diff:.2f} s")
        

def get_fwhm(f_grid, f_opt, efilter, verbose=False):
    """
    this code fits the 2.6 MeV peak using the gauss+step function
    and writes new columns to the df_grid "fwhm", "fwhmerr"
    """
    print("Grid file:",f_grid)
    print("DSP file:",f_opt)
    df_grid = pd.read_hdf(f_grid)

    f = h5py.File(f_opt,'r')
    for ged in f.keys():
        print("Detector:",ged)
        data =  f[ged]['data']
        
        # declare some new columns for df_grid
        cols = [f"fwhm_{ged}", f"fwhmerr_{ged}", f"rchi2_{ged}"]
        for col in cols:
            df_grid[col] = np.nan
        
        for i, row in df_grid.iterrows():
            try:
                energies = data[f"{efilter}_{i}"][()]
                mean = np.mean(energies)
                bins = 12000
                hE, xE, vE = ph.get_hist(energies,bins,(mean/2,mean*2))
            except:
                print("Energy not find in",ged,"and entry",i)
            
            # set histogram centered and symmetric on the peak
            try:
                mu = xE[np.argmax(hE)]
                imax = np.argmax(hE)
                hmax = hE[imax]
                idx = np.where(hE > hmax/2) # fwhm
                ilo, ihi = idx[0][0], idx[0][-1]
                sig = (xE[ihi] - xE[ilo]) / 2.355
                idx = np.where(((xE-mu) > -8 * sig) & ((xE-mu) < 8 * sig))
                idx0 = np.where(((xE-mu) > -4.5 * sig) & ((xE-mu) < 4.5 * sig))
                ilo, ihi = idx[0][0], idx[0][-1]
                ilo0, ihi0 = idx0[0][0], idx0[0][-1]
                xE, hE, vE = xE[ilo:ihi+1], hE[ilo:ihi], vE[ilo:ihi]
            except:
                continue
            
            # set initial guesses for the peakshape function
            hstep = 0
            tau = np.mean(hE[:10])
            bg0 = 1
            x0 = [hmax, mu, sig, bg0, hstep]
            
            try:
                xF, xF_cov = pf.fit_hist(pf.gauss_step, hE, xE, var=vE, guess=x0)
                xF_err = np.sqrt(np.diag(xF_cov))
                
                # goodness of fit
                chisq = []
                for j, h in enumerate(hE):
                    model = pf.gauss_step(xE[j], *xF)
                    diff = (model - h)**2 / model
                    chisq.append(abs(diff))
                # update the master dataframe
                fwhm = xF[2] * 2.355 * 2614.5 / mu
                fwhmerr = xF_err[2] * 2.355 * 2614.5 / mu 
                rchi2 = sum(np.array(chisq) / len(hE))
                
                df_grid.at[i, f"fwhm_{ged}"] = fwhm
                df_grid.at[i, f"fwhmerr_{ged}"] = fwhmerr
                df_grid.at[i, f"rchi2_{ged}"] = rchi2
            except:
                print("Fit not computed for detector",ged,"and entry",i)
                
            if verbose:
                plt.cla()
                plt.plot(xE, pf.gauss_step(xE, *xF), c='r', label='peakshape')
                gaus, step = pf.gauss_step(xE, *xF, components=True)
                gaus = np.array(gaus)
                step = np.array(step)
                plt.plot(xE, gaus, ls="--", lw=2, c='g', label="gaus")
                plt.plot(xE, step, ls='--', lw=2, c='m', label='step + bg')
                plt.plot(xE[1:], hE, lw=1, c='b', label=f"data {ged}")
                plt.xlabel(f"ADC channels", ha='right', x=1)
                plt.ylabel("Counts", ha='right', y=1)
                plt.legend(loc=2, fontsize=10,title=f"FWHM = {fwhm:.2f} $\pm$ {fwhmerr:.2f} keV")
                plt.show()
                
                
            # write the updated df_grid to the output file.  
            if not verbose:
                df_grid.to_hdf(f_grid, key="pygama_optimization")

    if not verbose:
        print("Update grid file:",f_grid,"with detector",ged)
        print(df_grid)
            
def plot_fwhm(f_grid,f_opt,d_out,efilter, verbose=False):
    """
    select the best energy resolution, plot best result fit and fwhm vs parameters
    """
    print("Grid file:",f_grid)
    df_grid = pd.read_hdf(f_grid)
    f_res = f"{d_out}/{efilter}_results.h5"
    if 'trapE' in efilter: df = pd.DataFrame(columns=['ged','rise','flat','rc','fwhm','fwhmerr'])
    if efilter == 'zacE' or efilter == 'cuspE': df = pd.DataFrame(columns=['ged','sigma','flat','decay','fwhm','fwhmerr'])
    f = h5py.File(f_opt,'r')
    for chn, ged in enumerate(f.keys()):
        d_det = f"{d_out}/{ged}"
        try: os.mkdir(d_det)
        except: pass
        d_det = f"{d_det}/{efilter}"
        try: os.mkdir(d_det)
        except: pass
        
        data =  f[ged]['data']
        # find fwhm minimum values
        df_grid = df_grid.loc[(df_grid[f"rchi2_{ged}"]<20)&(df_grid[f"fwhm_{ged}"]>0)]
        minidx = df_grid[f'fwhm_{ged}'].idxmin()
        df_min = df_grid.loc[minidx]
        #try:
        #plot best result fit
        energies = data[f"{efilter}_{minidx}"][()]
        mean = np.mean(energies)
        bins = 12000
        hE, xE, vE = ph.get_hist(energies,bins,(mean/2,mean*2))
        mu = xE[np.argmax(hE)]
        hmax = hE[np.argmax(hE)]
        idx = np.where(hE > hmax/2)
        ilo, ihi = idx[0][0], idx[0][-1]
        sig = (xE[ihi] - xE[ilo]) / 2.355
        idx = np.where(((xE-mu) > -8 * sig) & ((xE-mu) < 8 * sig))
        ilo, ihi = idx[0][0], idx[0][-1]
        xE, hE, vE = xE[ilo:ihi+1], hE[ilo:ihi], vE[ilo:ihi]
        x0 = [hmax, mu, sig, 1, 0]
        xF, xF_cov = pf.fit_hist(pf.gauss_step, hE, xE, var=vE, guess=x0)
        xF_err = np.sqrt(np.diag(xF_cov))
        fwhm = xF[2] * 2.355 * 2614.5 / mu
        fwhmerr = xF_err[2] * 2.355 * 2614.5 / mu 
        plt.plot(xE, pf.gauss_step(xE, *xF), c='r', label='peakshape')
        gaus, step = pf.gauss_step(xE, *xF, components=True)
        gaus = np.array(gaus)
        step = np.array(step)
        plt.plot(xE, gaus, ls="--", lw=2, c='g', label="gaus")
        plt.plot(xE, step, ls='--', lw=2, c='m', label='step + bg')
        plt.plot(xE[1:], hE, lw=1, c='b', label=f"data {ged}")
        plt.xlabel(f"ADC channels", ha='right', x=1)
        plt.ylabel("Counts", ha='right', y=1)
        plt.legend(loc=2, fontsize=10,title=f"FWHM = {fwhm:.2f} $\pm$ {fwhmerr:.2f} keV")
        plt.savefig(f"{d_det}/Fit_{ged}-{efilter}.pdf")
        plt.cla()
        #except: continue
        if efilter=='zacE' or efilter=='cuspE':
            #try:
            sigma, flat, decay = df_min[:3]
            results = [ged, f'{sigma:.2f}', f'{flat:.2f}', f'{decay:.2f}', f'{fwhm:.2f}', f'{fwhmerr:.2f}']
            # 1. vary the sigma cusp
            df_sigma = df_grid.loc[(df_grid.flat==flat)&(df_grid.decay==decay)&(df_grid.decay==decay)]
            x, y, err =  df_sigma['sigma'], df_sigma[f'fwhm_{ged}'], df_sigma[f'fwhmerr_{ged}']
            plt.errorbar(x,y,err,fmt='o')
            plt.xlabel("Sigma Cusp ($\mu$s)", ha='right', x=1)
            plt.ylabel(r"FWHM (keV)", ha='right', y=1)
            plt.savefig(f"{d_det}/FWHM_vs_Sigma_{ged}-{efilter}.pdf")
            plt.cla()
            # 2. vary the flat time
            df_flat = df_grid.loc[(df_grid.sigma==sigma)&(df_grid.decay==decay)]
            x, y, err =  df_flat['flat'], df_flat[f'fwhm_{ged}'], df_flat[f'fwhmerr_{ged}']
            plt.errorbar(x,y,err,fmt='o')
            plt.xlabel("Flat Top ($\mu$s)", ha='right', x=1)
            plt.ylabel("FWHM (keV)", ha='right', y=1)
            plt.savefig(f"{d_det}/FWHM_vs_Flat_{ged}-{efilter}.pdf")
            plt.cla() 
            # 3. vary the rc constant
            df_decay = df_grid.loc[(df_grid.sigma==sigma)&(df_grid.flat==flat)]
            x, y, err =  df_decay[f'decay'], df_decay[f'fwhm_{ged}'], df_decay[f'fwhmerr_{ged}']
            plt.errorbar(x,y,err,fmt='o')
            plt.xlabel("Decay constant ($\mu$s)", ha='right', x=1)
            plt.ylabel(r"FWHM (keV)", ha='right', y=1)
            plt.savefig(f"{d_det}/FWHM_vs_Decay_{ged}-{efilter}.pdf")
            plt.cla()
            #except:
            #print("")
        if 'trapE' in efilter:
            rise, flat, rc = df_min[:3]
            results = [ged, f'{rise:.2f}', f'{flat:.2f}', f'{rc:.2f}', f'{fwhm:.2f}', f'{fwhmerr:.2f}']
            # 1. vary the rise time
            df_rise = df_grid.loc[(df_grid.flat==flat)&(df_grid.rc==rc)]
            x, y, err =  df_rise['rise'], df_rise[f'fwhm_{ged}'], df_rise[f'fwhmerr_{ged}']
            #plt.plot(x,y,".b")
            plt.errorbar(x,y,err,fmt='o')
            plt.xlabel("Ramp time ($\mu$s)", ha='right', x=1)
            plt.ylabel(r"FWHM (kev)", ha='right', y=1)
            # plt.ylabel(r"FWHM", ha='right', y=1)
            plt.savefig(f"{d_det}/FWHM_vs_Rise_{ged}-{efilter}.pdf")
            plt.cla()
            
            # 2. vary the flat time
            df_flat = df_grid.loc[(df_grid.rise==rise)&(df_grid.rc==rc)]
            x, y, err =  df_flat['flat'], df_flat[f'fwhm_{ged}'], df_flat[f'fwhmerr_{ged}']
            #plt.plot(x,y,'.b')
            plt.errorbar(x,y,err,fmt='o')
            plt.xlabel("Flat time ($\mu$s)", ha='right', x=1)
            plt.ylabel("FWHM (keV)", ha='right', y=1)
            plt.savefig(f"{d_det}/FWHM_vs_Flat_{ged}-{efilter}.pdf")
            plt.cla() 
            # 3. vary the rc constant
            df_rc = df_grid.loc[(df_grid.rise==rise)&(df_grid.flat==flat)]
            x, y, err =  df_rc['rc'], df_rc[f'fwhm_{ged}'], df_rc[f'fwhmerr_{ged}']
            #plt.plot(x,y,'.b')
            plt.errorbar(x,y,err,fmt='o')
            plt.xlabel("RC constant ($\mu$s)", ha='right', x=1)
            plt.ylabel(r"FWHM (keV)", ha='right', y=1)
            plt.savefig(f"{d_det}/FWHM_vs_RC_{ged}-{efilter}.pdf")
            plt.cla()
        df.loc[chn] = results
    print("Results file:",f_res)
    df.to_hdf(f_res, key='results',mode='w')
    print(df)
    dets = range(len(df['fwhm']))
    fwhm = np.array([float(df['fwhm'][i]) for i in dets])
    fwhm_err = np.array([float(df['fwhmerr'][i]) for i in dets])
    plt.cla()
    plt.errorbar(dets,fwhm,fwhm_err,fmt='o',c='red',label=f'{efilter} filter')
    plt.xlabel("detector number", ha='right', x=1)
    plt.ylabel("FWHM (keV)", ha='right', y=1)
    plt.legend()
    plt.savefig(f"{d_out}/FWHM_{efilter}.pdf")

def compare_fwhm( d_out, efilter1, efilter2, efilter3, verbose=False):
    print(f'Comparing FWHM using {efilter1} and {efilter2} filters')
    df_1 = pd.read_hdf(f'{d_out}/{efilter1}_results.h5',key='results')
    df_2 = pd.read_hdf(f'{d_out}/{efilter2}_results.h5',key='results')
    dets = range(len(df_1['fwhm']))
    fwhm_1 = np.array([float(df_1['fwhm'][i]) for i in dets])
    fwhm_1_err = np.array([float(df_1['fwhmerr'][i]) for i in dets])
    fwhm_2 = np.array([float(df_2['fwhm'][i]) for i in dets])
    fwhm_2_err = np.array([float(df_2['fwhmerr'][i]) for i in dets])
    fwhm_diff = 100*(fwhm_1 - fwhm_2)/fwhm_1
    fwhm_diff_err =  100*np.sqrt(np.square(fwhm_1_err)  + np.square(fwhm_2_err) )/fwhm_1
    plt.figure(1)
    plt.errorbar(dets,fwhm_1,fwhm_1_err,fmt='o',c='red',label=f'{efilter1} filter')
    plt.errorbar(dets,fwhm_2,fwhm_2_err,fmt='o',c='blue',label=f'{efilter2} filter')
    if efilter3 is not 0:
        df_3 = pd.read_hdf(f'{d_out}/{efilter3}_results.h5',key='results')
        fwhm_3 = np.array([float(df_3['fwhm'][i]) for i in dets])
        fwhm_3_err = np.array([float(df_3['fwhmerr'][i]) for i in dets])
        fwhm_2_diff = 100*(fwhm_1 - fwhm_3)/fwhm_1
        fwhm_2_diff_err =  100*np.sqrt(np.square(fwhm_1_err)  + np.square(fwhm_3_err) )/fwhm_1
        plt.errorbar(dets,fwhm_3,fwhm_3_err,fmt='o',c='green',label=f'{efilter3} filter')
    plt.xlabel("detector number", ha='right', x=1)
    plt.ylabel("FWHM (keV)", ha='right', y=1)
    plt.legend()
    plt.savefig(f"{d_out}/FWHM_compare_{efilter1}-{efilter2}-{efilter3}.pdf")
    if verbose: plt.show(block=False)
    plt.figure(2)
    plt.cla()
    plt.errorbar(dets,fwhm_diff,fwhm_diff_err,fmt='o',c='green',label=f'FWHM {efilter1} - {efilter2}')
    plt.axhline(0, c='black', linestyle=":")
    plt.xlabel("detector number", ha='right', x=1)
    plt.ylabel("FWHM difference (%)", ha='right', y=1)
    plt.legend()
    plt.savefig(f"{d_out}/FWHM_diff_{efilter1}-{efilter2}.pdf")
    if verbose: plt.show()
    
    
if __name__=="__main__":
    main()
