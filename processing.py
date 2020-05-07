#!/usr/bin/env python3
import sys, os, io
import json
import argparse
import numpy as np
import pandas as pd
from pprint import pprint
from pygama import DataSet
import pygama.utils as pu
from pygama.dsp.ProcessingChain import ProcessingChain
from pygama.dsp.processors import *
from pygama.dsp.units import *
from pygama.dsp.processors import *
from pygama.io import io_base as io
from pygama.io import lh5
import h5py, sys

def main(argv):
    """
    Uses pygama's amazing DataSet class to process runs
    for different data sets and arbitrary configuration options
    defined in a JSON file.
    """
    #datadir = os.environ["CAGEDATA"]
    run_db, cal_db = f'./meta/runDB.json', f'./meta/calDB.json'

    # -- parse args --
    par = argparse.ArgumentParser(description="data processing suite for CAGE")
    arg, st, sf = par.add_argument, "store_true", "store_false"
    arg("-ds", nargs='*', action="store", help="load runs for a DS")
    arg("-r", "--run", nargs=1, help="load a single run")
    arg("-d2r", "--daq_to_raw", action=st, help="run daq_to_raw on list")
    arg("-r2d", "--raw_to_dsp", action=st, help="run raw_to_dsp on list")
    arg("-t", "--test", action=st, help="test mode, don't run")
    arg("-n", "--nevt", nargs='?', default=np.inf, help="limit max num events")
    arg("-i", "--ioff", nargs='?', default=0, help="start at index [i]")
    arg("-o", "--ovr", action=st, help="overwrite existing files")

    arg('-v', '--verbose', default=2, type=int,
        help="Verbosity level: 0=silent, 1=basic warnings, 2=verbose output, 3=debug. Default is 2.")

    arg('-b', '--block', default=8, type=int,
        help="Number of waveforms to process simultaneously. Default is 8")

    arg('-g', '--group', default='',
        help="Name of group in LH5 file. By default process all base groups. Supports wildcards.")

    # -- declare the DataSet --
    args = par.parse_args()
    d_args = vars(par.parse_args())
    #ds = pu.get_dataset_from_cmdline(d_args, run_db, cal_db)
    # -- declare the DataSet --
    if d_args["ds"]:
        ds_lo = int(d_args["ds"][0])
        try:
            ds_hi = int(d_args["ds"][1])
        except:
            ds_hi = None
            ds = DataSet(1,ds_lo, ds_hi, md=run_db, v=d_args["verbose"])
            
    if d_args["run"]:
        ds = DataSet(1,run=int(d_args["run"][0]), md=run_db, v=d_args["verbose"])
    
    #print(ds.runs)
    #pprint(ds.paths)

    # -- start processing --
    if args.daq_to_raw:
        daq_to_raw(ds, args.ovr, args.nevt, args.verbose, args.test)

    if args.raw_to_dsp:
        raw_to_dsp(ds, args.ovr, args.nevt, args.test, args.verbose, args.block,
                   args.group)


def daq_to_raw(ds, overwrite=False, nevt=np.inf, v=False, test=False):
    """
    Run daq_to_raw on a set of runs.
    [raw file] ---> [raw_run{}.lh5] (basic info & waveforms)
    """
    from pygama.io.daq_to_raw import daq_to_raw

    for run in ds.runs:
        daq_file = "/lfs/l1/legend/users/dandrea/pygama/pgt/raw/pgt_longtrace_run0117-20200110-105115-calib.fcio"
        raw_file = "pgt_longtrace_run0117-20200110-105115-calib_raw.lh5"
        #daq_file = ds.paths[run]["daq_path"]
        #raw_file = ds.paths[run]["raw_path"]
        if raw_file is not None and overwrite is False:
            print("file exists, overwrite flag isn't set.  continuing ...")
            continue

        #conf = ds.paths[run]["build_opt"]
        #opts = ds.config["build_options"][conf]["daq_to_raw_options"]
        opts = ds.config["daq_to_raw"]
        print("opts",opts)
        if test:
            print("test mode (dry run), processing DAQ file:", daq_file)
            print("output file:", raw_file)
            continue

        # # old pandas version
        # daq_to_raw(daq_file, run, verbose=v, output_dir=ds.raw_dir,
        #              overwrite=overwrite, n_max=nevt, config=ds.config)#settings=opts)

        # new lh5 version
        daq_to_raw(daq_file, raw_filename=raw_file, run=run, chan_list=None,
                   prefix=ds.rawpre, n_max=nevt, verbose=v, output_dir=ds.raw_dir,
                   overwrite=overwrite, config=ds.config)


def raw_to_dsp(ds, overwrite=False, nevt=None, test=False, verbose=2, block=8, group=''):
    """
    Run raw_to_dsp on a set of runs.
    [raw file] ---> [dsp_run{}.lh5] (digital signal processing results)
    """
    for run in ds.runs:
        raw_file = "/lfs/l1/legend/users/dandrea/pygama/pgt/tier1/pgt_longtrace_run0117-20200110-105115-calib_raw.lh5"
        dsp_file = "/lfs/l1/legend/users/dandrea/pygama/pgt/tier2/pgt_longtrace_run0117-20200110-105115-calib_dsp.lh5"
        #raw_file = ds.paths[run]["raw_path"]
        #dsp_file = ds.paths[run]["dsp_path"]
        print("raw_file: ",raw_file)
        print("dsp_file: ",dsp_file)
        if dsp_file is not None and overwrite is False:
            continue

        if dsp_file is None:
            # declare new file name
            dsp_file = raw_file.replace('raw_', 'dsp_')
            
        if test:
            print("test mode (dry run), processing raw file:", raw_file)
            continue
            
        print("Definition of new LH5 version")
        #f_lh5 = lh5.Store()
        #data = f_lh5.read_object("raw", raw_file)
        #wf_in = data['waveform']['values'].nda
        #dt = data['waveform']['dt'].nda[0] * unit_parser.parse_unit(data['waveform']['dt'].attrs['units'])
        
        lh5_in = lh5.Store()
        #groups = lh5_in.ls(raw_file, group)
        f = h5py.File(raw_file,'r')
        print("File info: ",f.keys())
        for group in f.keys():
            print("Processing: " + raw_file + '/' + group)
            #data = lh5_in.read_object(group, raw_file)
            data =  f[group]['raw']
            
            #wf_in = data['waveform']['values'].nda
            #dt = data['waveform']['dt'].nda[0] * unit_parser.parse_unit(data['waveform']['dt'].attrs['units'])
            wf_in = data['waveform']['values'][()]
            dt = data['waveform']['dt'][0] * unit_parser.parse_unit(data['waveform']['dt'].attrs['units'])
            
            # Parameters for DCR calculation
            dcr_trap_int = 200
            dcr_trap_flat = 1000
            dcr_trap_startSample = 1200
            
            # Set up processing chain
            proc = ProcessingChain(block_width=block, clock_unit=dt, verbosity=verbose)
            proc.add_input_buffer("wf", wf_in, dtype='float32')
            
            # Basic Filters
            proc.add_processor(mean_stdev, "wf[0:1000]", "bl", "bl_sig")
            proc.add_processor(np.subtract, "wf", "bl", "wf_blsub")
            proc.add_processor(pole_zero, "wf_blsub", 145*us, "wf_pz")
            proc.add_processor(trap_norm, "wf_pz", 10*us, 5*us, "wf_trap")
            proc.add_processor(asymTrapFilter, "wf_pz", 0.05*us, 2*us, 4*us, "wf_atrap")
            
            # Timepoint calculation
            proc.add_processor(np.argmax, "wf_blsub", 1, "t_max", signature='(n),()->()', types=['fi->i'])
            proc.add_processor(time_point_frac, "wf_blsub", 0.95, "t_max", "tp_95")
            proc.add_processor(time_point_frac, "wf_blsub", 0.8, "t_max", "tp_80")
            proc.add_processor(time_point_frac, "wf_blsub", 0.5, "t_max", "tp_50")
            proc.add_processor(time_point_frac, "wf_blsub", 0.2, "t_max", "tp_20")
            proc.add_processor(time_point_frac, "wf_blsub", 0.05, "t_max", "tp_05")
            proc.add_processor(time_point_thresh, "wf_atrap[0:2000]", 0, "tp_0")
            
            # Energy calculation
            proc.add_processor(np.amax, "wf_trap", 1, "trapEmax", signature='(n),()->()', types=['fi->f'])
            proc.add_processor(fixed_time_pickoff, "wf_trap", "tp_0+(5*us+9*us)", "trapEftp")
            proc.add_processor(trap_pickoff, "wf_pz", 1.5*us, 0, "tp_0", "ct_corr")
            
            # Current calculation
            proc.add_processor(avg_current, "wf_pz", 10, "curr(len(wf_pz)-10, f)")
            proc.add_processor(np.amax, "curr", 1, "curr_amp", signature='(n),()->()', types=['fi->f'])
            proc.add_processor(np.divide, "curr_amp", "trapEftp", "aoe")

            # DCR calculation: use slope using 1000 samples apart and averaging 200
            # samples, with the start 1.5 us offset from t0
            proc.add_processor(trap_pickoff, "wf_pz", 200, 1000, "tp_0+1.5*us", "dcr_unnorm")
            proc.add_processor(np.divide, "dcr_unnorm", "trapEftp", "dcr")
            
            # Tail slope. Basically the same as DCR, except with no PZ correction
            proc.add_processor(linear_fit, "wf_blsub[3000:]", "wf_b", "wf_m")
            proc.add_processor(np.divide, "-wf_b", "wf_m", "tail_rc")            
            
            #add zac filter energy calculation
            sigma = 10*us
            flat = 1*us
            decay = 160*us
            proc.add_processor(zac_filter, "wf", sigma, flat, decay, "wf_zac(101, f)")
            proc.add_processor(np.amax, "wf_zac", 1, "zacE", signature='(n),()->()', types=['fi->f'])
            
            # Set up the LH5 output
            lh5_out = lh5.Table(size=proc._buffer_len)
            lh5_out.add_field("zacE", lh5.Array(proc.get_output_buffer("zacE"), attrs={"units":"ADC"}))
            lh5_out.add_field("trapEmax", lh5.Array(proc.get_output_buffer("trapEmax"), attrs={"units":"ADC"}))
            lh5_out.add_field("trapEftp", lh5.Array(proc.get_output_buffer("trapEftp"), attrs={"units":"ADC"}))
            lh5_out.add_field("ct_corr", lh5.Array(proc.get_output_buffer("ct_corr"), attrs={"units":"ADC*ns"}))
            lh5_out.add_field("bl", lh5.Array(proc.get_output_buffer("bl"), attrs={"units":"ADC"}))
            lh5_out.add_field("bl_sig", lh5.Array(proc.get_output_buffer("bl_sig"), attrs={"units":"ADC"}))
            lh5_out.add_field("A", lh5.Array(proc.get_output_buffer("curr_amp"), attrs={"units":"ADC"}))
            lh5_out.add_field("AoE", lh5.Array(proc.get_output_buffer("aoe"), attrs={"units":"ADC"}))
            lh5_out.add_field("dcr", lh5.Array(proc.get_output_buffer("dcr"), attrs={"units":"ADC"}))
            
            lh5_out.add_field("tp_max", lh5.Array(proc.get_output_buffer("tp_95", unit=us), attrs={"units":"us"}))
            lh5_out.add_field("tp_95", lh5.Array(proc.get_output_buffer("tp_95", unit=us), attrs={"units":"us"}))
            lh5_out.add_field("tp_80", lh5.Array(proc.get_output_buffer("tp_80", unit=us), attrs={"units":"us"}))
            lh5_out.add_field("tp_50", lh5.Array(proc.get_output_buffer("tp_50", unit=us), attrs={"units":"us"}))
            lh5_out.add_field("tp_20", lh5.Array(proc.get_output_buffer("tp_20", unit=us), attrs={"units":"us"}))
            lh5_out.add_field("tp_05", lh5.Array(proc.get_output_buffer("tp_05", unit=us), attrs={"units":"us"}))
            lh5_out.add_field("tp_0", lh5.Array(proc.get_output_buffer("tp_0", unit=us), attrs={"units":"us"}))
            lh5_out.add_field("tail_rc", lh5.Array(proc.get_output_buffer("tail_rc", unit=us), attrs={"units":"us"}))
            
            print("Processing:\n",proc)
            proc.execute()
            
            #groupname = group[:group.rfind('/')+1]+"data"
            groupname = group+"/data"
            print("Writing to: " + dsp_file + "/" + groupname)
            lh5_in.write_object(lh5_out, groupname, dsp_file)



if __name__ == "__main__":
    main(sys.argv[1:])
