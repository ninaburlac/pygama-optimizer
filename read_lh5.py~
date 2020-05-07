import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pygama.analysis.histograms as pgh
import pygama.analysis.peak_fitting as pgpf
import pygama.io.lh5 as lh5
import h5py, sys
import argparse

mode     = 1
filename = './tier1/t1_run118.lh5'
nevts    = 1
fevt     = 0
nsamp    = 5592
chn      = 60

par = argparse.ArgumentParser(description="lh5 FlashCam data parser")
arg, st, sf = par.add_argument, "store_true", "store_false"
arg("-f", "--file" , nargs=1, help="tier filename")
arg("-m", "--mode" , nargs=1, help="[1=open tier1-2=open tier2")
arg("-n", "--nevts", nargs=1, help="number of events to show")
arg("-e", "--fevt" , nargs=1, help="first event to show")
arg("-s", "--samples" , nargs=1, help="number of waveform samples")
arg("-c", "--channel" , nargs=1, help="detector number")
arg("-z", "--zac" , action=st, help="use ZAC energy")
args = vars(par.parse_args())
if args['file']: filename = args['file'][0]
if args['mode']: mode     = int(args['mode'][0])
if args['nevts']: nevts   = int(args['nevts'][0])
if args['fevt']: fevt     = int(args['fevt'][0])
if args['samples']: nsamp = int(args['samples'][0])
if args['channel']: chn = int(args['channel'][0])
zac = 0
if args['zac']: zac = 1
if mode == 0:
    f = h5py.File(filename,'r')
    #for ged in f.keys():
        
    dset = f['g0%d/raw'%chn]
    print("Plot Waveforms for channel:",chn)
    #old_evt = fevt

    #print(dset['tracelist'].keys())
    #for evt,t,chn,bl,ene,nt in zip(dset['ievt'],dset['timestamp'],dset['channel'],dset['baseline'],dset['energy'],dset['numtraces']):
    #print("  nwf %d evt %d time %f chn %d bl %d energy %d" % (nwf,evt,t,chn,bl,ene))
    
    #    if nwf < fevt:
    #        nwf = nwf+1
    #        continue
    #    if nwf > fevt + nevts: break
    for nwf in range(0,1000):
        #waveform = dset['waveform']['values'][nwf*nsamp:(nwf+1)*nsamp]
        #plt.plot(np.array(waveform,dtype=np.int32)-bl, label="ch %d - bl = %d" %(chn,bl))
        #print(dset['numtraces'][nwf])
        #print(dset['ievt'][nwf])
        try:
            waveform = dset['waveform']['values'][nwf*nsamp:(nwf+1)*nsamp][0]#[dset['ievt'][()]<10]
            plt.plot(np.array(waveform), label="chn %d" % chn)
        except: continue

    plt.xlabel("samples")
    plt.ylabel("ADC")
    plt.legend()#title=("run117 - chn %d" % chn))
    #plt.savefig("./waveform_chn%d.pdf" % chn )
    plt.show()
    #nwf = nwf+1

     

if mode == 1:
    print("Open tier1: ", filename)
    f = h5py.File(filename,'r')
    print("File info: ",f.keys())
    chn = 0
    for ged in f.keys():
        try: 
            dset = f[ged]['raw']
            print("key: ",ged,"Data info: ",dset.keys())
        except:     
            conf = f[ged]
            print("Header info: ",conf.keys())
            
        energies = dset['energy'][()]
        #energies = f[ged]['raw/energy'][()]
        #energies = dset['energy'][dset['channel'][()]==chn]
        maxe = np.amax(energies)
        h, b, v = pgh.get_hist(energies, bins=3500, range=(maxe/4,maxe))
        pgh.plot_hist(h, b, label=ged )
        #bin_max = b[np.where(h == h.max())][0]
        #print("chn %d Raw energy max %d, histogram max %d at %d " % (chn,maxe,h.max(),bin_max ))
        #min_ene = int(bin_max*0.95)
        #max_ene = int(bin_max*1.05)
        #hist, bins, var = pgh.get_hist(energies, bins=500, range=(min_ene, max_ene))
        #pgh.plot_hist(hist, bins, label="chn %d" % chn )
        #chn = chn + 1
        
    plt.xlabel("uncalibrated energy", ha='right', x=1)
    plt.ylabel("counts", ha='right', y=1)
    plt.yscale('log')
    plt.legend()
    plt.show()
    #plt.savefig("./peak_chn%d.pdf" % chn )
    #plt.cla()
        

if mode == 2:
    print("Open tier2: ", filename)
    #st2 = lh5.Store()
    #st2.ls(filename)
    #df2 = st2.read_object('data', filename).get_dataframe()
    f = h5py.File(filename,'r')
    print("File info: ",f.keys())
    chn = 0
    for ged in f.keys():
        data = f[ged]['data']
        if chn==0: print(ged," Tier2 info: ", data.keys())
        try:
            if zac==1:
                energy = data["zacE"][()]
                print("Using ZAC energy")
            else:
                energy = data["trapEmax"][()]
            #print(energy)
            histene, binse, vare = pgh.get_hist(energy,3500,(0,50000))
            pgh.plot_hist(histene,binse,label=ged,show_stats=True)    
            #plt.step(np.concatenate(([binse[0]], binse)), np.concatenate(([0], histene, [0])), where="post",label="trap energy")
            #plt.plot(binse, pgpf.gauss(binse,coeffe[0],coeffe[1],coeffe[2]), "-r",label="FWHM = {0:.2f} keV".format(fwhm_std))
        except:
            print("Energy not find in",ged)
        # fit 2.6 MeV peak
        #coeffe, cov_matrixe = pgpf.fit_hist(pgpf.gauss,histene,binse,vare,(40880,10,500))
        #pgpf.gauss(coeffe[0],coeffe[1],coeffe[2])
        #fwhm_std = coeffe[1]/coeffe[0]*2614.5*2.35
        chn = chn + 1
    
    #plt.figure(1)
    #plt.cla()
    plt.xlabel("uncalibrated energy", ha='right', x=1)
    plt.ylabel("counts", ha='right', y=1)
    plt.legend()
    plt.yscale('log')
    #plt.tight_layout()
    #plt.show(block=False)
    plt.show()
    #plt.grid(True)                      
    #plt.pause(0.0001)

"""
hist, bins, var = pgh.get_hist(zacenergy,300,(40250,41250)) #zac energy spectrum

# fit 2.6 MeV peak
coeff, cov_matrix = pgpf.fit_hist(pgpf.gauss,hist,bins,var,(40750,10,500))
pgpf.gauss(coeff[0],coeff[1],coeff[2])
fwhm_zac = coeff[1]/coeff[0]*2614.5*2.35

# plot

plt.figure(2)
plt.cla()
plt.step(np.concatenate(([bins[0]], bins)), np.concatenate(([0], hist, [0])), where="post",label="zac energy")
plt.plot(bins, pgpf.gauss(bins,coeff[0],coeff[1],coeff[2]), "-r",label="FWHM = {0:.2f} keV".format(fwhm_zac))
plt.xlabel("uncalibrated energy", ha='right', x=1)
plt.ylabel("counts", ha='right', y=1)
#plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.show()
plt.pause(0.0001)
"""
