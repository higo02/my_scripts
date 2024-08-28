#Author: Higo de Araujo Oliveira    Nov 6,2023
#------------------------------------------Imports----------------------------------
import numpy as np
import matplotlib.pyplot as plt
from thermo.gpumd.data import load_kappa,load_shc
from thermo.gpumd.calc import running_ave, hnemd_spectral_kappa
import sys
#------------------------------------------Functions--------------------------------
def export_to_file(filepath,td,t,kappa_dict,ra='yes'):
    if ra=='yes': #ra==yes means this function is being used to export single run running_ave
        kappa_ra={}
        kappa_ra['kxi_ra']=running_ave(kappa_dict['kxi'],t)
        kappa_ra['kxo_ra']=running_ave(kappa_dict['kxo'],t)
        kappa_ra['kyi_ra']=running_ave(kappa_dict['kyi'],t)
        kappa_ra['kyo_ra']=running_ave(kappa_dict['kyo'],t)
        kappa_ra['kz_ra']=running_ave(kappa_dict['kz'],t)
    else: #if ra!=yes, the input dict itself is exported
        kappa_ra=kappa_dict

    with open(filepath,'w') as text:
        text.write("#kappa is in units of W/m/K \n")
        text.write("#transport direction is: "+str(td)+"\n")
        text.write("#t(ns)   kappa_"+str(td)+"x_in  kappa_"+str(td)+"x_out  kappa_"+str(td)+"x_total    kappa_"+str(td)+"y_in    kappa_"+str(td)+"y_out kappa_"+str(td)+"y_total    kappa_"+str(td)+"z_total \n")
        for idx in range(len(t)):
            kappa_xtotal=kappa_ra['kxi_ra'][idx]+kappa_ra['kxo_ra'][idx]
            kappa_ytotal=kappa_ra['kyi_ra'][idx]+kappa_ra['kyo_ra'][idx]
            text.write(f"{t[idx]:.4f}"+"    "+f"{kappa_ra['kxi_ra'][idx]:.4f}"+"   "+f"{kappa_ra['kxo_ra'][idx]:.4f}"+"   "+f"{kappa_xtotal:.4f}"+"   "+f"{kappa_ra['kyi_ra'][idx]:.4f}"+"   "+f"{kappa_ra['kyo_ra'][idx]:.4f}"+"   "+f"{kappa_ytotal:.4f}"+"   "+f"{kappa_ra['kz_ra'][idx]:.4f}"+"\n")

#--------------------------------------------------------------Main-----------------------------------------------------------
if len(sys.argv)!=3:
    print('main.py transp_dir=x or y num_of_directories')
    exit()
td=sys.argv[1]
ndir=int(sys.argv[2])
cf=1# volume correction factor for z_dir (due to vacuum, if present), cf=1 means no correction

#initialize averages dict with files in folder run0:
kappa_ave={}
ktdi_ra={} #dict for kappa in plane for direction td (x or y)
ktdo_ra={} #dict for kappa out of plane for direction td (x or y)

kappa=load_kappa(directory='run0')
t=np.arange(1,len(kappa['kxi'])+1)*1E-3 #convert to ns
export_to_file('run0/kappa_ra_run0.dat',td,t,kappa)

kappa_ave['kyi_ra'] = running_ave(kappa['kyi'],t)
kappa_ave['kyo_ra'] = running_ave(kappa['kyo'],t)
kappa_ave['kxi_ra'] = running_ave(kappa['kxi'],t)
kappa_ave['kxo_ra'] = running_ave(kappa['kxo'],t)
kappa_ave['kz_ra'] = running_ave(kappa['kz'],t)

ktdi_ra['0'] = running_ave(kappa['k'+td+'i'],t)
ktdo_ra['0'] = running_ave(kappa['k'+td+'o'],t)

#runs from run(1) to the run(ndir-1) folder:
for idx in range(1,ndir):
    print(idx)
    kappa=load_kappa(directory='run'+str(idx))
    print(kappa.keys()) 
    export_to_file('run'+str(idx)+'/kappa_ra_run'+str(idx)+'.dat',td,t,kappa)

    kappa_ave['kyi_ra'] += running_ave(kappa['kyi'],t)
    kappa_ave['kyo_ra'] += running_ave(kappa['kyo'],t)
    kappa_ave['kxi_ra'] += running_ave(kappa['kxi'],t)
    kappa_ave['kxo_ra'] += running_ave(kappa['kxo'],t)
    kappa_ave['kz_ra'] += running_ave(kappa['kz'],t)

    ktdi_ra[str(idx)] = running_ave(kappa['k'+td+'i'],t)
    ktdo_ra[str(idx)] = running_ave(kappa['k'+td+'o'],t)

#applies the correction factor and the normalization constant (for the average)
for key in kappa_ave:
    kappa_ave[key]=kappa_ave[key]*cf
    kappa_ave[key]/=ndir


#-------------------------------------------print final time statiscs-------------------------------
k_in=np.array([ktdi_ra[str(i)][-1] for i in range(ndir)])
k_out=np.array([ktdo_ra[str(i)][-1] for i in range(ndir)])
k_total=k_in+k_out


k_total_mean=np.mean(k_total)
k_in_mean=np.mean(k_in)
k_out_mean=np.mean(k_out)

sigma_total=np.std(k_total)
sigma_in=np.std(k_in)
sigma_out=np.std(k_out)

print('k_total_mean='+str(k_total_mean))
print('sigma_total='+str(sigma_total))

print('\n k_in_mean='+str(k_in_mean))
print('sigma_in='+str(sigma_in))

print('\n k_in_out='+str(k_out_mean))
print('sigma_out='+str(sigma_out))

#-------------------------------export ave over different runs to a file----------------------------
export_to_file('kappa_ave.dat',td,t,kappa_ave,ra='no')
