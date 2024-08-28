import numpy as np
import matplotlib.pyplot as plt
from thermo.gpumd.data import load_kappa,load_shc
from thermo.gpumd.calc import running_ave, hnemd_spectral_kappa
import sys

if len(sys.argv)!=3:
    print('main.py transp_dir=x or y num_of_directories')
    exit()
td=sys.argv[1]
ndir=int(sys.argv[2])
cf=(5.431/5.431)# volume correction factor for z_dir (due to vacuum, if present)

#initialize averages dict with files in folder run0:
kappa_ave={}
ktdi_ra={} #dict for kappa in plane for direction td (x or y)
ktdo_ra={} #dict for kappa out of plane for direction td (x or y)

kappa=load_kappa(directory='run0')
t=np.arange(1,len(kappa['kxi'])+1)*1E-3 #convert to ns
kappa_ave['kyi_ra'] = running_ave(kappa['kyi'],t)
kappa_ave['kyo_ra'] = running_ave(kappa['kyo'],t)
kappa_ave['kxi_ra'] = running_ave(kappa['kxi'],t)
kappa_ave['kxo_ra'] = running_ave(kappa['kxo'],t)
kappa_ave['kz_ra'] = running_ave(kappa['kz'],t)

ktdi_ra['0'] = running_ave(kappa['k'+td+'i'],t)
ktdo_ra['0'] = running_ave(kappa['k'+td+'o'],t)
#runs from run1 to the runN folder:
for idx in range(1,ndir):
    print(idx)
    kappa=load_kappa(directory='run'+str(idx))
    print(kappa.keys()) 
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


#------------------------------statiscs-------------------------------
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

#--------------------------------save total heat conductivity in a dat file-----------------------
k_in_t=np.zeros(len(t))
k_out_t=np.zeros(len(t))
for i in range(ndir):
    k_in_t+=np.array(ktdi_ra[str(i)])
    k_out_t+=np.array(ktdo_ra[str(i)])
k_total_t=(k_in_t+k_out_t)/ndir
print(k_total_t)

with open('kappa_hnemd_70x70.dat','w') as arquivo:
    arquivo.write('#t    kappa (W/m/K) \n')
    for idx in range(len(t)):
        arquivo.write(str(t[idx])+'     '+str(k_total_t[idx])+'\n')
#exit()
#---------------------------------plot figure 1--------------------------------------
#plt.rc('text',usetex=True)
#plt.rc('font',family='serif')

#units for plot
wmk=' (Wm$^{-1}$K$^{-1}$) '

fig=plt.figure(1,figsize=(10,8))
ax1=fig.add_subplot(221)
#ax1.plot(t,kappa_ave['k'+td+'i'],color='C8',alpha=0.3)
ax1.plot(t,kappa_ave['k'+td+'i_ra'],color='red',linewidth=2)
ax1.grid(True,linestyle='--',color='gray',alpha=0.5)
#ax1.set_xlim([0,20])
ax1.set_ylim([-20,25])
ax1.set_xlabel('time (ns)')
ax1.set_ylabel(r'$\kappa_{'+td+',in}$'+wmk)
ax1.set_title('a')

ax1=fig.add_subplot(222)
#ax1.plot(t,kappa_ave['k'+td+'o'],color='C8',alpha=0.3)
ax1.plot(t,kappa_ave['k'+td+'o_ra'],color='blue',linewidth=2)
ax1.grid(True,linestyle='--',color='gray',alpha=0.5)
#ax1.set_xlim([0,20])
ax1.set_ylim([-20,25])
ax1.set_xlabel('time (ns)')
ax1.set_ylabel(r'$\kappa_{'+td+',out}$'+wmk)
ax1.set_title('b')


ax1=fig.add_subplot(223)
ax1.plot(t,kappa_ave['k'+td+'i_ra'],color='red',linewidth=2)
ax1.plot(t,kappa_ave['k'+td+'o_ra'],color='blue',linewidth=2)
ax1.plot(t,kappa_ave['k'+td+'i_ra']+kappa_ave['k'+td+'o_ra'],color='black',linewidth=2)
ax1.grid(True,linestyle='--',color='gray',alpha=0.5)
#ax1.set_xlim([0,20])
ax1.set_ylim([-20,25])
#ax1.set_yticks(range(-20,80,10))
ax1.set_xlabel('time (ns)')
ax1.set_ylabel(r'$\kappa$'+wmk)
ax1.legend([td+' in',td+' out',td+' total'],loc='best')
ax1.set_title('c')

ax1=fig.add_subplot(224)
ax1.plot(t,kappa_ave['kxi_ra']+kappa_ave['kxo_ra'],color='red',linewidth=2)
ax1.plot(t,kappa_ave['kyi_ra']+kappa_ave['kyo_ra'],color='blue',linewidth=2)
ax1.plot(t,kappa_ave['kz_ra'],color='black',linewidth=2)
ax1.grid(True,linestyle='--',color='gray',alpha=0.5)
#ax1.set_xlim([0,20])
ax1.set_ylim([-20,25])
#ax1.set_yticks(range(-70,70,10))
ax1.set_xlabel('time (ns)')
ax1.set_ylabel(r'$\kappa$'+wmk)
ax1.legend(['x'+td,'y'+td,'z'+td],loc='best')
ax1.set_title('d')

#---------------second figure------------------------------
fig=plt.figure(2,figsize=(10,8))
ax1=fig.add_subplot(221)
legend=[]
for key in ktdi_ra:
    ax1.plot(t,ktdi_ra[key])
    legend.append(key)
ax1.set_xlabel('time (ns)')
ax1.set_ylabel(r'$\kappa_{'+td+',in}$'+wmk)
ax1.set_ylim([-80,80])
#ax1.legend(legend,loc='lower left')
ax1.set_title('a')

ax1=fig.add_subplot(222)
legend=[]
for key in ktdo_ra:
    ax1.plot(t,ktdo_ra[key])
    legend.append(key)
ax1.set_xlabel('time (ns)')
ax1.set_ylabel(r'$\kappa_{'+td+',out}$'+wmk)
ax1.set_ylim([-80,80])
#ax1.legend(legend,loc='lower left')
ax1.set_title('b')

ax1=fig.add_subplot(223)
legend=[]
for key in ktdi_ra:
    ax1.plot(t,ktdi_ra[key]+ktdo_ra[key])
    legend.append(key)
ax1.set_xlabel('time (ns)')
ax1.set_ylabel(r'$\kappa_{x,total}$'+wmk)
ax1.set_ylim([-80,80])
#ax1.legend(legend,loc='lower left')
ax1.set_title('c')


plt.tight_layout()
plt.show()

