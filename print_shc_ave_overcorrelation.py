import numpy as np
import matplotlib.pyplot as plt
from thermo.gpumd.data import load_kappa,load_shc
from thermo.gpumd.calc import running_ave, hnemd_spectral_kappa
import sys

if len(sys.argv)!=3:
    print('main.py transp_dir=x or y ndir')
    exit()
td=sys.argv[1]
ndir=int(sys.argv[2])

shc_ave={}
shc=load_shc(Nc=500,num_omega=1000,directory='run0')['run0']
for key in shc:
    shc_ave[key]=shc[key]
#print(shc_ave.keys())
#---------------read box lenghts from model.xyz only Lx and Ly , due to vacuum in Lz
arquivo=open('run0/model.xyz','r')
line=arquivo.readline()
while(line.find('pbc')==-1):
    line=arquivo.readline()
line_elements=line.split('\"')
line_elements=line_elements[3].split()
Lx=float(line_elements[0])
Ly=float(line_elements[4])
Lz=5.431
print(Lx,Ly,Lz)
arquivo.close()

for idx in range(1,ndir):
    shc=load_shc(Nc=500,num_omega=1000,directory='run'+str(idx))['run0']
    for key in shc:
        shc_ave[key]+=shc[key]
    #print(shc_ave.keys())
    #---------------read box lenghts from model.xyz only Lx and Ly , due to vacuum in Lz
    arquivo=open('run'+str(idx)+'/model.xyz','r')
    line=arquivo.readline()
    while(line.find('pbc')==-1):
        line=arquivo.readline()
    line_elements=line.split('\"')
    line_elements=line_elements[3].split()
    Lx+=float(line_elements[0])
    Ly+=float(line_elements[4])
    Lz+=5.431
    print(Lx,Ly,Lz)
    print(idx)
    arquivo.close()
#----averages----------------------------------
Lx/=ndir
Ly/=ndir
Lz/=ndir
for key in shc_ave:
    shc_ave[key]/=ndir

V_ave=Lx*Ly*Lz
T=300
Fe=1e-4
hnemd_spectral_kappa(shc_ave,Fe,T,V_ave)
shc_ave['kw']=shc_ave['kwi']+shc_ave['kwo']
shc_ave['K']=shc_ave['Ki']+shc_ave['Ko']
#Gc=np.load('Gc.npy')

#------------------------------HNEMD parameters for lenght dependence
#lambda_i=shc_ave['kw']/Gc
#L=np.logspace(1,8,100)
#k_Lhnemd=np.zeros_like(L)
#for idx,el in enumerate(L):
#    k_Lhnemd[idx]=np.trapz(shc_ave['kw']/(1+(lambda_i/el)),shc_ave['nu'])

#------------------------------langevin NEMD data-------------
#L_langev=[214.533]
#kappa_langev=[20.132]


#------------------------------MP NEMD data----------------------------
#L_MP=[238.55,317.479,395.490,474.816]
#kappa_MP=[13.662,14.785,16.893,18.592]
#-----------------------------write to file------------------------

with open('kappa_w.dat','w') as f:
    f.write('#nu(THz)       kappa(W/m/K/THz)\n')
    for idx in range(len(shc_ave['nu'])):
        f.write(str(shc_ave['nu'][idx])+'   '+str(shc_ave['kw'][idx])+'\n')

#with open('kappa_lambda_w.dat','w') as f:
#    f.write('#nu(THz)       kappa(W/m/K/THz)        lambda(nm/THz) \n')
#    for idx in range(len(shc_ave['nu'])):
#        f.write(str(shc_ave['nu'][idx])+'   '+str(shc_ave['kw'][idx])+'     '+str(lambda_i[idx])+'\n')

with open('kappa_in_out_w.dat','w') as f:
    f.write('#nu(THz)       kappa_in        kappa_out   (W/m/K) \n')
    for idx in range(len(shc_ave['nu'])):
        f.write(str(shc_ave['nu'][idx])+'   '+str(shc_ave['kwi'][idx])+'     '+str(shc_ave['kwo'][idx])+'\n')


#-------------------------plt settings------------
fig=plt.figure(1,figsize=(10,8))
ax1=fig.add_subplot(221)
ax1.plot(shc_ave['t'],shc_ave['K'],linewidth=2)
#ax1.set_xlim([-0.5,0.5])
ax1.set_ylabel('K (eV.Ang/ps)')
ax1.set_xlabel('Correlation time (ps)')
ax1.set_title('a')

ax1=fig.add_subplot(222)
ax1.plot(shc_ave['nu'],shc_ave['kw'],linewidth=2)
ax1.set_xlabel(r'$\nu$ (THz)')
ax1.set_ylabel(r'$\kappa$($\omega$) (W/m/K/THz)')
ax1.set_title('b')

#ax1=fig.add_subplot(223)
#ax1.plot(shc_ave['nu'],lambda_i,linewidth=2)
#ax1.set_xlabel(r'$\nu$ (THz)')
#x1.set_ylabel(r'$\lambda$($\omega$) (nm)')
#ax1.set_title('c')

#ax1=fig.add_subplot(224)
#ax1.semilogx(L,k_Lhnemd,color='blue',linewidth=3,label='HNEMD')
#ax1.semilogx(L_langev,kappa_langev,marker='o',linewidth=0,markersize=5,color='red',label='langevin')
#ax1.semilogx(L_MP,kappa_MP,marker='x',linewidth=0,markersize=5,color='orange',label='MP')
#ax1.set_xlabel(r' L (nm)')
#ax1.set_ylabel(r'$\kappa$ (W/m/K)')
#ax1.legend(loc='best')
#ax1.set_title('d')

plt.show()

