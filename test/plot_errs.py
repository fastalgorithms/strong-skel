from numpy import *
from pylab import *

nu = ['7','10','20','40','60','80']
norders = ['3','5','7']


tf = zeros((6,3))
ts = zeros((6,3))
np = zeros((6,3))
tq = zeros((6,3))
mf = zeros((6,3))
errs = zeros((6,3))

for inu in range(6):
    for ino in range(3):
    
        fname = 'data_nov15_2021/diary_ik1_np'+str(nu[inu])+'_norder'+str(norders[ino])+'.dat'
        f = open(fname,'r')
        flines = f.readlines()
        errs[inu,ino] = float(flines[-1].split(':')[1])
        ts[inu,ino] = float(flines[-2].split(':')[1])
        tf[inu,ino] = float(flines[-3].split(':')[1])
        tq[inu,ino] = float(flines[-4].split(':')[1])
        mf[inu,ino] = float(flines[-10].split(':')[1][0:7])
        np[inu,ino] = float(nu[inu])**2*2
        f.close()

figure(1)
clf()
loglog(np[:,0],errs[:,0],'ko',label='p=4')
loglog(np[:,1],errs[:,1],'ks',label='p=6')
loglog(np[:,2],errs[:,2],'k^',label='p=8')
loglog(np[:,2],5e-7*ones(size(np[:,2])),'k--')
title('Fixed k')
ylim((1e-9,1e-2))
legend()


tf2 = zeros(6)
ts2 = zeros(6)
np2 = zeros(6)
tq2 = zeros(6)
mf2 = zeros(6)
errs2 = zeros(6)

for inu in range(1,6):
    
    fname = 'data_nov15_2021/diary_ik2_np'+str(nu[inu])+'_norder5.dat'
    f = open(fname,'r')
    flines = f.readlines()
    errs2[inu] = float(flines[-1].split(':')[1])
    ts2[inu] = float(flines[-2].split(':')[1])
    tf2[inu] = float(flines[-3].split(':')[1])
    tq2[inu] = float(flines[-4].split(':')[1])
    mf2[inu] = float(flines[-10].split(':')[1][0:7])
    np2[inu] = float(nu[inu])**2*2
    f.close()

figure(2)
clf()
loglog(np2,errs2,'ks',label='p=6')
loglog(np[:,2],5e-7*ones(size(np[:,2])),'k--')
title('Fixed ppw')
ylim((1e-9,1e-2))



tf = zeros((6,3))
ts = zeros((6,3))
np = zeros((6,3))
mf = zeros((6,3))
errs = zeros((6,3))

for inu in range(6):
    for ino in range(3):
    
        fname = 'diary_ik1_np'+str(nu[inu])+'_norder'+str(norders[ino])+'_ss.dat'
        f = open(fname,'r')
        flines = f.readlines()
        errs[inu,ino] = float(flines[-1].split(':')[1])
        ts[inu,ino] = float(flines[-2].split(':')[1])
        tf[inu,ino] = float(flines[-3].split(':')[1])
        mf[inu,ino] = float(flines[-10].split(':')[1][0:7])
        np[inu,ino] = float(nu[inu])**2*2
        f.close()

figure(3)
clf()
loglog(np[:,0],errs[:,0],'ko',label='p=4')
loglog(np[:,1],errs[:,1],'ks',label='p=6')
loglog(np[:,2],errs[:,2],'k^',label='p=8')
loglog(np[:,2],5e-7*ones(size(np[:,2])),'k--')
title('Fixed k')
ylim((1e-9,1e-2))
legend()


tf2 = zeros(6)
ts2 = zeros(6)
np2 = zeros(6)
mf2 = zeros(6)
errs2 = zeros(6)

for inu in range(6):
    
    fname = 'diary_ik2_np'+str(nu[inu])+'_norder5_ss.dat'
    f = open(fname,'r')
    flines = f.readlines()
    errs2[inu] = float(flines[-1].split(':')[1])
    ts2[inu] = float(flines[-2].split(':')[1])
    tf2[inu] = float(flines[-3].split(':')[1])
    mf2[inu] = float(flines[-10].split(':')[1][0:7])
    np2[inu] = float(nu[inu])**2*2
    f.close()

figure(4)
clf()
loglog(np2,errs2,'ks',label='p=6')
loglog(np[:,2],5e-7*ones(size(np[:,2])),'k--')
title('Fixed ppw')
ylim((1e-9,1e-2))
show()


