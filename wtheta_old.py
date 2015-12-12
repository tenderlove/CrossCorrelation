import glob, pickle, pdb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import ion

hfilt = 'psw'
pscale = 30.

#pickledir = '/data-2/cross/helms/pickle/nomask/pickle_nvss_' + hfilt + '/'
pickledir = '/data-2/cross/helms/pickle/check/'
dd_pickle = glob.glob(pickledir + 'DD*.pickle')
rr_pickle = glob.glob(pickledir + 'RR*.pickle')
dr_pickle = glob.glob(pickledir + 'DR*.pickle')

rlen, dlen, drlen = len(rr_pickle), len(dd_pickle), len(dr_pickle)
if dlen != drlen:
    print('What the fuuuuck?')
    raise Exception(rlen)

goto = dlen
goto_RR = rlen


#goto = rlen if rlen < dlen else dlen
print('Using %i data files.' %goto)
print('Using %i random files.' %goto_RR)


# Check the binning and params.
for i, dpick in enumerate([dd_pickle[0],rr_pickle[0], dr_pickle[0]]):
    thisp = pickle.load(open(dpick,'r'))
    nbins = len(thisp['r'])
    if i == 0: 
        print 'DD object has these keys:' 
    elif i == 1: 
        print 'RR object has these keys:'
    else:
        print 'DR object has these keys:'

    print thisp.keys()



# Compute w(theta) from DD and RR pickles.
dd_sum = []#np.zeros(goto)
rr_sum = []#np.zeros(goto)
dr_sum = []
dd_n = []
rr_n = []
dr_n = []
for i in range(goto):
    dname = pickledir + 'DD_'+str(i)+'.pickle'
    drname= pickledir + 'DR_'+str(i)+'.pickle'

    this_dd = pickle.load(open(dname,'r'))
    this_dr = pickle.load(open(drname,'r'))

    dd_sum.append(this_dd['sum'])
    dr_sum.append(this_dr['sum'])

    #dd_n.append(this_dd['n'])
    #dr_n.append(this_dr['n'])

for i n range(goto_RR):
    rname = pickledir + 'RR_'+str(i)+'.pickle'
    this_rr = pickle.load(open(rname,'r'))
    rr_sum.append(this_rr['sum'])
    #rr_n.append(this_rr['n'])


DD = zip(*dd_sum)
RR = zip(*rr_sum)
DR = zip(*dr_sum)
#DD_n = zip(*dd_n)
#RR_n = zip(*rr_n)
#DR_n = zip(*dr_n)

DDi, RRi, DRi = [], [], []

fgal, fran = [], [] # Ngal, Nrandom factors for each
#for tRR, tDD, tDR, tRR_n, tDD_n, tDR_n in zip(RR, DD, DR, RR_n, DD_n, DR_n):
for tDD, tDR in zip(DD, DR):
    #n1 = np.asarray(tDD_n).astype(int)
    #n2 = np.asarray(tRR_n).astype(int)
    #n3 = np.asarray(tDR_n).astype(int)
    
    s1 = np.asarray(tDD)
    s2 = np.asarray(tRR)

    #g = np.where((n1 != 0) & (n2 != 0) & (n3 != 0))[0]
    #g1 = np.where(n1 != 0)[0]
    #g2 = np.where((n2 != 0) & (n3 != 0))[0]

    #g1, g2 = g, g
    #pdb.set_trace()
    DD_avg = s1#[g]#/n1[g1]#[g]/n1[g]
    RR_avg = s2#[g]#/n2[g2]#/n2[g]
    DR_avg = s3#[g]#/n3[g2]#/n3[g] 
    
    #*** need to normalize by (Ng/Nr)^2, so two if statements (g1 (for dd !=0), 
    #and g2 (RR,DR != 0 inclusive).
#    Ngal = len(g1)
#    Nran = len(g2)
#    fgal.append(Ngal)
#    fran.append(Nran)
    #pdb.set_trace()
    
    #keepers = np.where(~np.isnan(DD_avg) & ~np.isnan(RR_avg) & ~np.isnan(DR_avg))

    DDi_ = np.mean(DD_avg)#[keepers])#~np.isnan(DD_avg)])
    DRi_ = np.mean(DR_avg)#[keepers])#~np.isnan(DR_avg)])

    DDi.append(DDi_)
    DRi.append(DRi_)


for tRR in zip(RR):
    s3 = np.asarray(tDR)
    RRi_ = np.mean(RR_avg)#[keepers])#~np.isnan(RR_avg)])
    RRi.append(RRi_)


DD = np.asarray(DDi)
RR = np.asarray(RRi)
DR = np.asarray(DRi)

fuckit = np.where(~np.isnan(DD) & ~np.isnan(RR) & ~np.isnan(DR))
if len(fuckit) <= 0:
    print 'fucking NaNs'
    pdb.set_trace()


#pdb.set_trace()
#len(data_R) * 1. / len(data)
factor = np.asarray(fran).astype(float) / np.asarray(fgal).astype(float)    
#factor = factor.mean()**2 
#factor = 1.
#corr = (factor **2 * DD - 2 * factor * DR + RR) / RR 
corr = ( DD - 2 * DR + RR) / RR 
#corr = factor.mean()**2*DD/RR - 1.

#Nr, Nd = 1., 1.
#corr = Nr/Nd * DD/RR - 1   # DD/DR - 1 wtf
theta = thisp['r'] * 30./3600
print factor

plt.ion()
plt.figure()

plt.xscale('log'); plt.yscale('log')
#plt.errorbar(r, wthetas, dw,lw=3,fmt='o',color='gray')

plt.plot(theta, corr, 'o',color='gray',lw=3)
plt.xlabel(r'$\theta\ (degrees)$', fontsize=18)
plt.ylabel(r'$w(\theta)$', fontsize=18)
plt.xlim([.01,10]); plt.ylim([1e-4,1])

#plt.ticklabel_format(style='sci', axis='y', scilimits=(1,4))
plt.rcParams.update({'font.size': 11})
plt.gcf().subplots_adjust(bottom=0.15)

xx=pickle.load(open("check/xx.pkl",'r'))
yy=pickle.load(open("check/yy.pkl",'r'))
er=pickle.load(open("check/er.pkl",'r'))
plt.errorbar(xx,yy,er,lw=3,fmt='o',color='red')


#newx, newy = np.loadtxt('/home/ketron/cross/check/sdss_digitized.txt',unpack=True)
#plt.plot(newx, newy, 'o', color = 'blue')
pdb.set_trace()

#plt.errorbar(r,wtheta,yerr=dw)

