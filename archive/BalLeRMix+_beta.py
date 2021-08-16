#Jul, 2020: Expanded xGrid for positive selection
#Feb 20, 2020: Fixed the bug in B2maf.
#Jul 22, 2019: Modified the published BallerMix script for beta-binomial distribution, updates include:
    ##- replace optparse with argparse
    ##- use more numpy/scipy
import sys,argparse
from math import log,exp,floor # natural log
from datetime import datetime
from scipy.special import loggamma, betaln

def getAd(alpha):
    return(-log(alpha))

''' Read genome-wide substitution:polymorphism ratio
#format: N sub poly; no header
#0 for sub, 1 for poly
'''
def getPolyDen(spectfile):
    g=[]; G=[]
    with open(spectfile,'r') as spect:
        l=next(spect)
        l=l.strip().split('\t')
        N=int(l[0]); s=float(l[1]); p=float(l[2])
        print(('Substitutions: %s ; polymorphisms: %s' %(s,p)))
    try:
        assert s+p == 1
    except:
        print(('s+p= %s + %s = %s != 1' % (s,p,s+p)))
        sys.exit()
    g={(0,N):s,(1,N):p}
    G={(0,N):log(s), (1,N):log(p)}
    return(g,G,N)

'''Read unormalized neutral spect
##format: x n freq; no header
'''
def getLogSpect(spectfile,nosub,MAF):
    g={}; G={}; N=[]; checksum=0.
    with open(spectfile,'r') as spect:
        for l in spect:
            l=l.strip().split("\t")
            x=int(l[0]); n=int(l[1])
            f=float(l[2])
            g[(x,n)]=f
            G[(x,n)]=log(f)
            checksum += f
            N.append(n)
    N = list(set(N))
    print(('Total prob is %f in the spectrum file.' %(checksum)))
    if len(N) >= 2:
        print("Current implementation only supports uniform sample size. Please try again.")
        sys.exit()
    else:
        N=N[0]
    #Make sure there's no missing keys:
    fold=False
    for k in range(N+1):
        if (k,N) not in g:
            g[(k,N)] = 0
        elif MAF and k > N/2:
            fold = True
    #if fold;     
    if MAF and fold:
        fg = {}; fG = {}; checksum=0
        for i in range(N/2):
            if nosub and i==0:
                print('Skip substitutions')
                continue
            fg[(i,N)] = g[(i,N)] + g[(N-i,N)]
            try:
                fG[(i,N)] = log(fg[(i,N)])
            except:
                print("fg:", i, N, fg)
            checksum += fg[(i,N)]
        if N%2 == 0:
            fg[(N/2,N)] = g[(N/2,N)]
            fG[(N/2,N)] = log(fg[(N/2,N)])
            checksum += fg[(N/2,N)]
        print(('After folding, total probability is %s.'%(checksum)))
        return(fg,fG,N)
    #if not folded:
    return(g,G,N) #spect, logspect, size

''' Precompute binomial coefficients for sample size N
#log of Binom(n choose k)
# B = N!/k!(N-k)!
# logB = sum_{n-k+1}^n{logi} - sum_1^{n-k}{logi}
'''
def logBinCoeff(N):
    logIs = [log(i) for i in range(1,N+1)]
    logBs={}
    logBs[0]=0 ; logBs[N]=0
    for k in range(1,N):
        logBs[k] = sum(logIs[N-k:])-sum(logIs[:k])
    return logBs


''' Generate look-up table for beta function
# beta(a,b) = gamma(a)*gamma(b) / gamma(a+b)
# logbeta(k+a, N-k+b) = loggamma(k+a) + loggamma(N-k+b) - loggamma(a+b+N)
'''
def getab(x,a):
    b = a/x -a 
    return(a,b)

def logBeta(N,a,b):
    lnBeta = {}
    for k in range(N+1):
        #lnBeta[(k,N,a,b)] = logA[k] + logB[N-k] - logABN
        lnBeta[(k,N,a,b)] = betaln(a+k, b+N-k)
    return(lnBeta)
   

''' Precompute beta-binomial probability for sample size N, chance x for B2, B0, and B1
# logP[k] = logBs[k] + log(beta(k+a,N-k+b) ) - log(beta(a,b))
# Ancestral sites (k=0) not considered. Normalize probability across k=1...N
# For each x, return a dash table of k:prob[k]
'''
def getBetaBinom(N,x,a,logBs,nofreq,nosub):
    a,b = getab(x,a)
    logBab = betaln(a,b)
    if nofreq:
        logP0 = betaln(a,b+N) - logBab
        norm = log( 1-exp(logP0) )
        sub = exp( betaln(a+N,b) - logBab - norm)
        poly = 1-sub
        return(poly,sub)
    else:
        lgBeta = logBeta(N,a,b)
        logProb = {}; Prob={}
        P0 = exp(betaln(a,b+N) - logBab)
        if nosub:
            Pn = exp(betaln(a+N,b) -logBab)
            norm = log(1-P0-Pn)
        else:
            norm = log(1-P0)
        if nosub:
            Ncount = N
        else:
            Ncount = N+1
        for k in range(1,Ncount): #make sure to exclude k=0
            logProb[k] = logBs[k] + lgBeta[(k,N,a,b)] - logBab - norm
            Prob[k] = exp(logProb[k])
            try:
                assert logProb[k] < 0
            except:
                print(('logProb[k] >= 0; k=%s, logBs[k]=%s, logProb[k]=%s, lgBeta[(k,N,a,b)]=%s, logBab=%s, norm=%s, a=%s, b=%s, P0=%s, N=%s, 1-P0=%s, log(1-P0)=%s' % ( k,logBs[k], logProb[k], lgBeta[(k,N,a,b)], logBab, norm, a, b, P0, N, 1-P0, log(1-P0) )))
        return(Prob)#,logProb

''' Return folded probability distribution for B0maf or B2maf'''
def foldBinom(Prob,N,nosub):
    foldedProb={}
    #print Prob.keys()
    for k in range(N/2):
        if nosub and k==0:
            #print('Skip substitutions')
            continue
        elif k not in Prob:
            Prob[k] = 0
            #print('Assign f[%s] = 0'%(k))
            foldedProb[k] = Prob[k] + Prob[N-k]
        elif N-k not in Prob:
            #print('Assign f[%s] = 0'%(N-k))
            foldedProb[k] = Prob[k] 
        else:
            foldedProb[k] = Prob[k] + Prob[N-k]
        try:
            assert foldedProb[k] < 1
        except:
            print(Prob[k] , Prob[N-k])
    if N%2 == 0:
        foldedProb[N/2] = Prob[N/2]
    #print(foldedProb)
    return(foldedProb)


''' Generate the grids of x and A values to optimize over. The grid of alpha is fixed.'''
def getGrids(x, abeta, bal, pos, MAF, seqA, listA):
    #define the default list of mixing proportion alpha:
    aGrid = [0,1e-8,1e-7,1e-6,1e-5]+[m*1e-4 for m in range(1,10001)]
    #generate sorted AdGrid based on the dense grid of alpha:
    AdGrid = [-log(1e-32)]+[ getAd(a) for a in aGrid[1:] ] 
    assert len(aGrid) == len(AdGrid)
    #get xGrid
    if x:
        xGrid = [float(x)]
    else:
        xGrid = [.05*i for i in range(1,11)]
    #print(xGrid)
    #get abetaGrid
    #print abeta
    if abeta > 0 :
        abetaGrid = [float(abeta)]
    elif bal:
        abetaGrid = [i for i in range(1,10)] + [5*i for i in range(1,20)] + [10*i for i in range(10,21)] + [300,500, 1e3, 1e4,1e6,1e9]
    elif pos:
        abetaGrid = [0.001,0.01, 0.05,0.1,0.2,0.5,0.8]
        xGrid = [.1*i for i in range(1,11)]
    else:
        abetaGrid = [0.001,0.01, 0.05,0.1,0.2,0.5,0.8] + [i for i in range(1,10)] + [5*i for i in range(1,20)] + [10*i for i in range(10,21)] + [300,500,1e3, 1e4,1e6,1e9]
    #print abetaGrid
    if not seqA and not listA:
        AGrid = [100*i for i in range(1,12)] + [200*i for i in range(6,13)] + [500*i for i in range(5,10)] + [1000*i for i in range(5,11)] + [1e6,1e8]
        #The grid is: 100-1100, 1200-2400, 2500-4000, 5000-9000, 1e4,1e6,1e8
        #total of 29 values
    elif listA:
        AGrid = [float(x) for x in listA.split(',')]
    else:
        Amin,Amax,Astep = [float(x) for x in seqA.split(',')]
        n=(Amax-Amin)/Astep
        AGrid=[Amin+Atep*i for i in range(n+1)]
    return(xGrid,AGrid,AdGrid,aGrid,abetaGrid)

''' Generate the look-up table for g(k) and f(k,d; x,a,A) distributions '''
#By default, the bin width of alpha values is 1e-4
def initialize(spectfile,xGrid,abetaGrid,aGrid,MAF,nofreq,nosub):
    Fx={}
    if nofreq: #B1
        print('Only consider polymorphism density...')
        Spec,logSpec,N = getPolyDen(spectfile)
        print(('Sample size: %s'%(N)))
        logBs = logBinCoeff(N)
        #Spec[0]:subs, Spec[1]:poly
        Fx={}
        for X in xGrid:
            Fx[X] = {}
            for abeta in abetaGrid:
                Fx[X][abeta] = {}
                poly,sub = getBetaBinom(N,X,abeta,logBs, nofreq, nosub)
                Hx = [sub, poly]
                for a in aGrid:
                    Fx[X][abeta][a] = {}
                    for k in [0,1]: #0 for sub, 1 for poly
                        gx = Spec[(k,N)]
                        fx = a*Hx[k] + (1-a)*gx
                        Fx[X][abeta][a][k] = log(fx)
    else: 
        Spec,logSpec,N = getLogSpect(spectfile,nosub,MAF)
        print(('Sample size: %s'%(N)))
        logBs = logBinCoeff(N)
        if not MAF:
            print('Using derived allele frequency...')
        else:
            print('Using minor allele frequency...')
        Prob={}; oriProb={}; 
        for x in xGrid:
            Fx[x]={}
            for abeta in abetaGrid:
                #if abeta < 1 and x not in [0.05,0.25,0.45]:
                #    continue
                Fx[x][abeta]={}
                oriProb[x] = getBetaBinom(N,x,abeta,logBs,nofreq,nosub)
                oriProb[1-x] = getBetaBinom(N,1-x,abeta,logBs,nofreq,nosub)
                if MAF:
                    Prob[x] = foldBinom(oriProb[x],N,nosub)
                    Prob[1-x] = foldBinom(oriProb[1-x],N,nosub)
                else:
                    Prob[x] = oriProb[x]
                    Prob[1-x] = oriProb[1-x]
                for a in aGrid:
                    Fx[x][abeta][a]={}
                    #set the list of k
                    if MAF: #B2maf
                        klist = list(range(N/2+1))
                        if nosub:
                            klist = list(range(1,N/2+1))
                    else: #B
                        klist = list(range(1,N+1))
                        if nosub:#remove k=N
                            klist = list(range(1,N))
                    for k in klist:
                        gx = Spec[(k,N)]
                        fx = a*(.5*Prob[x][k]+.5*Prob[1-x][k]) + (1-a)*gx
                        #print('%g\t%f\t%f\t%f'%(k,a,gx,fx))
                        try:
                            assert fx <1
                        except:
                            print('fx=%s; abeta=%s; gx=%s; a=%s; k=%s;\nProb[x][k]=%s; Prob[1-x][k]=%s'%(fx, abeta,gx,a,k,Prob[x][k],Prob[1-x][k]))
                            sys.exit()
                        try:
                            Fx[x][abeta][a][k] = log(fx)
                        except:
                            print(fx, abeta, a, k)
                            sys.exit()
    print(('%s Initialization finished.\n' % (datetime.now())))
    #print Prob[0.5]
    #print sum(Prob[0.5].values())
    #print(Fx[0.5][1e3])
    #print Fx.keys()
    #print Fx[0.5].keys()
    #print Fx[0.5][1e6].keys()
    #total=0
    #for k,v in Fx[0.05][0.001][0.5].items():
    #    print k,v,exp(v)
    #    total += exp(v)
    #print total
    #print logSpec
    return(Spec,logSpec,Fx,N)

''' Return parsed data
# all sitePos are genetic positions (by cM)
# default value of Rrate is 1e-6 cM/site
# Input format: physPos, genPos, x, n
'''
def readInput(infile,nosub,nofreq,MAF,phys=False,Rrate=1e-6):
    phys_pos = []
    Ks=[]; Ns = []; pos_list=[]; numSites=0 
    postype=1-int(phys) #index; 0 is physical position, 1 is genetic position
    with open(infile,'r') as sites:
        l=next(sites)#skipping header
        if nosub:
            #make sure to only read frequencies
            for l in sites:
                l=l.strip().split('\t')
                if float(l[2])/float(l[3]) not in [0.,1.]:
                    numSites += 1
                    physPos,k,n = [float(l[0]),int(l[2]),int(l[3])]
                    #genPos = float(l[1])
                    physPos = int(float(l[0]))
                    sitepos = float(l[postype])*phys*Rrate + float(l[postype])*(1-int(phys))
                    k = MAF*min(k,n-k) + (1-MAF)*k
                    Ks.append(k); Ns.append(n); pos_list.append(sitepos)
                    phys_pos.append(physPos)
        elif nofreq: # Take all none N cases as 1
            for l in sites:
                l=l.strip().split('\t')
                numSites += 1
                physPos,k,n = [float(l[0]),int(l[2]),int(l[3])]
                #genPos = float(l[1])
                physPos = int(float(l[0]))
                sitepos = float(l[postype])*phys*Rrate + float(l[postype])*(1-int(phys))
                k = (k!=n) #k=1 if not sub, k=0 if sub
                Ks.append(k); Ns.append(n); pos_list.append(sitepos)
                phys_pos.append(physPos)
        else:
            for l in sites:
                l=l.strip().split('\t')
                numSites += 1
                try:
                    physPos,k,n = [float(l[0]),int(l[2]),int(l[3])]
                except:
                    print(l, numSites)
                    sys.exit()
                #genPos = float(l[1])
                physPos = int(float(l[0]))
                sitepos = float(l[postype])*phys*Rrate + float(l[postype])*(1-int(phys))
                k = MAF*min(k,n-k) + (1-MAF)*k
                Ks.append(k); Ns.append(n); pos_list.append(sitepos)
                phys_pos.append(physPos)
        if(MAF):
            print((Ks[:5]))
    return(phys_pos,pos_list,Ks,Ns,numSites)

''' Return log likelihood for the given window
# testSite is genetic position (in recomb unit)
# logP_neut = logSpect[(k,N)]; logP_alt = Fx[x][a][k]'''
def calcBaller(pos_list,Ks,Ns, testSite, testSite_i, logSpect, Fx, xGrid, abetaGrid, AGrid, AdGrid, aGrid):
    #note that the test site is included in Ks
    #Optimize over xGrid and AGrid
    L = len(AdGrid) ; numSites = len(pos_list)
    #testSite = pos_list[testSite_i]
    Tmax = [-100,0,-1,0] #LR, x, abeta, A
    nSites=0
    # Go through all A values
    for A in AGrid:
        # i indexes positions. pos_list ascending
        # c indexes Ad value in the pre-computed descendingly sorted list
        i_index=[] ; c_index=[]
        # going through sites from the center to outward
        i = max(testSite_i-1,0) ; c = len(AdGrid)-1
        # Leftward:
        while i >= 0 and c > 0 and testSite_i != 0:
            pos = pos_list[i]
            dist = testSite - pos
            Ad = A*dist
            while AdGrid[c] < Ad and c >= 0:
                c -= 1
            #now AdGrid[c] >= Ad or c=0
            if c >=0:
                i_index.append(i)
                c_index.append(c)
            i -= 1
        # now i==0. Starting rightward from center
        #i=min(testSite_i+1,pos_list[-1]) 
        i = min(testSite_i+1,numSites) ; c = len(AdGrid)-1
        while i < len(pos_list) and c > 0:
            pos = pos_list[i]
            dist = pos - testSite
            Ad = A*dist
            while AdGrid[c] < Ad and c >=0:
                c -= 1
            if c >= 0:
                i_index.append(i)
                c_index.append(c)
            i += 1
        # if noCenter, consider testSite_i too
        if testSite != pos_list[testSite_i]:
            pos = pos_list[testSite_i]
            dist = abs(pos - testSite)
            Ad = A*dist
            while AdGrid[c] < Ad and c >=0:
                c -= 1
            if c>=0:
                i_index.append(testSite_i)
                c_index.append(c)
        #go through the grid of x
        #note that this grid only spans 0-0.5 
        for x in xGrid:
            for abeta in abetaGrid:
                #if abeta < 1 and x not in [0.05,0.25,0.45]:
                #    continue
                La=0; L0=0
                #go through all the sites
                for j in range(len(i_index)):
                    i = i_index[j] #index in the pos_list
                    c = c_index[j]
                    alpha = aGrid[c]
                    k = Ks[i] ; N = Ns[i]
                    try:
                        La += Fx[x][abeta][alpha][k]
                    except(KeyError):
                        #print x, abeta, alpha, k
                        print(( list(Fx[x].keys()) ))
                        sys.exit()
                    L0 += logSpect[(k,N)]
                #get LR
                T = 2*(La - L0)
                if T >=Tmax[0]:
                    Tmax = [T,x,abeta,A]
                    nSites = len(i_index)
    return(Tmax,nSites)


'''Perform the scan with fixed window size of w (bp), with step size s (bp).'''
def scan_fixSize_noCenter(xGrid,abetaGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,Rrate,w,s,MAF=False):
    #w = float(r) 
    print(("writing output to %s" % (outfile)))
    with open(outfile,'w') as scores:
        scores.write('physPos\tgenPos\tLR\tx_hat\ts_hat\tA_hat\tnSites\n')#\tnumSites
        start = int( floor(2*float(phys_pos[0])/w)*(w/2) )
        end = start + s ; midpos = start + s/2
        start_i=0; end_i=0; pos_i=0
        while midpos <= phys_pos[-1]:
            #define the window
            while phys_pos[start_i] < start:
                start_i +=1 
            while phys_pos[pos_i] < midpos:
                pos_i += 1
            while end_i < numSites:
                if phys_pos[end_i] < end:
                    end_i += 1
                else:
                    break
            if start_i >= end_i:
                scores.write('%g\t%s\t0\tNA\tNA\t0\n' % (midpos,midpos*Rrate))#
                start+=s; midpos+=s; end+=s
            else:
                #Window = Ks[start_i:end_i+1]
                WinKs = Ks[start_i:pos_i] + [0] + Ks[pos_i:end_i+1]
                WinNs = Ns[start_i:pos_i] + Ns[pos_i] + Ns[pos_i:end_i+1]
                winPosList = pos_list[start_i:pos_i] + [midpos*Rrate] + pos_list[pos_i:end_i+1]
                #calcBaller args: pos_list,Ks,Ns, testSite, testSite_i, logSpect, biFx, xGrid, AGrid, AdGrid, aGrid,MAF
                #pos_list and phys_pos has matching indice
                Tmax,winSites = calcBaller(winPosList,WinKs, WinNs, midpos*Rrate, pos_i-start_i, logSpec, Fx, xGrid,abetaGrid, AGrid, AdGrid, aGrid)
                scores.write('%g\t%s\t%s\t%s\t%g\t%g\t%d\n' % (midpos, midpos*Rrate, Tmax[0], Tmax[1], Tmax[2], Tmax[3], winSites))#
                #read in, take next step
                start+=s; midpos+=s; end+=s 
    scores.close()
    return(0)

'''Perform the scan with fixed window size of w (bp), centered on every s informative sites'''
def scan_fixSize_siteCenter(xGrid,abetaGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,Rrate,w,s,MAF=False):
    #w = float(r) 
    print(("writing output to %s" % (outfile)))
    with open(outfile,'w') as scores:
        scores.write('physPos\tgenPos\tLR\tx_hat\ts_hat\tA_hat\tnSites\n')#\tnumSites
        i=0; 
        start_i=0; end_i=0
        while i < numSites:
            testSite=phys_pos[i]
            start =  max(0, testSite-r/2); end = min(testSite+r/2,phys_pos[-1])
            while phys_pos[start_i] < start:
                start_i += 1
            while end_i < numSites:
                if phys_pos[end_i] < end:
                    end_i += 1
                else:
                    break
            assert end_i >= start_i
            end_i = min(end_i, numSites-1)
            #Window = Ks[start_i:end_i+1]
            Tmax,winSites = calcBaller(pos_list[start_i:end_i+1], Ks[start_i:end_i+1], Ns[start_i:end_i+1], pos_list[i], i-start_i , logSpec, Fx, xGrid,abetaGrid, AGrid, AdGrid, aGrid)
            scores.write('%g\t%s\t%s\t%s\t%g\t%g\t%d\n' % (phys_pos[i], pos_list[i], Tmax[0], Tmax[1], Tmax[2], Tmax[3], winSites))#
            i+=int(s)
    scores.close()
    return(0)

'''Perform the scan with site-based window, with s sites on either side'''
def scan_siteBased(xGrid,abetaGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,phys,Rrate,r,s=1,MAF=False):
    print(("writing output to %s" % (outfile)))
    with open(outfile,'w') as scores:
        scores.write('physPos\tgenPos\tLR\tx_hat\ts_hat\tA_hat\tnSites\n')#\tnumSites
        i=0
        while i < numSites:
            testSite = pos_list[i]
            start_i = max(0, i-r) ; end_i = min(numSites,i+r+1)
            Tmax,winSites = calcBaller(pos_list[start_i:end_i], Ks[start_i:end_i], Ns[start_i:end_i], testSite, i-start_i , logSpec, Fx, xGrid,abetaGrid, AGrid, AdGrid, aGrid)
            scores.write('%g\t%s\t%s\t%s\t%g\t%g\t%d\n' % (phys_pos[i], pos_list[i], Tmax[0], Tmax[1], Tmax[2], Tmax[3], winSites))#
            i+=s
    scores.close()
    return(0)

def scan_alpha(xGrid,abetaGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec,logSpec,Fx,N,s=1,MAF=False):
    print(("writing output to %s" % (outfile)))
    with open(outfile,'w') as scores:
        scores.write('physPos\tgenPos\tLR\tx_hat\ts_hat\tA_hat\tnSites\n')#\tnumSites\tLa\tL0
        i=0
        while i < numSites:
            testSite = pos_list[int(i)]
            Tmax,winSites = calcBaller(pos_list, Ks, Ns, testSite, int(i), logSpec, Fx, xGrid, abetaGrid, AGrid, AdGrid, aGrid)
            scores.write('%g\t%s\t%s\t%s\t%g\t%g\t%d\n' % (phys_pos[i], pos_list[i], Tmax[0], Tmax[1], Tmax[2], Tmax[3], winSites))#\t%s\t%s,Tmax[4],Tmax[5]
            i+=int(s)
    scores.close()
    return(0)


def scan(xGrid,abetaGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,size=False,r=0,s=1,phys=False,nofreq=False,MAF=False,noCenter=False,Rrate=1e-6):
    if size:
        print('You\'ve chosen to fix the scan window size.')
        if r==0:
            print('Please set a window width in bp with \"-w\" or \"--window\" command.')
            sys.exit()
        if not phys:
            print(('Please make sure to use physical positions as coordinates if fixed-length windows are chosen. Scan will continue with physical positions with a rec rate of %f cM/nt.'%(Rrate)))
            phys = True
            #sys.exit()
        if noCenter:
            w = float(r) 
            print(('Computing LR on %.3f kb windows on every %s bp. Using physical positions by default.' % (w/1e3, s)))
            scan_fixSize_noCenter(xGrid,abetaGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec,logSpec,Fx,N,Rrate,w,s)
        else: #site-centered
            w = float(r) 
            print(('Computing LR on %.2f kb windows on every %g informative sites. Using physical positions by default.' % (w/1e3,s)))
            scan_fixSize_siteCenter(xGrid,abetaGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,Rrate,w,s)
    # When fixSize == False, and radius (-r) provided. Scan with fixed number of sites
    elif r != 0: 
        print(('Computing LR on every %s site/s, with %s informative sites on either side.' % (s, r)))
        scan_siteBased(xGrid,abetaGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,phys,Rrate,r,s)
    #window size not given 
    #then use all data (but the test site)
    else:
        print(('Computing LR on every %s site/s, using all the data with alpha >= 1e-8.' % (s)))
        scan_alpha(xGrid,abetaGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec,logSpec,Fx,N,s)
    print((str(datetime.now())+'. Scan finished.'))


'''#generate the config file given the concatenated input'''
def getConfig(infile,configfile):
    Config={}; numSites=0# N: [s,p]
    with open(infile,'r') as sites:
        l=next(sites)#skip the header by default
        for l in sites:
            x,n = [int(x) for x in l.strip().split('\t')[2:] ]
            if x==0:
                print('Please make sure the input has derived allele frequency. Sites with 0 observed allele count (k=0) will be ignored.\n')
                continue
            if n not in Config:
                Config[n] = [0,0]
            Config[n][0] += int(x==n)
            Config[n][1] += 1-int(x==n)
            numSites+=1
    sizes = sorted(Config.keys())
    with open(configfile,'w') as config:
        for N in sizes:
            config.write('%s\t%s\t%s\n' % ( N, Config[N][0]/float(numSites) , Config[N][1]/float(numSites) ))
    sites.close(); config.close()
    print('Done')

'''#generate spectrum file given the concatenated input'''
def getSpect(infile,spectfile,nosub=False,MAF=False):
    Spect={}; numSites=0
    if nosub:
        print('Generating spectrum for all polymorphic sites. Substitutions (x=n or x=0) won\'t be considered.')
        with open(infile,'r') as sites:
            l=next(sites)
            for l in sites:
                (x,n)=[ int(i) for i in l.strip().split('\t')[2:] ]
                if MAF:
                    x = min(x, n-x)
                if x == 0 or x == n:
                    continue
                if (x,n) in Spect:
                    Spect[(x,n)] += 1
                else:
                    Spect[(x,n)] = 1
                numSites+=1
    else:
        with open(infile,'r') as sites:
            l=next(sites)
            for l in sites:
                (x,n)=[ int(i) for i in l.strip().split('\t')[2:] ]
                if MAF:
                    x = min(x, n-x)
                elif x==0:
                    print('Please make sure the input has derived allele frequency. Sites with 0 observed allele count (k=0) should not be included.\n')
                    sys.exit()

                if (x,n) in Spect:
                    Spect[(x,n)] += 1
                else:
                    Spect[(x,n)] = 1
                numSites+=1
    #write out
    pairs = sorted(Spect.keys())
    with open(spectfile,'w') as spec:
        for x,n in pairs:
            spec.write('%s\t%s\t%s\n' % (x,n,float(Spect[(x,n)])/float(numSites) ))
    sites.close(); spec.close()
    print('Done.')

#main function to scan through the input
#only work with minor allele frequency as equilibrium frequency
def main(): 
    #parsing arguments
    parser = argparse.ArgumentParser()#usage='python {} -i <input file> -o <output file> --spect <spect/config file> [--help] [--nofreq] [--nosub] [--MAF] [--getSpect] [--getConfig] [--fixSize] [--physPos] [--rec <recomb rate>] [-w <window size>] [--noCenter] [-s <step size>] [--fixX <x>] [--fixAlpha] [--findBal] [--findPos] [--rangeA <min,max,step>] [--listA <A1,A2,..,Ak>]'.format(sys.argv[0])
    parser.add_argument('-i','--input', dest='infile', help = 'Path and name of your input file.\n', required=True)
    parser.add_argument('-o','--output', dest='outfile', help = 'Path and name of your output file.\n')
    parser.add_argument('--spect',dest='spectfile', help = 'Path and name of the allele frequency spectrum file or configuration file.\n', required=True)
    parser.add_argument('--getSpect', dest='getSpec', action='store_true', default=False, help='Option to generate frequency spectrum file from the concatenated input file. Use \"-i\" and \"--spect\" commands to provide names and paths to input and output files, respectively. Indicate the input type with \"--MAF\" and/or \"--nosub\".\n')
    parser.add_argument('--getConfig', dest='getConfig', action='store_true', default=False, help='Option to generate configuration file from the concatenated input file. Use \"-i\" and \"--spect\" commands to provide names and paths to input and output files, respectively.\n\n')

    parser.add_argument('--nofreq', dest='nofreq', action='store_true', default=False, help = 'Option to ignore allele frequency information. All polymorphic sites will be considered as equivalent.')
    parser.add_argument('--nosub', dest='nosub', action='store_true', default=False, help = 'Option to not include substitution in input data.')
    parser.add_argument('--MAF', dest='MAF', action='store_true', default=False, help = 'Option to use minor allele frequency, instead of polarized allele frequency. The latter is default.')
    parser.add_argument('--findBal',dest='bal', action='store_true', default=False, help="Option to only look for footprints of balancing selection.\n")
    parser.add_argument('--findPos',dest='pos', action='store_true', default=False, help="Option to only look for footprints of positive selection.\n")
    parser.add_argument('--physPos', action='store_true', dest = 'phys', default = False, help = 'Option to use physical positions instead of genetic positions (in cM). Default is using genetic positions.\n')
    parser.add_argument('--rec', dest='Rrate', default = 1e-6 , help='The uniform recombination rate in cM/nt. Default value is 1e-6 cM/nt. Only useful when choose to use physical positions as coordinates.\n\n')

    parser.add_argument('--fixSize', action='store_true', dest = 'size', default = False, help = 'Option to fix the size of scanning windows. When true, provide the length of window in neucleotide (nt) with \"-w\" or \"--window\" command.\n')
    parser.add_argument('-w','--window', dest='w', type = int, default=0, help='Number of sites flanking the test locus on either side. When choose to fix window size (\"--fixSize\"), input the length of window in bp.\n')
    parser.add_argument('--noCenter', action='store_true', dest='noCenter', default=False, help = 'Option to have the scanning windows not centered on informative sites. Require that the window size (\"-w\") in physical positions (\"--physPos\") is provided. Default is True.\n')
    parser.add_argument('-s','--step', dest='step', type = float, default=1, help='Step size in bp (when using \"--noCenter\") or the number of informative sites. Default value is one site or one nucleotide.\n\n')

    parser.add_argument('--fixX', dest='x', help='Option to fix the presumed equilibrium frequency.\n')
    parser.add_argument('--fixAlpha', dest='abeta', type=float, default=-1, help='Option to fix the alpha parameter in the beta-binomial distribution.\n')
    parser.add_argument('--rangeA', dest='seqA', help='Range of the values of the linkage parameter A to optimize over. Format should follow <Amin>,<Amax>,<Astep> with no space around commas.\n')
    parser.add_argument('--listA', dest='listA', help='Manually provide a list of A values to optimize over. Please separate the values with comma, no space.\n')

    
    #if len(sys.argv[1:]) == 0:
    #    parser.print_help()
    #    sys.exit()

    opt = parser.parse_args(sys.argv[1:])
    #print opt.abeta
    if opt.getSpec:
        print('You\'ve chosen to generate site frequency spectrum...')
        print(('Concatenated input: %s \nSpectrum file: %s' % (opt.infile, opt.spectfile )))
        getSpect(opt.infile, opt.spectfile, opt.nosub, opt.MAF)
        sys.exit()
    elif opt.getConfig:
        print('You\'ve chosen to generate the substitution-polymorphism configuration...')
        print(('Concatenated input: %s \nConfiguration file: %s' % (opt.infile, opt.spectfile )))
        getConfig(opt.infile, opt.spectfile)
        sys.exit()

    #generate the grids to optimize over/with
    xGrid,AGrid,AdGrid,aGrid,abetaGrid = getGrids(opt.x, opt.abeta, opt.bal, opt.pos, opt.MAF, opt.seqA, opt.listA)
    #print abetaGrid
    print(('\nOptimizing over x= '+', '.join([str(x) for x in xGrid])))
    print(('\nOptimizing over alpha= '+', '.join([str(a) for a in abetaGrid])))
    print(('\nOptimizing over A= '+', '.join([str(x) for x in AGrid])))

    #initialization for a grid of x
    print(('\n%s. Initializing...'%(datetime.now())))
    (Spec, logSpec, Fx,N) = initialize(opt.spectfile,xGrid,abetaGrid,aGrid,opt.MAF,opt.nofreq,opt.nosub)#initialize(spectfile,xGrid,aGrid,MAF,nofreq,nosub)

    #start reading data
    print(("\n%s. Reading input file: %s" % (datetime.now(),opt.infile) ))
    (phys_pos,pos_list,Ks,Ns,numSites) = readInput(opt.infile,opt.nosub,opt.nofreq,opt.MAF,opt.phys,opt.Rrate)

    #finished reading file, start writing output
    print(("\n%s. Start computing likelihood raito..." %(datetime.now())))
    scan(xGrid,abetaGrid,AGrid,AdGrid,aGrid,opt.outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec, Fx,N, opt.size,opt.w,opt.step,opt.phys,opt.nofreq,opt.MAF,opt.noCenter,opt.Rrate)#, ProbsStar

    #pipeline finished
    print(('\n%s. Pipeline finished.'%(datetime.now())))



if __name__ == '__main__':
    main()
