'''
This script compute the B_1, B_2, and B_2maf statistics described in Cheng and DeGiorgio (2021) for jointly detecting both balancing selection and positive selection.
'''
import sys, os, re
from datetime import datetime

#module for reading input
class InputData:

	def __init__(self, infile, nofreq=False, MAF=False, minCount = 1, phys=False, Rrate=1e-6):
		self.position = []
		self.genPos = []
		self.count = []
		self.total = []
		self.numSites = 0 
		self.minCount = minCount
		self.Rrate = Rrate
		pos_type=1-int(phys) #index; 0 is physical position, 1 is genetic position

		#read input
		if nofreq:
			self.readPolyCalls(infile, pos_type, Rrate)
		else:
			self.readCounts(infile, pos_type, Rrate, MAF)
			if MAF:
				self.minCount = min([k for k in self.count if k > 0])
			else:
				self.minCount = min(self.count)

		self.sampSizes = set(self.total)
		self.count = np.array(self.count)
		self.total = np.array(self.total)
		self.genPos = np.array(self.genPos)

	#take all non-N cases as polymorphisms
	##k=1 if polymorphism, k=0 if substitution
	def readPolyCalls(self, infile, pos_type, Rrate):
		translate = False
		with open(infile, 'r') as sites:
			l = next(sites) #skip the header
			for l in sites:
				l = l.strip().split('\t')
				self.numSites += 1
				#in case physical positions carry decimal points
				physPos,k,n = [int(float(l[0])),int(l[2]),int(l[3])]

				#check for input format
				if not translate:
					try:
						assert k in set([0, 1])
					except:
						print( f'Input includes different variant counts despite choosing not to use allele frequencies (with --noFreq). All sites with counts smaller than substitutions will be considered as polymorphic. All sites with identical counts as sample sizes will be substitutions.')
						translate = True
						k = (k!=n)
				else:
					#k=1 if polymorphism, k=0 if substitution
					k = (k!=n) 

				#store physPos*Rrate (column 0) if use physical position, else store genetic position (column 1)
				sitepos = float(l[pos_type])*(1-pos_type)*Rrate + float(l[pos_type])*(pos_type)
				
				self.count.append(k)
				self.total.append(n)
				self.genPos.append(sitepos)
				self.position.append(physPos)
		#just to make sure it closes
		sites.close()

	def readCounts(self, infile, pos_type, Rrate, MAF):
		translate = False
		with open(infile , 'r') as sites:
			l = next(sites)
			for l in sites:
				l = l.strip().split('\t')
				self.numSites += 1
				#in case physical positions carry decimal points
				physPos,k,n = [int(float(l[0])),int(l[2]),int(l[3])]

				#check for B2maf
				if MAF and not translate:
					try:
						assert k <= n/2
					except:
						print(f'Input data includes non-MAF site/s (frequency >= 0.5) despite choosing to use B_2maf (with --MAF). These frequencies will be folded for following analyses.')
						translate = True
						k = n-k
				elif MAF and translate:
					k = min(k, n-k)
				#store physPos*Rrate (column 0) if use physical position, else store genetic position (column 1)
				sitepos = float(l[pos_type])*(1-pos_type)*Rrate + float(l[pos_type])*(pos_type)

				self.count.append(k)
				self.total.append(n)
				self.genPos.append(sitepos)
				self.position.append(physPos)
		#just to make sure it closes
		sites.close()


class Grids:

	def __init__ (self, x, abeta, bal, pos, seqA, listA) :
		#print(x, abeta, bal, pos, seqA, listA)
		#get _xGrid
		if x is not None:
			_xGrid = [float(x)]
		else:
			_xGrid = [.05*i for i in range(1,11)]
		#print(_xGrid)
		#get _abetaGrid
		if abeta is not None :
			try:
				_abetaGrid = [float(abeta)]
			except:
				print (f'The value for \"a\" provided ({abeta}) is not legitimate. Using the default grid instead.')
				_abetaGrid = [0.001,0.01, 0.05,0.1,0.2,0.5,0.8] + [i for i in range(1,10)] + [5*i for i in range(1,20)] + [10*i for i in range(10,21)] + [300,500,1e3, 1e4,1e6,1e9]
		elif bal:
			_abetaGrid = [i for i in range(1,10)] + [5*i for i in range(1,20)] + [10*i for i in range(10,21)] + [300,500, 1e3, 1e4,1e6,1e9]
		elif pos:
			_abetaGrid = [0.001,0.01, 0.05,0.1,0.2,0.5,0.8]
			_xGrid = [.1*i for i in range(1,11)]
		else:
			_abetaGrid = [0.001,0.01, 0.05,0.1,0.2,0.5,0.8] + [i for i in range(1,10)] + [5*i for i in range(1,20)] + [10*i for i in range(10,21)] + [300,500,1e3, 1e4,1e6,1e9]
		#print(f'Optimizing over alpha in {_abetaGrid}')

		#get _AGrid
		if not seqA and not listA:
			_AGrid = [100*i for i in range(1,12)] + [200*i for i in range(6,13)] + [500*i for i in range(5,10)] + [1000*i for i in range(5,11)] + [1e6,1e8]
			#The grid is: 100-1100, 1200-2400, 2500-4000, 5000-9000, 1e4,1e6,1e8
			#total of 29 values
		elif listA:
			_AGrid = [float(x) for x in listA.split(',')]
			#print(f'Optimizing over the grid for A: {_AGrid}.')
		else:
			Amin,Amax,Astep = [float(x) for x in seqA.split(',')]
			n=(Amax-Amin)/Astep
			_AGrid=[Amin+Atep*i for i in range(n+1)]
			#print(f'Optimizing over the grid for A: {_AGrid}.')
		self.x = _xGrid
		self.A = _AGrid
		self.abeta = _abetaGrid



#module for initialization
class NeutralSFS:

	#read spect for B2 or B2maf. For B2maf, check whether polarized frequency exist, and if so, fold.
	def readSpect (self, spectfile, MAF):
		g = {}; N = []
		checksum=0.
		with open(spectfile,'r') as spect:
			for l in spect:
				l=l.strip().split("\t")
				x=int(l[0]); n=int(l[1])
				f=float(l[2])
				if MAF:
					try:
						assert x < (n/2 + 1)
						g[(x,n)]=f
					except:
						print('You have indicated to use minor allele frequencies (--MAF) but provided SFS based on polarized allele frequency. This SFS will be folded accordingly.')
						if (n-x, n) in g:
							g[ (n-x, n)] += f
						else:
							g[ (n-x, n)] = f
				else:
					g[(x, n)] = f

				checksum += f
				N.append(n)
				if n not in self.sampProps:
					self.sampProps[n] = 0.
				self.sampProps[n] += f
		spect.close()
		#sanity check
		try:
			assert np.isclose(checksum, 1.)
		except:
			print(f'Fraction of sites do not add up to 1! Sum = {checksum}. Please double-check your inputs.' )
			sys.exit()

		self.spect = g ; self.sampSizes = set(N)
		#return(g, N)

	#@staticmethod
	def readConfig(self, spectfile):
		N = [] ; checksum = 0.
		with open(spectfile,'r') as spect:
			for l in spect:
				#l=next(spect) #only read first line
				l=l.strip().split('\t')
				n=int(l[0]); s=float(l[1]); p=float(l[2])
				print(('Substitutions: %s ; polymorphisms: %s' %(s,p)))
				checksum += (s+p)
				g={(0,n):s,(1,n):p}
				N.append(n)
				if n not in self.sampProps:
					self.sampProps[n] = 0
				self.sampProps[n] += (s+p)
		#in case
		spect.close()
		#sanity check
		try:
			assert checksum == 1.
		except:
			print(f'Fraction of sites do not add up to 1! Sum = {checksum}. Please double-check your inputs.' )
			sys.exit()
		
		self.spect = g ; self.sampSizes = set(N)
		#return(g, N)

	def __init__(self, spectfile, nofreq, MAF, nosub):
		self.spect = {}
		self.sampSizes = set()
		self.probs = []
		self.logProbs = []
		self.sampProps = {}
		self.propSizes = []

		#read helper file and get sample size
		if nofreq: #B1
			#read GW/neut
			self.readConfig(spectfile)

		else: #B2 and B2maf
			self.readSpect(spectfile, MAF, nosub)

	#get individual probs for poly-sub data
	def _get_one_neut_prob_B1(self, k, n):
		return ( self.spect[(k, n)] )

	#get individual probs for allele freq data
	def _get_one_neut_prob_B2s(self, k, n):
		return ( self.spect[(k, n)] )

	#get the proportion of sample sizes
	def _get_one_prop_samp_size(self, n):
		return ( self.sampProps[n] )

	#get a list of per-site neut probs 
	def get_neut_probs(self, InputData, nofreq, MAF):
		if self.probs == []:
			#sanity check
			try:
				assert InputData.sampSizes.issubset( self.sampSizes )
			except:
				print(f'Input data includes sample sizes not included in the helper file. Please double-check your inputs.')
				sys.exit()

			#vectorize
			if nofreq:
				self._get_neut_probs = np.vectorize( self._get_one_neut_prob_B1 )
			else:
				self._get_neut_probs = np.vectorize( self._get_one_neut_prob_B2s )

			self.probs = self._get_neut_probs( InputData.count, InputData.total )

		#sanity check
		assert len(self.probs) == InputData.numSites

		if self.logProbs == []:
			#get log for convenience
			self.logProbs = np.log( self.probs )

		#get frac of sample size for future use
		if self.propSizes == []:
			self._get_prop_samp_sizes = np.vectorize( self._get_one_prop_samp_size )
			self.propSizes = self._get_prop_samp_sizes( InputData.total )


import numpy as np
from scipy.stats import betabinom

class NormalizedBetaBinom:

	#convert (x, a) to (a, b)
	#compatible with numpy
	@staticmethod
	def _get_b(x,a):
		b = a/x -a 
		return(b)

	def __init__(self, InputData, Grids, nofreq, MAF):
		#create a look-up library such that each configuration only get computed once
		self.betabinomprobs = {}
		self.rawProbs = {}
		self.normBase = {}
		self.normProbs = {}

		for x in Grids.x :
			for a in Grids.abeta :
				#container for probs
				full_probs_list = np.zeros( InputData.numSites )
				#deal with different n separately
				for n in InputData.sampSizes:
					indice = np.where( InputData.total == n )
					subset_counts = InputData.count[indice]
					#get raw probs & fold dist'n
					if nofreq:
						raw_probs = 0.5 * ( self.get_raw_probs(subset_counts , n, x, a, 'B1') + self.get_raw_probs(subset_counts, n, 1.-x, a, 'B1') )
						probs_list = raw_probs / self.get_B1_normBase(n, x, a, InputData.minCount)
					elif MAF:
						raw_probs = 0.5 * ( self.get_raw_probs(subset_counts , n, x, a, 'B2maf') + self.get_raw_probs(subset_counts, n, 1.-x, a, 'B2maf') )
						probs_list = raw_probs / self.get_B2maf_normBase(n, x, a, InputData.minCount)
					else:
						raw_probs = 0.5 * ( self.get_raw_probs(subset_counts , n, x, a, 'B2') + self.get_raw_probs(subset_counts, n, 1.-x, a, 'B2') )
						probs_list = raw_probs / self.get_B2_normBase(n, x, a, InputData.minCount)

					#save
					full_probs_list[indice] = probs_list

				#store them. They should be indexed identically to SNP positions
				#these probs are all conditional on the given n. Will be multiplied by the proportion of n's in total input before mixing.
				self.normProbs[(x, a)] = full_probs_list

	#return a full list of normalized probs (for selection component) for each locus
	def get(self, x, a):
		return ( self.normProbs[(x, a)] )

	#return betabinom probs (for B2 and B2maf only)
	def _get_betabinom_probs(self, k, n, x, a):
		b = self._get_b(x, a)
		if (n, x, a) not in self.betabinomprobs:
			self.betabinomprobs[(n,x,a)] = betabinom(n, a, b)	

		return ( self.betabinomprobs[(n,x,a)].pmf(k) ) 

	#compatible with when k is np.array
	#do not fold the distn yet.
	def get_raw_probs(self, k, n, x, a, Stat):
		if (n, x, a) not in self.rawProbs:
			#for B1, polymorphism is 1, substitution is 0
			if Stat == "B1":
				#sanity check
				assert set(k) == set([0,1])
				b = self._get_b(x, a)
				probs = np.where( k == 0 , betabinom(n,a,b).pmf(n) , (1. - betabinom(n,a,b).pmf(n) - betabinom(n,a,b).pmf(n)) )
				self.rawProbs[(n, x, a)] = probs

			elif Stat == "B2":
				self.rawProbs[(n, x, a)] = self._get_betabinom_probs(k, n, x, a)

			elif Stat == "B2maf":
				probs = self._get_betabinom_probs(k, n, x, a) + self._get_betabinom_probs(n-k, n, x, a)
				#n/2 double counted if n is even
				if n % 2 == 0:
					probs = np.where( k == int(n/2), probs/2 , probs)
				
				self.rawProbs[(n, x, a)] = probs

		return ( self.rawProbs[(n, x, a)] )

	#folded
	def get_B2_normBase(self, n, x, a, minCount):
		if (n, x, a) not in self.normBase:
			excluded_probs = 0.5 * ( self._get_betabinom_probs( np.arange(minCount) , n, x, a) + self._get_betabinom_probs( np.arange(minCount) , n, 1. - x, a) )
			self.normBase[(n,x,a)] = 1. - np.sum(excluded_probs)
		
		return ( self.normBase[(n,x,a)] )

	def get_B2maf_normBase(self, n, x, a, minCount):
		if (n, x, a) not in self.normBase:
			#make sure the excluded counts do not include n bc it's substitution (and shows up as 0 in data)
			excluded_counts = np.concatenate( (np.arange(minCount), np.arange(n - minCount + 1, n) ) )
			excluded_probs = 0.5 * ( self._get_betabinom_probs( excluded_counts , n, x, a) + self._get_betabinom_probs( excluded_counts , n, 1. - x, a) )
			self.normBase[(n, x, a)] = 1. - np.sum(excluded_probs)

		return ( self.normBase[(n,x,a)] )

	def get_B1_normBase(self, n, x, a, minCount):
		#this is identical to B2's
		return ( self.get_B2_normBase(n, x, a, minCount) )

#import time

def calcBaller(window_indice, testSite, InputData, NeutralSFS, NormalizedBetaBinom, Grids):
	#define the window
	#window = np.array([InputData.genPos[i] for i in window_indice])
	window = InputData.genPos[ window_indice ]
	#dist = np.abs(window - testSite)
	dist = np.abs( InputData.genPos - testSite )

	#initiate optimization
	#CLRs will not be saved for now. Can be made to be written out in the future, be there interests in likelihood surfaces.
	## [CLR, x, abeta, A, numSites]
	Tmax = [0., 0., 0., 0., 0.]
	#last = time.time()
	for A in set(Grids.A):
		alphas = np.exp(- A*dist)
		effective_indice = np.where( (alphas >= 1e-8) & (InputData.genPos != testSite) )[0]
		# get intersect with  window_indice
		subwindow_indice = list( set(window_indice) & set(effective_indice) )
		if len(subwindow_indice) == 0:
			continue
		subset_alphas = alphas[subwindow_indice]
		
		#extract proportion of samples with size n
		prop_sampSize = NeutralSFS.propSizes[ subwindow_indice ]

		#loop through other Grids
		for x in set(Grids.x):
			for abeta in set(Grids.abeta):
				#get per-site neutral probs
				Neut_probs = NeutralSFS.probs[ subwindow_indice ]

				#extract pre-computed selection probs
				Sel_probs = NormalizedBetaBinom.get(x, abeta)[subwindow_indice]
				#sanity check
				try:
					assert len(Neut_probs) == len(Sel_probs)
					assert len(subwindow_indice) == len(prop_sampSize)
					assert len(Sel_probs) == len(prop_sampSize)
				except:
					print(f'len(Neut_probs) = {len(Neut_probs)}, len(Sel_probs) = {len(Sel_probs)}, len(subwindow_indice) = {len(subwindow_indice)}, len(prop_sampSize) = {len(prop_sampSize)}.\nA = {A}, x = {x}, abeta = {abeta}.')
					print(subwindow_indice, list(subwindow_indice), subwindow_indice[0], A, x, abeta)
					sys.exit()
				#account for conditional probs
				Sel_probs = Sel_probs * prop_sampSize
				#mixture
				mix_probs = subset_alphas * Sel_probs + (1. - subset_alphas) * Neut_probs
				#obtain the CLs
				CL_sel = np.sum( np.log(mix_probs) )
				CL_neut = np.sum( NeutralSFS.logProbs[subwindow_indice] )
				#get CLR
				T = 2 * (CL_sel - CL_neut)
				#pick
				if T > Tmax[0]:
					Tmax = [T, x, abeta, A, len(subwindow_indice)]
			#done for abeta loop
		#done for x loop
	#done for A loop
	#print(f'Looping through 3 grids in {time.time()-last}s.')
	return ( Tmax )


class Scan:

	'''Perform the scan with fixed window size of w (bp), with step size s (bp).'''
	def _fixSize_noCenter(self, InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile , w, s):
		print(("writing output to %s" % (outfile)))
		with open(outfile,'w') as scores:
			scores.write('physPos\tgenPos\tCLR\tx_hat\ts_hat\tA_hat\tnSites\n')#\tnumSites
			start = int( floor(2*float(position[0])/w)*(w/2) )
			end = start + s ; midpos = start + s/2
			start_i=0; end_i=0; pos_i=0
			while midpos <= InputData.position[-1]:
				#define the window
				while InputData.position[start_i] < start:
					start_i +=1 
				while InputData.position[pos_i] < midpos:
					pos_i += 1
				while end_i < numSites:
					if InputData.position[end_i] < end:
						end_i += 1
					else:
						break
				gen_site = midpos * InputData.Rrate
				#in case of large genomics gap
				if start_i >= end_i:
					scores.write('%g\t%s\t0\tNA\tNA\t0\n' % (midpos, gen_site))#
					start+=s; midpos+=s; end+=s
				else:
					window_indice = np.arange(start_i, end_i+1)
					Tmax, xhat, aHat, Ahat, winSites = calcBaller(window_indice, gen_site, InputData, NeutralSFS, NormalizedBetaBinom, Grids)
					scores.write( f'{midpos}\t{midpos * Rrate}\t{Tmax}\t{xhat}\t{aHat}\t{Ahat}\t{winSites}\n' )
					#scores.write('%g\t%s\t%s\t%s\t%g\t%g\t%d\n' % (midpos, midpos*Rrate, Tmax, xhat, aHat, Ahat, winSites))#
					#read in, take next step
					start+=s; midpos+=s; end+=s 
		scores.close()
		return(0)


	'''Perform the scan with fixed window size of w (bp), centered on every s informative sites'''
	def _fixSize_siteCenter(self, InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile , w, s):
		print(("writing output to %s" % (outfile)))
		with open(outfile,'w') as scores:
			scores.write('physPos\tgenPos\tCLR\tx_hat\ts_hat\tA_hat\tnSites\n')#\tnumSites
			i=0; 
			start_i=0; end_i=0
			while i < InputData.numSites:
				testSite = InputData.position[i]
				gen_site = InputData.genPos[i]
				start =  max(0, testSite-r/2); end = min(testSite+r/2,position[-1])
				while InputData.position[start_i] < start:
					start_i += 1
				while end_i < numSites:
					if InputData.position[end_i] < end:
						end_i += 1
					else:
						break
				try:
					assert end_i >= start_i
				except:
					print(start, start_i, end, end_i)
					sys.exit(1)
				end_i = min(end_i, numSites-1)
				window_indice = np.arange(start_i, end_i+1)
				Tmax, xhat, aHat, Ahat, winSites = calcBaller(window_indice, gen_site, InputData, NeutralSFS, NormalizedBetaBinom, Grids)
				scores.write( f'{testSite}\t{gen_site}\t{Tmax}\t{xhat}\t{aHat}\t{Ahat}\t{winSites}\n' )
				i+=int(s)
		scores.close()
		return(0)

	'''Perform the scan with windows centered on every s sites, with a radius of r sites on either side'''
	def _siteBased(self, InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile , r, s):
		print(("writing output to %s" % (outfile)))
		with open(outfile,'w') as scores:
			scores.write('physPos\tgenPos\tLR\tx_hat\ts_hat\tA_hat\tnSites\n')#\tnumSites
			i=0
			while i < InputData.numSites:
				testSite = InputData.genPos[int(i)]
				phys_site = InputData.position[int(i)]
				start_i = max(0, i-r) ; end_i = min(numSites,i+r+1)
				window_indice = np.arange(start_i, end_i+1)
				Tmax, xhat, aHat, Ahat, winSites = calcBaller(window_indice, gen_site, InputData, NeutralSFS, NormalizedBetaBinom, Grids)
				scores.write( f'{phys_site}\t{testSite}\t{Tmax}\t{xhat}\t{aHat}\t{Ahat}\t{winSites}\n' )
				i+=s
		scores.close()
		return(0)


	'''Perform the scan with windows centered on every s informative sites, taking peripheral SNPs with alpha>=1e-8, which is already a built-in feature for calcBaller'''
	def _alpha(self, InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile , s):
		print(("writing output to %s" % (outfile)))
		with open(outfile,'w') as scores:
			scores.write('physPos\tgenPos\tCLR\tx_hat\ts_hat\tA_hat\tnSites\n')#\tnumSites\tLa\tL0
			i=0
			while i < InputData.numSites:
				testSite = InputData.genPos[int(i)]
				phys_site = InputData.position[int(i)]
				Tmax, xhat, aHat, Ahat, winSites = calcBaller(np.arange(InputData.numSites), testSite, InputData, NeutralSFS, NormalizedBetaBinom, Grids)
				scores.write( f'{phys_site}\t{testSite}\t{Tmax}\t{xhat}\t{aHat}\t{Ahat}\t{winSites}\n' )
				i+=int(s)
		scores.close()
		return(0)


	def __init__(self, InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile, fixSize=False,r=0,s=1,phys=False,noCenter=False):
		if fixSize:
			print('You\'ve chosen to fix the size (in nt) of sliding window for scanning.')
			if r == 0:
				print('Please set a window width in nt with \"-w\" or \"--window\" command.')
				sys.exit()
			if not phys:
				print(f'Please make sure to use physical positions as coordinates if fixed-length windows are chosen (--fixSize). Scan will continue with physical positions with a rec rate of {InputeData.Rrate} cM/nt.')
				phys = True
				#sys.exit()
			if noCenter:
				w = float(r) 
				print ( ('Computing LR on %.3f kb windows on every %s nt. Using physical positions by default.' % (w/1e3, s)) )
				self._fixSize_noCenter(InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile , w, s)
			else: #site-centered
				w = float(r) 
				print(('Computing LR on %.3f kb windows on every %g informative sites. Using physical positions by default.' % (w/1e3,s)))
				self._fixSize_siteCenter(InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile , w, s)
		# When fixSize == False, and radius (-r) provided. Scan with fixed number of sites
		elif r != 0: 
			print(('Computing LR on every %s site/s, with a radius of %g informative sites on either side.' % (s, r)))
			self._siteBased(InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile , r, s)
		#window size not given 
		#then use all data (but the test site)
		else:
			print(('Computing LR on every %s site/s, using informative sites with exp(-A*dist) >= 1e-8.' % (s)))
			self._alpha(InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile , s)
		print(f'{datetime.now()}. Scan finished.')



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
def getSpect(infile,spectfile,MAF=False):
	Spect={}; numSites=0; translate = False
	with open(infile,'r') as sites:
		l=next(sites)
		for l in sites:
			(x,n)=[ int(i) for i in l.strip().split('\t')[2:] ]
			if MAF:
				try:
					assert x <= n/2
				except:
					if not translate:
						print(f'Input data includes non-MAF site/s (frequency >= 0.5) despite choosing to use B_2maf (with --MAF). These frequencies will be folded for following analyses.')
						translate = True
					x = n - x
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

'''
usage='python {} -i <input file> -o <output file> --spect <spect/config file> [--help] [--nofreq] [--nosub] [--MAF] [--getSpect] [--getConfig] [--fixSize] [--physPos] [--rec <recomb rate>] [-w <window size>] [--noCenter] [-s <step size>] [--fixX <x>] [--fixAlpha] [--findBal] [--findPos] [--rangeA <min,max,step>] [--listA <A1,A2,..,Ak>]'.format(sys.argv[0])
'''
def main():
	import argparse
	#parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input', dest='infile', help = 'Path and name of your input file.\n', required=True)
	parser.add_argument('-o','--output', dest='outfile', help = 'Path and name of your output file.\n')
	parser.add_argument('--spect',dest='spectfile', help = 'Path and name of the allele frequency spectrum file or configuration file.\n', required=True)
	parser.add_argument('--minCount', dest='minCount', default = 1, help = 'If rare variants are removed from the input, please provide the smallest allele count included in the input. Default value is 1.')

	parser.add_argument('--getSpect', dest='getSpec', action='store_true', default=False, help='Option to generate frequency spectrum file from the concatenated input file. Use \"-i\" and \"--spect\" commands to provide names and paths to input and output files, respectively. Indicate the input type with \"--MAF\".\n')
	parser.add_argument('--getConfig', dest='getConfig', action='store_true', default=False, help='Option to generate configuration file from the concatenated input file. Use \"-i\" and \"--spect\" commands to provide names and paths to input and output files, respectively.\n\n')

	parser.add_argument('--noFreq', dest='nofreq', action='store_true', default=False, help = 'Option to compute B_1 statistic and ignore allele frequency information (if given). All polymorphic sites (non-zero counts) will be considered as equivalent. Substitutions should be represented as having zero count in the input.')
	#deprecated --noSub bc B0 do not perform well with betabinomial
	#parser.add_argument('--noSub', dest='nosub', action='store_true', default=False, help = 'Option to not include substitution in input data.')
	parser.add_argument('--MAF', dest='MAF', action='store_true', default=False, help = 'Option to compute B_2maf statistic and use minor allele frequency instead of polarized allele frequency (if given). The latter is default (B_2 statisitc).')
	parser.add_argument('--findBal',dest='bal', action='store_true', default=False, help="Option to only look for footprints of balancing selection.\n")
	parser.add_argument('--findPos',dest='pos', action='store_true', default=False, help="Option to only look for footprints of positive selection.\n")

	#args concerning recombination rate / genetic map
	parser.add_argument('--usePhysPos', action='store_true', dest = 'phys', default = False, help = 'Option to use physical positions instead of genetic positions (in cM). Default is using genetic positions.\n')
	parser.add_argument('--rec', dest='Rrate', default = 1e-6 , help='The uniform recombination rate in cM/nt. Default value is 1e-6 cM/nt. Only useful when choose to use physical positions as coordinates.\n\n')

	#args to customize sliding windows
	parser.add_argument('--fixWinSize', action='store_true', dest = 'size', default = False, help = 'Option to fix the size (in nt) of sliding windows during scan. When true, please also provide the length of window in neucleotide (nt) with \"-w\" or \"--window\" command.\n')
	parser.add_argument('-w','--window', dest='w', type = int, default=0, help='Number of sites flanking the test locus on either side. When choose to fix window size (\"--fixSize\"), input the length of window in bp.\n')
	parser.add_argument('--noCenter', action='store_true', dest='noCenter', default=False, help = 'Option to have the scanning windows not centered on informative sites. Require that the window size (\"-w\") in physical positions (\"--physPos\") is provided. Default is True.\n')
	parser.add_argument('-s','--step', dest='step', type = float, default=1, help='Step size in bp (when using \"--noCenter\") or the number of informative sites. Default value is one site or one nucleotide.\n\n')

	#args to modify search grids
	parser.add_argument('--fixX', dest='x', help='Option to fix the presumed equilibrium frequency.\n')
	parser.add_argument('--fixAlpha', dest='abeta', type=float, default=None, help='Option to fix the alpha parameter in the beta-binomial distribution.\n')
	parser.add_argument('--rangeA', dest='seqA', help='Range of the values of the linkage parameter A to optimize over. Format should follow <Amin>,<Amax>,<Astep> with no space around commas.\n')
	parser.add_argument('--listA', dest='listA', help='Manually provide a list of A values to optimize over. Please separate the values with comma, no space.\n')

	
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	opt = parser.parse_args(sys.argv[1:])

	#generate helper files
	if opt.getSpec:
		print('You\'ve chosen to generate site frequency spectrum...')
		print(('Concatenated input: %s \nSpectrum file: %s' % (opt.infile, opt.spectfile )))
		getSpect(opt.infile, opt.spectfile, opt.MAF)
		sys.exit()
	elif opt.getConfig:
		print('You\'ve chosen to generate the substitution-polymorphism configuration...')
		print(('Concatenated input: %s \nConfiguration file: %s' % (opt.infile, opt.spectfile )))
		getConfig(opt.infile, opt.spectfile)
		sys.exit()

	#read & save input data
	print( f"\n{datetime.now()}. Reading input from {opt.infile}" )
	
	data = InputData(opt.infile, opt.nofreq, opt.MAF, opt.minCount, phys=opt.phys, Rrate=opt.Rrate)

	##read helper files & initialize
	Neutral =  NeutralSFS(opt.spectfile, opt.nofreq, opt.MAF)

	print( f'\n{datetime.now()}. Initializing...')
	print('Retrieving per-site neutral probabilities...')
	Neutral.get_neut_probs(data, opt.nofreq, opt.MAF)

	#generate the grids to optimize over/with
	grid = Grids(opt.x, opt.abeta, opt.bal, opt.pos, opt.seqA, opt.listA)
	print('\nOptimizing over x= '+', '.join([str(x) for x in grid.x]) )
	print('\n \t alpha= '+', '.join([str(a) for a in grid.abeta]) ) 
	print('\n \t A= '+', '.join([str(A) for A in grid.A]) ) 

	#pre-compute per-site probs across the grids for mixing later
	Sel_Probs = NormalizedBetaBinom(data, grid, opt.nofreq, opt.MAF)

	#start writing output
	print(("\n%s. Start computing likelihood raito..." %(datetime.now())))
	#size=False,r=0,s=1,phys=False,noCenter=False
	#InputData, NeutralSFS, NormalizedBetaBinom, Grids, outfile, fixSize=False,r=0,s=1,phys=False,noCenter=False
	Scan(data, Neutral, Sel_Probs, grid, opt.outfile, fixSize = opt.size, r = opt.w, s = opt.step, phys = opt.phys, noCenter = opt.noCenter)

	#pipeline finished
	print(f'\n{datetime.now()}. Pipeline finished.')

if __name__ == '__main__':
	main()
