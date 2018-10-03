import sys
import math
import numpy as np

class Pheno:

    def __init__(self,size = 2550858,numbins = 100,pvThresh=1.,infoFilt=1.,rsPrune=False,useAbs=False):
    
        self.useAbs = useAbs
        self.size  = np.zeros(size)
        self.effects =  np.zeros(size)
        self.infoFilt = infoFilt
        self.rsPrune = rsPrune
        self.rsList = {}
        self.pvThresh = pvThresh
        self.removeAmbig = True
        self.ancest = {}
        self.freqs =  np.zeros(size)
        self.rsids =  np.array(np.zeros(size),dtype=object)
        self.rsidsLook = {}
        self.a1s =  np.array(np.zeros(size),dtype=object)
        self.a2s =  np.array(np.zeros(size),dtype=object)
        self.p  =  np.zeros(size)
        self.sum_stat = 0.
        self.numbins =numbins
        self.all_freqs = [[[0.,0.] for i in range(self.numbins)] for j in range(1703)]
        self.der_betas = []
        self.mean_beta = 0.
        self.null = []
        self.permBlocks = []
        self.pop_freq_prune = {}

    def getRs(self,f):
        rsf = open(f,'r')
        for line in rsf:
            self.rsList[line.strip().split()[0]] =0
        rsf.close()

    def getAncest(self,f):
        ah = open(f, 'r')
        for line in ah:
            data = line.strip().split()
            aa = data[5].split('=')[1][0]
            self.ancest[data[2]] = aa
        ah.close()

    def getData(self,f,t1,t2,pf=False): 
        fh = open(f, 'r')
 
        self.numbins = int(round(self.numbins*(t2-t1)))

        self.all_freqs = [[[0.,0.] for i in range(self.numbins)] for j in range(1703)]

        for line in fh:
  
            der_all_beta = 0
            der_all_freq = 0

            data = line.strip().split()
            
            # filter out lines with trailing data or pvalues exceeding the threshold
            if len(data) > 10 and float(data[10]) > self.pvThresh:
                continue

            # FILTER out alleles with low imputation quality
            if float(data[8]) < self.infoFilt:
                continue

            # if desired, prune to a list of alleles with population frequencies from TGP
            if pf:
                if data[0] not in self.pop_freq_prune:
                    continue

            # if desired, prune to a precomputed list of rsnumbers
            if self.rsPrune:
                if data[0] not in self.rsList:
                    continue
           
            # if desired, remove alleles that are strand ambigious, which could be mismatched for the ancestral allele calls
            if self.removeAmbig and ((data[3] == "A" and data[4] == "T") or (data[3] == "T" and data[4] == "A") or (data[3] == "G" and data[4] == "C") or (data[3] == "C" and data[4] == "G")):
                continue

            # some alleles are ecoded lower case, make upper for consistency
            data[5] = data[5].upper()
            data[3] = data[3].upper()
           
            # remove alleles that do not have ancestral calls
            if data[5] != "A" and data[5] != "T" and data[5] != "G" and data[5] != "C":
               continue

            # assign the derived allele beta
            if data[5] != data[3]:
                der_all_freq = float(data[6])
                der_all_beta = float(data[7])
            else:
                der_all_freq = 1.-float(data[6])
                der_all_beta = -1.*float(data[7])

            # filter out allele if it does not fall in desired range
            if der_all_freq <= t1 or der_all_freq >= t2:
                continue

            # filter out fixed alleles
            if der_all_freq == 1:
                continue

            if der_all_freq == 0:
                continue
     
            # set the frequency index 
            index = int(math.floor((der_all_freq-t1)*self.numbins*1./(t2-t1))) 

            self.all_freqs[int(data[9])][index][0] += der_all_beta
            self.all_freqs[int(data[9])][index][1] += 1

            self.der_betas.append(der_all_beta)

        self.mean_beta = np.nanmean(self.der_betas)
    
    def get_S_beta(self,arr):
        sum_stat = 0
        for i in xrange(len(arr)):
            if arr[i] == 'nan':
                continue
            if self.useAbs:
                sum_stat += abs(arr[i])
                continue
            sum_stat += arr[i]
        self.sum_stat = sum_stat
    
    def get_S_beta_arr(self,arr):
        retarr = np.zeros(len(arr))
        retarr[0] = arr[0]
        if self.useAbs:
            sum_stat += abs(arr[0])
        for i in xrange(len(arr)):
            if arr[i] == 'nan':
                continue
            if self.useAbs:
                retarr[i] = retarr[i-1]+abs(arr[i])
                continue
            retarr[i] = retarr[i-1]+arr[i]
        return retarr
 
    def get_pop_diff(self, thresh):
     
        fh = open('/Users/uricchio/projects/data/1KG/pop_diff.GBR.TSI.txt', 'r')

        for line in fh:
            data = line.strip().split()
            data[1] = float(data[1])       
            data[2] = float(data[2])       
            if abs(data[1]-data[2]) <  thresh:
                self.pop_freq_prune[data[0]] = True
        fh.close()

    def runIter(self,perm=False,pr=False):
        betas = [[0,0] for i in range(self.numbins)]
        flip = []
        if not perm:
            flip = [1 for i in range(1703)] 
        if perm:
            flip = np.random.choice([-1,1], size=1703) 
        j = 0
        for chunk in self.all_freqs:
            i = 0
            for thing in chunk:
                betas[i][0] += flip[j]*thing[0]
                betas[i][1] += thing[1] 
                i += 1
            j += 1
        j = 0
        
        allbetas = []
        sfs = []
        meanbeta = 0
        tot = 0
        for thing in betas:
            if thing[1] > 0:
                allbetas.append(thing[0]/thing[1]) #-self.mean_beta)
            else:
                allbetas.append('nan')
            sfs.append(thing[1])
            meanbeta += thing[0]
            tot += thing[1]
        self.get_S_beta(allbetas)
        if pr:
            k = 0
            for al in allbetas:
                if not self.useAbs:
                    print al, betas[k][1]
                else:
                    if al == 'nan':
                        print al, betas[k][1]
                    else:
                        print  abs(al), betas[k][1]
                k += 1
        if not perm:
            print "#", 
        print self.sum_stat 
        sys.stdout.flush()
        
    def getSize(self,f):
        numbins = self.numbins
        infoFilt = self.infoFilt
        pvThresh = self.pvThresh
        rsPrune = self.rsPrune
        useAbs = self.useAbs
  
        #print >> sys.stderr, numbins
        with open(f) as fh:
            for i, l in enumerate(fh):
                pass
        size = i + 1
        self.__init__(size=size,numbins=numbins,infoFilt=infoFilt,pvThresh=pvThresh,rsPrune=rsPrune,useAbs=useAbs)
