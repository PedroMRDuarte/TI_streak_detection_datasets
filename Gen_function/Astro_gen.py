import numpy as np
import skimage.io as io
import os
from skimage.morphology import disk,dilation
from skimage.draw import line



class Astro_streak_Image():
    def __init__(self,n,m):
        noise = 1.5
        self.n = n
        self.m = m

    def build8bit(self,I):
        I8bit = I
        #print(np.mean(I8bit))
        I8bit = (I8bit - I8bit.min()) / (I8bit.max() - I8bit.min()) #normalizes data in range 0 - 255
        I8bit = 255 * I8bit  
        return I8bit.astype('uint8')

    def ex3_4(self,ex1,ex2,x3,x4):
        m = (ex1[1]-ex2[1])/(ex1[0]-ex2[0])
        # b = y - m*x
        b = ex1[1] - m * ex1[0]
        ex3 = (x3, int(m*x3+b))
        ex4 = (x4, int(m*x4+b))
        return ex3,ex4

    def generate_image(self,ex1,ex2,streakint,Istars):
        # n,m size of image
        # ex1 -> extreme 1 as (x1,y1)
        # ex2 -> extreme 2 as (x2,y2)
        # thick -> int
        # seed_stars -> seed of the random generated stars.
        
        I = np.zeros((self.n,self.m))
        np.random.seed(seed=None)
        noise = 1.5
        noiseI = np.random.normal(loc = 11.44, scale = noise,size = I.shape)
        
        rr,cc = line(ex1[0],ex1[1],ex2[0],ex2[1])
        I[rr, cc] = 1
        #Ipred = I.copy()
        d = disk(2)
        I = dilation(I,d)
        SNRpower = streakint*np.mean(noiseI)
        I = I*SNRpower
        #Ipred = (Ipred == 1).astype(np.int)
        #print(I.max())
        bstreak = (I > 0)
        bstars = (Istars>0)
        bnoise = 1 - (bstreak+bstars)
        I = I + np.random.normal(loc = 0, scale = 3,size = I.shape)
    
        I = I*bstreak + Istars + noiseI*bnoise

        #add noise
        I8bit = I
        return I8bit



    def dil(self,d,I):
        di = disk(d)
        I = dilation(I,di) 
        return I

    def stars(self,seed_stars,numstrs,rn):
        n = self.n
        m = self.m
        Istars = np.zeros((n,m))
        np.random.seed(seed_stars)
        a = int(numstrs//5)
        ab = numstrs-int(4*numstrs//5)-rn
        strs1 = np.random.rand(a,2)
        strs2 = np.random.rand(a,2)
        strs3 = np.random.rand(a,2)
        strs4 = np.random.rand(a,2)
        strs5 = np.random.rand(ab,2)
        strs6 = np.random.rand(rn,2)

        l = 1
        Istarsaux1 = np.zeros((n,m))
        Istarsaux2 = np.zeros((n,m))
        Istarsaux3 = np.zeros((n,m))
        Istarsaux4 = np.zeros((n,m))
        Istarsaux5 = np.zeros((n,m))
        Istarsaux6 = np.zeros((n,m))
        vals = np.random.normal(loc = 11, scale = 1.5,size = (numstrs,1))*255 
        for i in range(a):
            
            Istarsaux1[int(strs1[i,0]*n),int(strs1[i,1]*m)] = vals[i]
            Istarsaux2[int(strs2[i,0]*n),int(strs2[i,1]*m)] = vals[i+a]
            Istarsaux3[int(strs3[i,0]*n),int(strs3[i,1]*m)] = vals[i+2*a]
            Istarsaux4[int(strs4[i,0]*n),int(strs4[i,1]*m)] = vals[i+3*a]
            if i < (ab):
                Istarsaux5[int(strs5[i,0]*n),int(strs5[i,1]*m)] = vals[i+4*a]
            if i < rn:
                if i == 1:
                    Istarsaux6[int(strs6[i,0]*n),int(strs6[i,1]*m)] = np.random.rand()*(0.19+0.8)*255
                else:
                    Istarsaux6[int(strs6[i,0]*n),int(strs6[i,1]*m)] = np.random.rand()*(0.49+0.5)*255
            #d = disk(np.floor(6/l))
            #Istarsaux = dilation(Istarsaux,d)
            # l = l+1
            # if l == 7:
            #     l=1
        Istarsaux2 = self.dil(2,Istarsaux2)
        Istarsaux3 = self.dil(4,Istarsaux3)
        Istarsaux4 = self.dil(4,Istarsaux4)
        Istarsaux5 = self.dil(6,Istarsaux5)
        Istarsaux6 = self.dil(10,Istarsaux6)

        Istars = Istars + Istarsaux1 + Istarsaux2 + Istarsaux3 + Istarsaux4 + Istarsaux5 + Istarsaux6
        return Istars.astype('uint8')

    def save(self,ex1,ex2,I8bit,SNR):
        #import matplotlib.image
        #from matplotlib import cm
        #matplotlib.image.imsave('Dataset\IN_%s_%s_%s_%s_%.2f.png' %(ex1[1],ex1[0],ex2[1],ex2[0],SNR),(I8bit).astype('uint8'),vmin= 0, vmax=255)#.astype(np.uint8),cmap = cm.gray)
        io.imsave('Dataset\IN_%s_%s_%s_%s_%.2f.png' %(ex1[1],ex1[0],ex2[1],ex2[0],SNR),(I8bit).astype('uint8'),check_contrast=False)

    
def main(n,m,ex1,ex2,seed = int(np.random.randint(low = 0, high =1e8)),numstars = 80, rn = 5, SNRlist = [5,4.5,4,3.5,3,2.5,2,1.5,1.25]):
    #call class
    img = Astro_streak_Image(n,m)
    #generate the stars mask
    Istars = img.stars(seed, numstars, rn)
    xlen = abs(ex1[0]-ex2[0])
    if os.path.isdir('Dataset') == False:
        os.makedirs('Dataset')

    a = int(0.1*xlen)
    if max(ex1[0],ex2[0])>=(n/2):
        x3 = int(min(ex1[0],ex2[0])-a)
        x4 = int(x3-xlen)
        x5 = int(x4-a)
        x6 = int(x5-xlen)
        x7 = int(x6-a)
        x8 = int(x7-xlen)
    else:
        x3 = int(max(ex1[0],ex2[0])+a)
        x4 = int(x3+xlen)
        x5 = int(x4+a)
        x6 = int(x5+xlen)
        x7 = int(x6+a)
        x8 = int(x7+xlen)

    ex3,ex4 = img.ex3_4(ex1,ex2,x3,x4)
    ex5,ex6 = img.ex3_4(ex1,ex2,x5,x6)
    ex7,ex8 = img.ex3_4(ex1,ex2,x7,x8)

    for i in SNRlist:
        SNR = i
        I8bit = img.generate_image(ex1,ex2,i,Istars)
        img.save(ex1,ex2,I8bit,SNR)

        I8bit = img.generate_image(ex3,ex4,i,Istars)
        img.save(ex3,ex4,I8bit,SNR)
        
        I8bit = img.generate_image(ex5,ex6,i,Istars)
        img.save(ex5,ex6,I8bit,SNR)
        I8bit = img.generate_image(ex7,ex8,i,Istars)
        img.save(ex7,ex8,I8bit,SNR)


ex1 = [1800,2700]
ex2 = [1825,2710]

## Change the input values:
# main function automatically generates 4 images with colinear streaks.
    #  n, m, ex1, ex2, seed = %int , numstars = 80(default), rn = 5(default), SNRlist = [5,4.5,4,3.5,3,2.5,2,1.5,1.25](default)
    # n --> x-length of the generated image
    # m --> y-length of the generated image
    # ex1 --> list with the 1st extreme coordinates of the first streak
    # ex2 --> list with the 2nd extreme coordinates of the first streak
    # seed --> seed of the stars
    # numstars --> number of stars in the image
    # rn --> number of bright stars
    # SNRlist --> list of the SNR values of the streak to generate
    

main(2000,3000, ex1, ex2,SNRlist=[7,5])




