# Saves numpy arrays for importing SPH data into SKIRT

import pynbody
import numpy as np
import math
from timeit import default_timer as timer
from os.path import expanduser
import os
import argparse
from scipy import spatial
from scipy.spatial.transform import Rotation as R

class galaxy:
    def __init__(self, z, name):
        if z == '2.0':
            filePath = '/scratch/ntf229/nihao2/n10/z'+z+'/'+name+'/'+name+'.00240'
        elif z == '3.6':
            filePath = '/scratch/ntf229/nihao2/n10/z'+z+'/'+name+'/'+name+'.00128'
        else:
            print('invalid redshift')
            exit()
        
        # Load NIHAO data
        self.data = pynbody.load(filePath)
        
        self.full_x_pos = np.float32(self.data.star['pos'].in_units('pc')[:][:,0])
        self.full_y_pos = np.float32(self.data.star['pos'].in_units('pc')[:][:,1])
        self.full_z_pos = np.float32(self.data.star['pos'].in_units('pc')[:][:,2])
        self.full_mass = np.float32(self.data.star['massform'].in_units('Msol')[:]) # in solar masses	
        self.full_current_mass = np.float32(self.data.star['mass'].in_units('Msol')[:])
        self.full_metals = np.float32(self.data.star['metals'][:])
        self.full_age = np.float32(self.data.star['age'].in_units('yr')[:])
        self.full_x_vel = np.float32(self.data.star['vel'].in_units('km s**-1')[:][:,0])
        self.full_y_vel = np.float32(self.data.star['vel'].in_units('km s**-1')[:][:,1])
        self.full_z_vel = np.float32(self.data.star['vel'].in_units('km s**-1')[:][:,2])
        self.full_smooth = 2*np.float32(self.data.star['smooth'].in_units('pc')[:]) # 2 times gravitational softening length
        
        self.full_x_pos_dust = np.float32(self.data.gas['pos'].in_units('pc')[:][:,0])
        self.full_y_pos_dust = np.float32(self.data.gas['pos'].in_units('pc')[:][:,1])
        self.full_z_pos_dust = np.float32(self.data.gas['pos'].in_units('pc')[:][:,2])
        self.full_smooth_dust = 2*np.float32(self.data.gas['smooth'].in_units('pc')[:]) # 2 times gravitational softening length
        self.full_mass_dust = np.float32(self.data.gas['mass'].in_units('Msol')[:]) # in solar masses	
        self.full_metals_dust = np.float32(self.data.gas['metals'][:])
        self.full_temp_dust = np.float32(self.data.gas['temp'][:])
        self.full_density_dust = np.float32(self.data.gas['rho'].in_units('Msol pc**-3')[:])

        print('gas density units:',self.data.gas['rho'].units)
        print('max stellar particle mass:',np.amax(self.full_mass))
        print('min stellar particle mass:',np.amin(self.full_mass))
        
        self.full_length_star = len(self.full_x_pos) # starting length
        print(self.full_length_star, 'full star') 
        
        self.full_length_dust = len(self.full_x_pos_dust) # starting length
        print(self.full_length_dust, 'full dust') 
                
        # Halo catalogue
        h = self.data.halos() # ordered by number of particles (starts at 1)
        haloNum = int(1)
        
        xMin = np.amin(h[haloNum]['x'].in_units('pc'))
        xMax = np.amax(h[haloNum]['x'].in_units('pc'))
        yMin = np.amin(h[haloNum]['y'].in_units('pc'))
        yMax = np.amax(h[haloNum]['y'].in_units('pc'))
        zMin = np.amin(h[haloNum]['z'].in_units('pc'))
        zMax = np.amax(h[haloNum]['z'].in_units('pc'))
        
        xLength = abs(xMax - xMin)
        yLength = abs(yMax - yMin)
        zLength = abs(zMax - zMin)
        
        diameter = np.amax([xLength,yLength,zLength]) / 6 # division by 6 here scales all galaxies (calibrated from images)  
        
        xCenter = (xMax + xMin)/2
        yCenter = (yMax + yMin)/2
        zCenter = (zMax + zMin)/2
        
        # Make cuts based on position
        self.xAbsMin = xCenter - (diameter/2)															
        self.xAbsMax = xCenter + (diameter/2)
        self.yAbsMin = yCenter - (diameter/2)
        self.yAbsMax = yCenter + (diameter/2)
        self.zAbsMin = zCenter - (diameter/2)
        self.zAbsMax = zCenter + (diameter/2)
        
    def starCut(self):
        xMin = self.full_x_pos > self.xAbsMin
        xMax = self.full_x_pos[xMin] < self.xAbsMax
        yMin = self.full_y_pos[xMin][xMax] > self.yAbsMin
        yMax = self.full_y_pos[xMin][xMax][yMin] < self.yAbsMax
        zMin = self.full_z_pos[xMin][xMax][yMin][yMax] > self.zAbsMin
        zMax = self.full_z_pos[xMin][xMax][yMin][yMax][zMin] < self.zAbsMax
        
        self.x_pos = self.full_x_pos[xMin][xMax][yMin][yMax][zMin][zMax]
        self.y_pos = self.full_y_pos[xMin][xMax][yMin][yMax][zMin][zMax]
        self.z_pos = self.full_z_pos[xMin][xMax][yMin][yMax][zMin][zMax]
        self.smooth = self.full_smooth[xMin][xMax][yMin][yMax][zMin][zMax]
        self.x_vel = self.full_x_vel[xMin][xMax][yMin][yMax][zMin][zMax]
        self.y_vel = self.full_y_vel[xMin][xMax][yMin][yMax][zMin][zMax]
        self.z_vel = self.full_z_vel[xMin][xMax][yMin][yMax][zMin][zMax]
        self.mass = self.full_mass[xMin][xMax][yMin][yMax][zMin][zMax]
        self.current_mass = self.full_current_mass[xMin][xMax][yMin][yMax][zMin][zMax]
        self.metals = self.full_metals[xMin][xMax][yMin][yMax][zMin][zMax]
        self.age = self.full_age[xMin][xMax][yMin][yMax][zMin][zMax]
        
        self.starLength = len(self.x_pos)
        print(self.starLength, 'stars')
        
    def dustCut(self):    
        xMin = self.full_x_pos_dust > self.xAbsMin
        xMax = self.full_x_pos_dust[xMin] < self.xAbsMax
        yMin = self.full_y_pos_dust[xMin][xMax] > self.yAbsMin
        yMax = self.full_y_pos_dust[xMin][xMax][yMin] < self.yAbsMax
        zMin = self.full_z_pos_dust[xMin][xMax][yMin][yMax] > self.zAbsMin
        zMax = self.full_z_pos_dust[xMin][xMax][yMin][yMax][zMin] < self.zAbsMax
        
        self.x_pos_dust = self.full_x_pos_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.y_pos_dust = self.full_y_pos_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.z_pos_dust = self.full_z_pos_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.smooth_dust = self.full_smooth_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.mass_dust = self.full_mass_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.metals_dust = self.full_metals_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.temp_dust = self.full_temp_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.density_dust = self.full_density_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        
        self.dustLength = len(self.x_pos_dust)
        print(self.dustLength, 'dust')
        
    def shift(self):    
        xCenter = (self.xAbsMax + self.xAbsMin)/2
        yCenter = (self.yAbsMax + self.yAbsMin)/2
        zCenter = (self.zAbsMax + self.zAbsMin)/2
        
        self.x_pos = self.x_pos - xCenter
        self.y_pos = self.y_pos - yCenter
        self.z_pos = self.z_pos - zCenter
        
        self.x_pos_dust = self.x_pos_dust - xCenter
        self.y_pos_dust = self.y_pos_dust - yCenter
        self.z_pos_dust = self.z_pos_dust - zCenter

    def ageSmooth(self):
        k = 100 # k-th nearest neighbor to calculate deltaT
        scale = 100 # scales width of lognormal distribution relative to deltaT 
        numSamples = 50 # number of samples per parent particle
        ageLim = np.amax(self.age)
        newAges = []
        newMasses = []
        new_x_pos = []
        new_y_pos = []
        new_z_pos = []
        new_smooth = []
        new_x_vel = []
        new_y_vel = []
        new_z_vel = []
        new_current_mass = []
        new_metals = []
        # sort star particles by age
        ind = np.argsort(self.age)
        self.x_pos = self.x_pos[ind]
        self.y_pos = self.y_pos[ind]
        self.z_pos = self.z_pos[ind]
        self.smooth = self.smooth[ind]
        self.x_vel = self.x_vel[ind]
        self.y_vel = self.y_vel[ind]
        self.z_vel = self.z_vel[ind]
        self.mass = self.mass[ind]
        self.current_mass = self.current_mass[ind]
        self.metals = self.metals[ind]
        self.age = self.age[ind]
        for i in range(len(self.age)):
            parentAge = self.age[i]
            parentMass = self.mass[i]
            if i < k: # no k-th younger neighbor
                deltaT = self.age[i+k]-parentAge
            elif i+k >= len(self.age): # no k-th older neighbor
                deltaT = parentAge-self.age[i-k]
            else: # can find k-th neighbor before and after
                deltaT = np.amax([self.age[i+k]-parentAge, parentAge-self.age[i-k]])
            sampledAges = np.random.normal(loc=parentAge, scale=deltaT*scale, size=numSamples)
            num = (np.asarray(sampledAges) > ageLim).sum() + (np.asarray(sampledAges) < 0).sum()
            while num > 0:
                newSampledAges = np.random.normal(loc=parentAge, scale=deltaT*scale, size=num)
                sampledAges = np.append(sampledAges[sampledAges <= ageLim][sampledAges[sampledAges <= ageLim] > 0], newSampledAges)
                num = (np.asarray(sampledAges) > ageLim).sum() + (np.asarray(sampledAges) < 0).sum()
            sampledMasses = np.zeros(numSamples)+(parentMass/numSamples)
            newAges.extend(sampledAges.tolist())
            newMasses.extend(sampledMasses.tolist())
            new_x_pos.extend((np.zeros(numSamples) + self.x_pos[i]).tolist())
            new_y_pos.extend((np.zeros(numSamples) + self.y_pos[i]).tolist())
            new_z_pos.extend((np.zeros(numSamples) + self.z_pos[i]).tolist())
            new_smooth.extend((np.zeros(numSamples) + self.smooth[i]).tolist())
            new_x_vel.extend((np.zeros(numSamples) + self.x_vel[i]).tolist())
            new_y_vel.extend((np.zeros(numSamples) + self.y_vel[i]).tolist())
            new_z_vel.extend((np.zeros(numSamples) + self.z_vel[i]).tolist())
            new_current_mass.extend((np.zeros(numSamples) + (self.current_mass[i]/numSamples)).tolist())
            new_metals.extend((np.zeros(numSamples) + self.metals[i]).tolist())
        # change from lists to numpy arrays
        self.age = np.asarray(newAges)
        self.mass = np.asarray(newMasses)
        self.x_pos = np.asarray(new_x_pos)
        self.y_pos = np.asarray(new_y_pos)
        self.z_pos = np.asarray(new_z_pos)
        self.smooth = np.asarray(new_smooth)
        self.x_vel = np.asarray(new_x_vel)
        self.y_vel = np.asarray(new_y_vel)
        self.z_vel = np.asarray(new_z_vel)
        self.current_mass = np.asarray(new_current_mass)
        self.metals = np.asarray(new_metals)
              
    def youngStars(self, tauClear):
        tau_clear = float(tauClear) * 1.e6 # convert from Myrs to years
        k = 1.3806485*10**(-19) # boltzmann constant in cm**2 kg s**-2 K**-1
        youngStarMask = self.age < 1.e7 # Mask for young stars (smooth_time + 10 Myrs)
        youngStarIndex = [] # Indices of young stars 
        for i in range(self.starLength):
            if youngStarMask[i]:
                youngStarIndex.append(i)
        # MAPPINGS-III arrays
        self.young_x_pos = []
        self.young_y_pos = []
        self.young_z_pos = []
        self.young_mass = []
        self.young_current_mass = []
        self.young_metals = []
        self.young_age = []
        self.young_x_vel = []
        self.young_y_vel = []
        self.young_z_vel = []
        self.young_smooth = []
        self.young_SFR = []
        self.young_logC = []
        self.young_p = []
        self.young_f_PDR = []
        def mass_prob(x): # star cluster mass probability distribution
            return x**(-1.8)
        mass_min = 700
        mass_max = 1e6
        prob_max = 7.57*10**(-6) # slightly larger than 700**(-1.8)
        N_MC = 10000000 # masked distribution is well resolved and still runs fast
        # accept / reject Monte Carlo:
        mass = np.random.uniform(mass_min,mass_max,N_MC)  # get uniform temporary mass values
        prob = np.random.uniform(0,prob_max,N_MC)  # get uniform random probability values
        mask = prob < mass_prob(mass) # accept / reject
        sampled_masses = mass[mask] # sample of star cluster masses following the desired distribution
        for i in range(len(youngStarIndex)):
            parent_index = youngStarIndex[i]
            parent_mass = self.mass[parent_index]
            parent_age = self.age[parent_index]
            ind = len(self.young_x_pos) # don't subtract 1, haven't appended yet
            # create MAPPINGS-III particle
            self.young_x_pos.append(self.x_pos[parent_index])
            self.young_y_pos.append(self.y_pos[parent_index])
            self.young_z_pos.append(self.z_pos[parent_index])
            self.young_current_mass.append(self.mass[parent_index]) # assume no mass loss 
            self.young_metals.append(self.metals[parent_index])
            self.young_age.append(self.age[parent_index])
            self.young_x_vel.append(self.x_vel[parent_index])
            self.young_y_vel.append(self.y_vel[parent_index])
            self.young_z_vel.append(self.z_vel[parent_index])
            self.young_smooth.append(self.smooth[parent_index])
            self.young_f_PDR.append(np.exp(-1*parent_age / tau_clear)) # Groves 2008 equation 16
            self.young_mass.append(parent_mass)
            self.young_SFR.append(1.e-7 * self.young_mass[ind]) # (units: yr**-1) assumes constant SFR over the last 10 Myrs
            self.young_logC.append(np.random.normal(5., 0.4))
            self.young_p.append(k*10**((5/2)*(self.young_logC[ind]-(3/5)*np.log10(np.random.choice(sampled_masses))))) # in kg cm**-1 s**-2
            self.young_p[ind] *= 100 # convert to Pascals
            # create ghost particle
            self.x_pos_dust = np.append(self.x_pos_dust, self.x_pos[parent_index])
            self.y_pos_dust = np.append(self.y_pos_dust, self.y_pos[parent_index])
            self.z_pos_dust = np.append(self.z_pos_dust, self.z_pos[parent_index])
            self.smooth_dust = np.append(self.smooth_dust, 3*self.smooth[parent_index]) # 3 times smoothing length of parent star
            self.mass_dust = np.append(self.mass_dust, -100.*self.young_mass[ind]*self.young_f_PDR[ind]) # 100 times negative mass in PDR (1% SF efficiency) 
            self.metals_dust = np.append(self.metals_dust, self.metals[parent_index]) 
            self.temp_dust = np.append(self.temp_dust, 7000.) # assume 7000K (doesn't make a difference as long as lower than maxTemp)
            self.density_dust = np.append(self.density_dust, 0.) # don't need density, set to 0
        # change MAPPINGS-III from lists to numpy arrays
        self.young_x_pos = np.asarray(self.young_x_pos)
        self.young_y_pos = np.asarray(self.young_y_pos)
        self.young_z_pos = np.asarray(self.young_z_pos)
        self.young_mass = np.asarray(self.young_mass)
        self.young_current_mass = np.asarray(self.young_current_mass)
        self.young_metals = np.asarray(self.young_metals)
        self.young_age = np.asarray(self.young_age)
        self.young_x_vel = np.asarray(self.young_x_vel)
        self.young_y_vel = np.asarray(self.young_y_vel)
        self.young_z_vel = np.asarray(self.young_z_vel)
        self.young_smooth = np.asarray(self.young_smooth)
        self.young_SFR = np.asarray(self.young_SFR)
        self.young_logC = np.asarray(self.young_logC)
        self.young_p = np.asarray(self.young_p)
        self.young_f_PDR = np.asarray(self.young_f_PDR)
        self.young_f_PDR0 = np.zeros(len(self.young_f_PDR)) # for unattenuated SEDs
        # delete parent particles
        self.x_pos = np.delete(self.x_pos, youngStarIndex)
        self.y_pos = np.delete(self.y_pos, youngStarIndex)
        self.z_pos = np.delete(self.z_pos, youngStarIndex)
        self.mass = np.delete(self.mass, youngStarIndex)
        self.current_mass = np.delete(self.current_mass, youngStarIndex)
        self.metals = np.delete(self.metals, youngStarIndex)
        self.age = np.delete(self.age, youngStarIndex)
        self.x_vel = np.delete(self.x_vel, youngStarIndex)
        self.y_vel = np.delete(self.y_vel, youngStarIndex)
        self.z_vel = np.delete(self.z_vel, youngStarIndex)
        self.smooth = np.delete(self.smooth, youngStarIndex)
        
if __name__=='__main__':
    start = timer()
    parser = argparse.ArgumentParser()
    parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (see galaxy.youngStars()) 
    parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs
    parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True)
    parser.add_argument("--z") # redshift
    parser.add_argument("--galaxy") # name of galaxy 
    args = parser.parse_args()
    # Directory structure stores important parameters
    particlePath = '/scratch/ntf229/nihao2/Particles/'
    if eval(args.ageSmooth):
        particlePath += 'ageSmooth/'
    else:
        particlePath += 'noAgeSmooth/'
    if eval(args.SF):
        particlePath += 'SF/tauClear'+args.tauClear+'/'
    else:
        particlePath += 'noSF/'
    particlePath += 'z'+args.z+'/'
    particlePath += args.galaxy+'/'
    os.system('mkdir -p '+particlePath)
    g = galaxy(args.z, args.galaxy)
    g.starCut()
    g.dustCut()
    g.shift()
    if eval(args.ageSmooth):
        g.ageSmooth()
    if eval(args.SF):
        g.youngStars(args.tauClear)
        np.save(particlePath+'youngStars.npy',np.float32(np.c_[g.young_x_pos, g.young_y_pos, g.young_z_pos, g.young_smooth, g.young_x_vel, g.young_y_vel, g.young_z_vel, g.young_SFR, g.young_metals, g.young_logC, g.young_p, g.young_f_PDR]))
        np.save(particlePath+'youngStars_f_PDR0.npy',np.float32(np.c_[g.young_x_pos, g.young_y_pos, g.young_z_pos, g.young_smooth, g.young_x_vel, g.young_y_vel, g.young_z_vel, g.young_SFR, g.young_metals, g.young_logC, g.young_p, g.young_f_PDR0]))
    np.save(particlePath+'gas.npy',np.float32(np.c_[g.x_pos_dust, g.y_pos_dust, g.z_pos_dust, g.smooth_dust, g.mass_dust, g.metals_dust, g.temp_dust]))
    np.save(particlePath+'stars.npy',np.float32(np.c_[g.x_pos, g.y_pos, g.z_pos, g.smooth, g.x_vel, g.y_vel, g.z_vel, g.mass, g.metals, g.age]))
    np.save(particlePath+'current_mass_stars.npy',np.float32(np.c_[g.current_mass]))
    np.save(particlePath+'gas_density.npy',np.float32(np.c_[g.density_dust]))
    end = timer()
    print('time: ', end - start)

