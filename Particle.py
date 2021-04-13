# Particle formation for the decoupled chemistry by Ali Naseri
import numpy as np
AV = 6.022E23 # Avogadro number
AMU = 1.0/AV # Atomic Mass Unit
BOLTZMANN = 1.3806488E-16 # Boltzman constant
C_MW = 12.011 # Molecular weight of carbon
C_MASS = C_MW * AMU # Mass of a single carbon atom
PAH_rad = 7.24 * 1E-8 /2 # PAH radius in cm
Ru = 8.314 # J/mol.K, Universal gas constant
NucEff = 1.0E-6 # Efficiency of successful collisions
PAHC = 18 # Number of C atoms in the colliding PAH, 18 belongs to A4R5
PI = 3.14 # Pi number
Particle_Density = 1.9 # g/cc, Soot density

# The function that calculates the inception rate
def inception_rate(Xpah,P,T):
    # mole concentration of the coliding species
    MolCon = (Xpah*P)/(Ru*T)/1e6 #1e6 converts mol/m3 to mol/cc
    # inception forward rate constant
    kfr = 2.2*0.1*AV*np.sqrt(8.0*PI*BOLTZMANN*(PAHC+ PAHC)/(C_MASS*PAHC*PAHC))*\
    (2*PAH_rad+2*PAH_rad)**2*np.sqrt(T) # cc/mol.s
    return NucEff*kfr*MolCon*MolCon*AV # #/cc.s, AV converts mol/cc.s to #/cc.s

# The function that calculates the particle volume fraction (FV)
def volume_fraction(particle_number):
    # particle_number is the total number of the particles, and each particle is assumed
    # to have 500 carbon atoms. The unit of particle_number is #/cc.
    FV = particle_number*500.0*C_MASS/Particle_Density*1e6 # 1e6 converts cc/cc to ppm
    return FV

# The funciton that calculates the amount of gas that needs to be removed from the system and updates the value of 
# precursor PAH based on that for the next iteration
def gas_scrub(particle_number, X_original, P, T):
    
    # particle_number * 500/18 #/cc is the total number of PAH molecules avilable in the particles
    # * AMU converts the number of PAHs to mole/cc
    # using the formula X_i = []_i*Ru*T/P, we can estimate the mole fraction of PAHs in the partilces
    X = particle_number * 500/18 * AMU * 1e6 * Ru * T / P # 1e6 converts mol/cc to mol/m3

    # Check point to make sure that the mole fraction of PAHs that should be reomved is smaller than the avilable
    # PAHs in the system 
    if (X_original-X)>=0:
       return X_original - X # X_original is coming from the gas-phpase modeling, and X is from the particle phase,
                             # X_original - X is the amount of the remaining PAH mole fraction for the next iteration
    else:
       return 1e-10          # if X_original - X < 0, to avoid a negative mole fraction, we use 1e-10. 

