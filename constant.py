import numpy as np
from scipy.integrate import quad

NAVO = 6.02204e23 #/mol
aB = 0.5217 #Angstron
q = 1.602192e-19 #C
m0 = 9.10938356e-31 #Kg
eV = 1.602192e-19 #J
R = 1.98719 #cal/mole-K
kB = 1.380662e-23 # J/K(R/NAVO)
k = 8.617e-5 #eV kB/q

#mn = 1.18*m0 #electron effective mass at 300K
#mp = 0.81*m0 #hole effective mass at 300K
#ni=1e10/cm^3 ni=np.sqrt(Nc(T)*Nv(T))*(T/300)**1.5*np.exp(-(Eg(T)-0.00743)/(2*k*T))
# Define T, ND, NA where main python code
# default value in parentheses
def mn(T=300):
    return 1.028+(6.11e-4)*T-(3.09e-7)*T**2
def mp(T=300):
    return 0.610+(7.83e-4)*T-(4.46e-7)*T**2
def Eg(T=300):
    return 1.17-4.73e-4*(T**2)/(T+636)
def Nc(T=300):
    return 2*(mn(T)*m0*kB*T/(2*np.pi*hbar**2))**(3/2)/1e6
def Nv(T=300):
    return 2*(mp(T)*m0*kB*T/(2*np.pi*hbar**2))**(3/2)/1e6
def ni(T=300):
    return np.sqrt(Nc(T)*Nv(T))*np.exp(-(Eg(T)-0.00743)/(2*k*T))
def ni_Ge(T=300):
    return 1.76e16*T**(3/2)*np.exp(-0.392/(k*T))

#Mobility
def mun(T=300,ND=1e16):
    NDrefT=1.3e17*(T/300)**2.4
    munminT=92*(T/300)**-0.57
    mun0T=1268*(T/300)**-2.33
    anT=0.91*(T/300)**-0.146
    return munminT+mun0T/(1+(ND/NDrefT))**anT
    
def mup(T=300,NA=1e16):
    NArefT=2.35e17*(T/300)**2.4
    mupminT=54.3*(T/300)**-0.57
    mup0T=406.9*(T/300)**-2.23
    apT=0.88*(T/300)**-0.146
    return mupminT+mup0T/(1+(NA/NArefT))**apT

def conductivityn(T=300,Nnet=1e16):
    return q*mun(T,Nnet)*Nnet
def conductivityp(T=300,Nnet=1e16):
    return q*mup(T,Nnet)*Nnet

def resistivityn(T=300,Nnet=1e16):
    return 1/(q*mun(T,Nnet)*Nnet)
def resistivityp(T=300,Nnet=1e16):
    return 1/(q*mup(T,Nnet)*Nnet)

def sheet_resistancen(T=300,Nnet=1e16,d=1e-4):
    return 1/(q*mun(T,Nnet)*Nnet)*1/d
def sheet_resistancep(T=300,Nnet=1e16,d=1e-4):
    return 1/(q*mup(T,Nnet)*Nnet)*1/d

def resistancen(T=300,Nnet=1e16,L=1e-4,W=1e-4,d=1e-4):
    return L/W*1/(q*mun(T,Nnet)*Nnet)*1/d
def resistancep(T=300,Nnet=1e16,L=1e-4,W=1e-4,d=1e-4):
    return L/W*1/(q*mup(T,Nnet)*Nnet)*1/d

def Resistance(T=300, NDB=1e14, NDI=1e16, L=1, W=1, d=1):
    def conductivities(x, NDI):
        ND_profile = NDB + NDI * np.exp(-5 * x)  
        mu_n = mun(T, ND_profile)  # Assuming `mun` is a function
        return q * mu_n * ND_profile
    
    conductivity = quad(conductivities, 0, d, args=NDI)
    Res = L / W / conductivity[0]
    
    return Res



h = 6.62617e-34 #J-s
hbar = 1.05458e-34 #J-s(h/(2*pi))
Mp = 1.67264e-27 #kg Proton rest mass
c = 2.99792e10 #cm/s speed of light in vacuum
kT_q = 0.0259 #V thermal voltage at 300K
Lambda = 1.23977 #um wavelength of 1-eV quantum

e0 = 8.85418e-14 #F/cm(1/(mu0c^2))
e_si = 11.7
e_ox = 3.9
