"""
    This file is to calculate the initial conditions for the jet and the ambient medium that need to be set in PLUTO.

    This has been adapted from Tyler's file.
"""

import numpy as np

# all constants in cgs units
ME = 9.1093826e-28 # electron mass in g
MP = 1.67262171e-24 # proton mass in g
MH = 1.673534e-24 # hydrogen atom mass in g
CL = 2.99792458e10 # speed of light in cm s^-1
GNEWT = 6.6742e-8 # gravitational constant in
KBOL = 1.3806505e-16 # Boltzmann constant
SIGMA_THOMSON = 0.665245873e-24 # Thomson cross-section
MSUN = 1.989e33 # solar mass in g
LSUN = 3.827e33 # solar luminosity
MINUTE = 60
HOUR = 60*MINUTE
DAY = 24*HOUR
WEEK = 7*DAY
YEAR = 365*DAY # seconds in a year
PC = 3.086e18 # parsec
KPC = 1000*PC # kpc

def gamma(velocity, divided_by_c=False):
    if divided_by_c:
        val=1/np.sqrt(1-velocity**2)
    else:
        val= 1 / np.sqrt(1 - (velocity/CL) ** 2)
    return val


def calcJetDensPres(lumi_total, gamma_inj, gamma_infinity, opening_angle, inner_jet_radius, unit_dens=1, unit_pres=CL**2):
    # calculates the density and pressure in real and PLUTO units (same as FLASH units where density
    # is normalized by 1 and pressure is normalized by c**2
    # opening angle should be in degrees
    #for 16TI sim do: calcJetDensPres(1.07e51, 5, 400, 10, 1e9)
    #for 20sp_down ic.calcJetDensPres(2*5.33e50, 5, 400, 10, 1e8) for first period and for FALLBACK_ENVELOPE_20SP_DOWN

    dens=lumi_total/(2*np.pi*CL**3*inner_jet_radius**2*gamma_inj*gamma_infinity*(1-np.cos(opening_angle*np.pi/180)))

    pres=(0.25*dens*CL**2)*((gamma_infinity/gamma_inj)-1)

    print("Density [g cm^-3] = %e"%(dens))
    print("Density [code] = %e"%(dens/unit_dens))
    print("Pressure [g cm^-1 s^-2] = %e"%(pres))
    print("Pressure [code] = %e"%(pres/unit_pres))

    return


def ambientMediumDens(domain_size, unit_dens=1, unit_pres=CL**2):
    #domain size must be in cm, prints the min density to have a optical depth of at most 1
    dens = MP/(SIGMA_THOMSON*domain_size)
    pres = dens*CL*CL*1.e-3

    print("Density [g cm^-3] = %e"%(dens))
    print("Density [code] = %e"%(dens/unit_dens))
    print("Pressure [g cm^-1 s^-2] = %e"%(pres))
    print("Pressure [code] = %e"%(pres/unit_pres))

    print("The density above is the maximum density so that the optical depth is ≤ 1.")
    print("The pressure is assumed to be 1000 times smaller than the density*c^2.")

#    print ('The ambient medium density has to be equal to or smaller than %e to have an optical depth of at \
#    most 1\nWith a ratio of 1e-3 between pressure and density times c^2, the pressure needs to be %e, or in code units \
#    %e'%(dens, dens*c_light**2*1e-3, dens*1e-3))

    return

def RIAF_NY_Quantities(f_adv, gamma_g, beta_m, alpha, m, mdot, rstar, unit_dens, unit_pres):
    # Abramowicz & Fragile 2013, eq. 100
    # m in solar masses
    # mdot in units of MdotEdd
    # rstar in units of Schwarzschild radius
    # f_adv < 1
    # gamma_g

    epsilonp = (1./f_adv)*((5./3. - gamma_g)/(gamma_g - 1))
    g = np.sqrt(1 + (18*alpha**2)/(5+2*epsilonp)**2) - 1

    c1 = (5+2*epsilonp)*g/(3.*alpha**2)
    c2 = np.sqrt(2*epsilonp*(5+2*epsilonp)*g/(9*alpha**2))
    c3 = 2*(5+2*epsilonp)*g/(9*alpha**2)

    v = -3.*10**10*alpha*c1*rstar**(-1./2.)
    Omega = 2.03*10**5*c2*(1./m)*rstar**(-3./2.)
    cs2 = 9.*10**20*c2*(1./rstar)
    rho = 1.07*10**-5*(1./alpha)*(1./c1)*c3**(-1./2.)*(1./m)*mdot*rstar**(-3./2.)
    prs = 9.67*10**15*(1./alpha)*(1./c1)*c3**(1./2.)*(1./m)*mdot*rstar**(-5./2.)
    B = 4.93*10**8*alpha**(-1./2.)*(1-beta_m)**(1./2.)*c1**(-1./2.)*c3**(1./4.)*m**(-1./2.)*mdot**(1./2.)*rstar**(-5./4.)
    qplus = 2.94*10**21*epsilonp*c3**(1./2.)*m**(-2.)*mdot*rstar**(-4.)
    tau_es = 1.75*(1./alpha)*(1./c1)*mdot*rstar**(-1./2.)

    print("Radial infall velocity [cm s^-1] = %e"%(v))
    print("Angular velocity [s^-1] = %e"%(Omega))
    print("Speed of sound, squared [cm s^-2] = %e"%(cs2))
    print("Density [g cm^-3] = %e"%(rho))
    print("Density [code] = %e"%(rho/unit_dens))
    print("Pressure [g cm^-1 s^-2] = %e"%(prs))
    print("Pressure [code] = %e"%(prs/unit_pres))
    print("Magnetic field [G] = %e"%(B))
    print("Viscous dissipation of energy per unit volume [erg cm^-3 s^-1] = %e"%(qplus))
    print("Tau [g cm^-1 s^-s] = %e"%(tau_es))


def RIAF_YN_PressureDensity(MBH, Mdot, R, alpha, s, unit_dens=1, unit_pres=CL**2):
    # mbh in solar masses
    # mdot in units of MdotEdd
    # r in units of Schwarzschild radius

    #m = MBH/MSUN
    #LEdd = 3.2*10000.*m*LSUN
    #MdotEdd = LEdd/0.1/CL/CL
    #mdot = Mdot/MdotEdd
    #RSch = 2.*GNEWT*MBH*MSUN/CL/CL # cm
    #r = R/RSch

    m = MBH
    r = R
    mdot = Mdot

    ne = 6.3*10**19*(1./alpha)*(1./m)*mdot*r**(-3./2. + s)
    dens = ne*ME
    dens_code = dens/unit_dens
    pres = 1.7*10**16*(1./alpha)*(1./m)*mdot*r**(-5./2. + s)
    pres_code = pres/unit_pres

    print("Density [cgs] = %e"%(dens))
    print("Density [code] = %e"%(dens/unit_dens))
    print("Pressure [cgs] = %e"%(pres))
    print("Pressure [code] = %e"%(pres/unit_pres))


def TQM_PressureDensity(r, MBH, sigma, cs, f_g, unit_dens, unit_pres):
    # radius r in Schwarzschild radii
    # velocity dispersion sigma in km/s
    # black hole mass MBH in MSUN
    # speed of sound cs in km/s

    Q = 1.

    r_Sch = 2.*GNEWT*MBH*MSUN/CL/CL # Schwarzschild radius in cm
    r_cm = r*r_Sch # r in cm
    r_km = r_cm*1.e-5 # r in km
    r_kpc = r_cm/KPC # r in kpc

    sigma200 = sigma/200
    Omega = np.sqrt(2)*sigma/r_km

    rho = 170.*sigma200**2*r_kpc**(-2.)/Q*MP # density in disc at r

    f_g_den = 2.*cs/Q/sigma
    f_g_actual = f_g/f_g_den

    #f_g = 2.*cs/Q/sigma # gas fraction
    h = f_g_actual*Q/2**(3./2.)*r_kpc # disc height in kpc
    h = h*KPC # disc height in cm

    prs = rho*h**2*Omega**2

    print("Density [g cm^-3] = %e"%(rho))
    print("Density [code] = %e"%(rho/unit_dens))
    print("Pressure [g cm^-1 s^-2] = %e"%(prs))
    print("Pressure [code] = %e"%(prs/unit_pres))


def SG_PressureDensity(MBH, Mdot, r, alpha, b, rmin, m):

    Omega = (GNEWT*MBH/r**3)**(-1./2.)
    MdotPrime = Mdot*(1. - np.sqrt(rmin/r))

    Teff = ((3./(8.*np.pi()))*MdotPrime*Omega/SIGMA_THOMSON)**(1./4.)
