#!/usr/bin/python
# constants.py
from numpy import pi

# Some physical constants and conversion factors
PI       = pi              # circle perimeter to diameter ratio
AMH      = 1.66e-24        # hydrogen mass (g)
KB       = 1.38e-16        # Boltzmann constant (cgs)
RG       = 8.3145e7        # Gas constant (cgs)
GGRAV    = 6.67259e-8      # Gravitational constant (cgs)
CLIGHT   = 2.99E10         # speed of light in vacuum (cm/s)
ECHARGE  = 4.8032e-10      # electron charge statcoulomb (cgs)
EMASS    = 9.10938e-28     # electron mass (g)
SIGMA_SB = 5.6704e-5       # Stephan Boltzmann constant (cm^2)
SIGMA_T  = 6.65245e-25     # Thompson-scat cross-section (cm^2)

MSUN     = 1.99E33         # solar radius (cgs)
RSUN     = 6.955e10        # solar mass (cgs)
GSUN     = 274.e2          # solar gravity (cgs)
MJUP     = 1.898E30        # Jupiter mass (cgs)
RJUP     = 7.1492E9        # Jupiter radius (cgs)

AU       = 1.496e13        # 1AU in cm
PC       = 3.0857E18       # 1pc in cm
KPC      = 3.0857E21       # 1Kpc in cm
DEG      = pi/180.         # conversion from deg to rad
HR       = 3600.           # 1hr in seconds
DAY      = 86400.          # 1day in seconds
YR       = 3.1536E7        # 1yr in seconds
MYR      = 3.1536E13       # 1Myr in seconds
EV       = 1.60218E-12     # 1 ev in ergs