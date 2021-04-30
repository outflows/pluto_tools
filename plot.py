import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyPLUTO as pypl
import pyPLUTO.pload as pp
import pyPLUTO.Image as img

#plutodir = os.environ['PLUTO_DIR']
#wdir = plutodir+'/My_Problems/RHD/GRB_in_AGN/2D_cart/2D_cart_01/'
wdir = os.getcwd()+'/'

def rd(filenumber,wdir):
    global D
    D = pp.pload(filenumber, w_dir=wdir)

def plot(Data, myvar, strvar, **kwargs):
    plt.clf()
    varscale = kwargs.pop('varscale', False)
    x1scale = kwargs.pop('x1scale', False)
    x2scale = kwargs.pop('x2scale', False)
    x3scale = kwargs.pop('x3scale', False)
    plane = kwargs.pop('plane', False)
    x1min = kwargs.pop('x1min', False)
    x1max = kwargs.pop('x1max', False)
    x2min = kwargs.pop('x2min', False)
    x2max = kwargs.pop('x2max', False)

    if plane == 'xz':
        label1 = r'$x\;[10^8\;\mathrm{cm}]$'
        label2 = r'$y\;[10^8\;\mathrm{cm}]$'
        if x1scale == 'lin':
            x1 = Data.x1
        elif x1scale == 'log':
            x1 = np.log10(Data.x1)
        elif x1scale == 'linlog':
            x1 = Data.x1
            for ind in x1:
                if ind >= 768:
                    x1[i] = np.log10(x1)
        if x2scale == 'lin':
            x2 = Data.x2
        elif x2scale == 'log':
            x2 = np.log10(Data.x2)
        elif x2scale == 'linlog':
            x2 = Data.x2
            for ind in x2:
                if ind >= 768:
                    x2[i] = np.log10(x2)

    if varscale == 'lin':
        myvar = myvar
    elif varscale == 'log':
        myvar = np.log10(myvar)

    cbarlabel = r'$\log(\rho)$'

    if strvar == 'rho':
        vmin = -2
        vmax = 8
    elif strvar == 'prs':
        vmin = -2
        vmax = 2
    I = img.Image()
    I.pldisplay(Data,
                myvar,
                x1 = x1,
                x2 = x2,
                vmin = vmin,
                vmax = vmax,
                label1 = label1,
                label2 = label2,
                title = r'$t=%.4g\;\mathrm{[code\,units]}$'%np.round(D.SimTime),
                cbar = (True,'vertical'),
                polar=[True,False],
                figsize=[3,15],
                x1min = x1min,
                x1max = x1max,
                x2min = x2min,
                x2max = x2max,)

    #plt.tight_layout()
    plt.savefig(wdir+'images/'+strvar+'_'+plane+'_'+str(myfilenumber).zfill(4)+'.png', dpi = 150) # Only to be saved as either .png or .jpg
    #plt.show()

first_snapshot = 0
last_snapshot = 200
plane = 'xz'

var = 'prs'
for myfilenumber in range(first_snapshot, last_snapshot+1):
    rd(myfilenumber,wdir)
    if var == 'rho':
        plot(D, D.rho, var, varscale='log', x1scale='lin', x2scale='lin', plane=plane, x1min=0, x1max=250, x2min=0, x2max=500)
    elif var == 'prs':
        plot(D, D.prs, var, varscale='log', x1scale='lin', x2scale='lin', plane=plane, x1min=0, x1max=250, x2min=0, x2max=500)

var = 'rho'
for myfilenumber in range(first_snapshot, last_snapshot+1):
    rd(myfilenumber,wdir)
    if var == 'rho':
        plot(D, D.rho, var, varscale='log', x1scale='lin', x2scale='lin', plane=plane, x1min=0, x1max=250, x2min=0, x2max=500)
    elif var == 'prs':
        plot(D, D.prs, var, varscale='log', x1scale='lin', x2scale='lin', plane=plane, x1min=0, x1max=250, x2min=0, x2max=500)
