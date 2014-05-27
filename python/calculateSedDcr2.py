import sys, os
import numpy as np
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
majorLocatorx   = MultipleLocator(3)
majorLocatory   = MultipleLocator(2)

# This one aggregates on SED, and the calculates DCR w.r.t. km15_5250.fits_g05_5470

contLevel  = 0.01
dcrLevel   = 5e-3
sedName    = "km15_5250.fits_g05_5470.gz"
airmass1   = 1.25                                                                                                                                       
dairmasses = np.arange(-0.25, 0.26, 0.05)                                                                                                               
dthetas    = np.arange(0, 181, 5)
zd1        = np.arccos(1./airmass1) * 180/np.pi

# Extracted from Winter2014 DCR3.py
def getOffset(wavelen, flux, zd):
    off = refract(wavelen*10**-3,  zd * np.pi / 180.) * 180. / np.pi * 3600. # arcsec
    return np.sum(off * flux) / np.sum(flux)
def refract(wavelength, zd, P=520.0, T=20.0, f=8.0):
    xnm1  = 10**-6 * ( 64.328 + 29498.1 / (146.0 - 1/wavelength**2) + 255.4 / (41.0 - 1/wavelength**2) )
    xnm1 *= P * (1 + (1.049 - 0.0157 * T) * 10**-6 * P) / (720.883 * (1 + 0.003661 * T))
    xnm1 -= 10**-6 * f * (0.0624 - 0.000680/wavelength**2) / (1 + 0.003661 * T)
    xn    = xnm1 + 1
    r0    = (xn**2 - 1) / (2 * xn**2)
    of    = r0 * np.tan(zd) # radians
    return of

def doit(args):
    sed, gmag, nsed = args
    if sed.startswith("bergeron"):
        pref = "wDs"
        suff = ".gz"
    elif sed.startswith("k"):
        pref = "kurucz"
        suff = ".gz"
    elif sed.startswith("burrows") or sed.startswith("L") or sed.startswith("l") or sed.startswith("m"):
        pref = "mlt"
        suff = ".gz"

    star = Sed()
    star.readSED_flambda(os.path.join(catDir, "SEDs/starSED/%s/%s%s" % (pref, sed, suff)))
    star.flambdaTofnu()
    Mg  = np.empty((len(dairmasses), len(dthetas)))
    Mr  = np.empty((len(dairmasses), len(dthetas)))
    Mi  = np.empty((len(dairmasses), len(dthetas)))
    Mz  = np.empty((len(dairmasses), len(dthetas)))
    
    # Resample array elements
    wavelen        = star.wavelen
    fnu            = star.fnu
    waveleng, fnug = star.resampleSED(wavelen, fnu, wavelen_match=gBand.wavelen)
    wavelenr, fnur = star.resampleSED(wavelen, fnu, wavelen_match=rBand.wavelen)
    waveleni, fnui = star.resampleSED(wavelen, fnu, wavelen_match=iBand.wavelen)
    wavelenz, fnuz = star.resampleSED(wavelen, fnu, wavelen_match=zBand.wavelen)
    fluxg = fnug * gBand.phi
    fluxr = fnur * rBand.phi
    fluxi = fnui * iBand.phi
    fluxz = fnuz * zBand.phi

    gr    = -2.5 * np.log10( np.sum(fluxg) / np.sum(fluxr) )
    gi    = -2.5 * np.log10( np.sum(fluxg) / np.sum(fluxi) )
    gz    = -2.5 * np.log10( np.sum(fluxg) / np.sum(fluxz) )
    rmag  = gmag - gr
    imag  = gmag - gi
    zmag  = gmag - gz

    # Find DCR at airmass zd1
    dg1   = getOffset(waveleng, fluxg, zd1) - refOffset["g"][int(zd1)]
    dr1   = getOffset(wavelenr, fluxr, zd1) - refOffset["r"][int(zd1)]
    di1   = getOffset(waveleni, fluxi, zd1) - refOffset["i"][int(zd1)]
    dz1   = getOffset(wavelenz, fluxz, zd1) - refOffset["z"][int(zd1)]
    for ia, dairmass in enumerate(dairmasses):
        # Find DCR at other zds
        zd2 = np.arccos(1./(airmass1 + dairmass)) * 180/np.pi
        dg2 = getOffset(waveleng, fluxg, zd2) - refOffset["g"][int(zd2)]
        dr2 = getOffset(wavelenr, fluxr, zd2) - refOffset["r"][int(zd2)]
        di2 = getOffset(waveleni, fluxi, zd2) - refOffset["i"][int(zd2)]
        dz2 = getOffset(wavelenz, fluxz, zd2) - refOffset["z"][int(zd2)]
        for it, dtheta in enumerate(dthetas):
            #print ia, it, (airmass1 + dairmass), zd1, zd2, dtheta, dg1, dg2
            Mg[ia, it] = np.sqrt(dg1**2 + dg2**2 - 2 * dg1 * dg2 * np.cos(dtheta * np.pi / 180))
            Mr[ia, it] = np.sqrt(dr1**2 + dr2**2 - 2 * dr1 * dr2 * np.cos(dtheta * np.pi / 180))
            Mi[ia, it] = np.sqrt(di1**2 + di2**2 - 2 * di1 * di2 * np.cos(dtheta * np.pi / 180))
            Mz[ia, it] = np.sqrt(dz1**2 + dz2**2 - 2 * dz1 * dz2 * np.cos(dtheta * np.pi / 180))
    #print Mg
    #print Mr
    if (gmag<16) or (gmag>25):
        Mg *= 0
        ng  = 0
    else:
        ng  = nsed

    if (rmag<16) or (rmag>24.7):
        Mr *= 0
        nr  = 0
    else:
        nr  = nsed

    if (imag<16) or (imag>24):
        Mi *= 0
        ni  = 0
    else:
        ni  = nsed

    if (zmag<16) or (zmag>23.3):
        Mz *= 0
        nz  = 0
    else:
        nz  = nsed
    
    return Mg, Mr, Mi, Mz, ng, nr, ni, nz
    
def integrate(sed, bp):
    wavelen = sed.wavelen
    fnu = sed.fnu
    if sed.needResample(wavelen=wavelen, wavelen_match=bp.wavelen):
        wavelen, fnu = sed.resampleSED(wavelen, fnu, wavelen_match=bp.wavelen)
    return np.sum(fnu * bp.phi)

if __name__ == "__main__":

    catDir = os.environ["CAT_SHARE_DATA"]
    filtDir = os.environ["LSST_THROUGHPUTS_BASELINE"]
    gBand  = Bandpass()
    rBand  = Bandpass()
    iBand  = Bandpass()
    zBand  = Bandpass()
    gBand.readThroughput(os.path.join(filtDir, "total_g.dat"))
    rBand.readThroughput(os.path.join(filtDir, "total_r.dat"))
    iBand.readThroughput(os.path.join(filtDir, "total_i.dat"))
    zBand.readThroughput(os.path.join(filtDir, "total_z.dat"))
    gBand.sbTophi()
    rBand.sbTophi()
    iBand.sbTophi()
    zBand.sbTophi()

    ref = Sed()
    ref.readSED_flambda(os.path.join(catDir, "SEDs/starSED/kurucz/%s"%(sedName)))
    ref.flambdaTofnu()
    wavelen        = ref.wavelen
    fnu            = ref.fnu
    waveleng, fnug = ref.resampleSED(wavelen, fnu, wavelen_match=gBand.wavelen)
    wavelenr, fnur = ref.resampleSED(wavelen, fnu, wavelen_match=rBand.wavelen)
    waveleni, fnui = ref.resampleSED(wavelen, fnu, wavelen_match=iBand.wavelen)
    wavelenz, fnuz = ref.resampleSED(wavelen, fnu, wavelen_match=zBand.wavelen)
    fluxg = fnug * gBand.phi
    fluxr = fnur * rBand.phi
    fluxi = fnui * iBand.phi
    fluxz = fnuz * zBand.phi

    refOffset = {}
    refOffset["r"] = {}
    refOffset["g"] = {}
    refOffset["i"] = {}
    refOffset["z"] = {}
    for ia, dairmass in enumerate(dairmasses):
        zd = np.arccos(1./(airmass1 + dairmass)) * 180/np.pi
        refOffset["g"][int(zd)] = getOffset(waveleng, fluxg, zd) 
        refOffset["r"][int(zd)] = getOffset(wavelenr, fluxr, zd) 
        refOffset["i"][int(zd)] = getOffset(waveleni, fluxi, zd) 
        refOffset["z"][int(zd)] = getOffset(wavelenz, fluxz, zd) 

    allg = np.zeros((len(dairmasses), len(dthetas)))
    allr = np.zeros((len(dairmasses), len(dthetas)))
    alli = np.zeros((len(dairmasses), len(dthetas)))
    allz = np.zeros((len(dairmasses), len(dthetas)))
    args = []
    seen = {}
    for line in open(sys.argv[1]).readlines()[1:]:
        sed, mag, nsed = line.split()

        # TEST: remove the M-dwarfs
        #if sed.startswith("m"):
        #    continue

        # OK, so I really want the mag in the respective passbands, so
        # don't reject here
        #mag = float(mag) 
        #if mag < 16 or mag > 25:
        #    continue

        # Also, dont group by SED then, just do it
        #if not sed in seen.keys():
        #    seen[sed]  = 0
        #seen[sed] += int(nsed)

    #for sed,nsed in seen.items():
    #    args.append((sed, mag, nsed))
    #    ntot += nsed
        
        mag = float(mag)
        nsed = int(nsed)
        args.append((sed, mag, nsed))

    import multiprocessing
    pool = multiprocessing.Pool(multiprocessing.cpu_count()-2)
    results = pool.map(doit, args)
    ntotg = ntotr = ntoti = ntotz = 0
    for i, result in enumerate(results):
        Mg, Mr, Mi, Mz, ng, nr, ni, nz = result
        allg[np.where(Mg>dcrLevel)] += ng
        allr[np.where(Mr>dcrLevel)] += nr
        alli[np.where(Mi>dcrLevel)] += ni
        allz[np.where(Mz>dcrLevel)] += nz
        ntotg += ng
        ntotr += nr
        ntoti += ni
        ntotz += nz
        print args[i][0], args[i][1], nr, nr, ni, nz
    #import pdb; pdb.set_trace()

    fig, ((sp1, sp2),(sp3, sp4)) = plt.subplots(2, 2)
    fig.suptitle("Fraction of stars showing DCR offsets of %.3f arcsec (contour at %.2f)\n w.r.t. %s" % (dcrLevel, contLevel, sedName), weight="bold", fontsize=14)

    # Contour at 1% of stars show dipoles
    hmap1 = sp1.pcolor(np.log10(allg/ntotg), cmap=plt.cm.Greys, shading="faceted", vmin=-2, vmax=0)
    sp1.contour(np.log10(allg/ntotg), levels=(np.log10(contLevel),), colors=("r"), linewidths=(3,), linestyles=("solid",))
    sp1.set_xticklabels(dthetas[::2], minor=False, weight="bold", fontsize=12)
    sp1.set_yticklabels(airmass1+dairmasses[::2], minor=False, weight="bold", fontsize=12)
    sp1.set_title("g-band (%.2e total)"%(ntotg), weight="bold", fontsize=14)
    sp1.xaxis.set_major_locator(majorLocatorx)
    sp1.yaxis.set_major_locator(majorLocatory)

    hmap2 = sp2.pcolor(np.log10(allr/ntotr), cmap=plt.cm.Greys, shading="faceted", vmin=-2, vmax=0)
    sp2.contour(np.log10(allr/ntotr), levels=(np.log10(contLevel),), colors=("r"), linewidths=(3,), linestyles=("solid",))
    sp2.set_xticklabels(dthetas[::2], minor=False, weight="bold", fontsize=12)
    sp2.set_yticklabels(airmass1+dairmasses[::2], minor=False, weight="bold", fontsize=12)
    cb2 = plt.colorbar(hmap2, ax=sp2, fraction=0.1)
    sp2.set_title("r-band (%.2e total)"%(ntotr), weight="bold", fontsize=14)
    sp2.xaxis.set_major_locator(majorLocatorx)
    sp2.yaxis.set_major_locator(majorLocatory)

    hmap3 = sp3.pcolor(np.log10(alli/ntoti), cmap=plt.cm.Greys, shading="faceted", vmin=-2, vmax=0)
    sp3.contour(np.log10(alli/ntoti), levels=(np.log10(contLevel),), colors=("r"), linewidths=(3,), linestyles=("solid",))
    sp3.set_xticklabels(dthetas[::2], minor=False, weight="bold", fontsize=12)
    sp3.set_yticklabels(airmass1+dairmasses[::2], minor=False, weight="bold", fontsize=12)
    sp3.set_title("i-band (%.2e total)"%(ntoti), weight="bold", fontsize=14)
    sp3.xaxis.set_major_locator(majorLocatorx)
    sp3.yaxis.set_major_locator(majorLocatory)

    hmap4 = sp4.pcolor(np.log10(allz/ntotz), cmap=plt.cm.Greys, shading="faceted", vmin=-2, vmax=0)
    sp4.contour(np.log10(allz/ntotz), levels=(np.log10(contLevel),), colors=("r"), linewidths=(3,), linestyles=("solid",))
    sp4.set_xticklabels(dthetas[::2], minor=False, weight="bold", fontsize=12)
    sp4.set_yticklabels(airmass1+dairmasses[::2], minor=False, weight="bold", fontsize=12)
    cb4 = plt.colorbar(hmap4, ax=sp4, fraction=0.1)
    sp4.set_title("z-band (%.2e total)"%(ntotz), weight="bold", fontsize=14)
    sp4.xaxis.set_major_locator(majorLocatorx)
    sp4.yaxis.set_major_locator(majorLocatory)

    sp1.set_ylabel("Airmass", weight="bold", fontsize=13)
    sp3.set_ylabel("Airmass", weight="bold", fontsize=13)
    sp3.set_xlabel("Delta Theta (degrees)", weight="bold", fontsize=13)
    sp4.set_xlabel("Delta Theta (degrees)", weight="bold", fontsize=13)

    plt.show()

    import pdb; pdb.set_trace()
