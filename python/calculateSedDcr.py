import sys, os
import numpy as np
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import matplotlib.pyplot as plt

def doit(args):
    line, gBand, rBand, catDir = args
    sed, mag, nsed = line.split()
    nsed = int(nsed)
    
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
    gflux = integrate(star, gBand)
    rflux = integrate(star, rBand)
    color = -2.5 * np.log10(gflux/rflux)
    #if int(color*100)==45:
    #    print sed, nsed
    return sed, color, nsed
    
def integrate(sed, bp):
    wavelen = sed.wavelen
    fnu = sed.fnu
    if sed.needResample(wavelen=wavelen, wavelen_match=bp.wavelen):
        wavelen, fnu = sed.resampleSED(wavelen, fnu, wavelen_match=bp.wavelen)
    return np.sum(fnu * bp.phi)

if __name__ == "__main__":
    import multiprocessing
    pool = multiprocessing.Pool(multiprocessing.cpu_count()-2)
    
    filtDir = os.environ["LSST_THROUGHPUTS_BASELINE"]
    gBand  = Bandpass()
    rBand  = Bandpass()
    gBand.readThroughput(os.path.join(filtDir, "total_g.dat"))
    rBand.readThroughput(os.path.join(filtDir, "total_r.dat"))
    gBand.sbTophi()
    rBand.sbTophi()
    
    catDir = os.environ["CAT_SHARE_DATA"]
    
    args = []
    for line in open(sys.argv[1]).readlines()[1:]:
        sed, mag, nsed = line.split()
        mag = float(mag) 
        if mag < 16 or mag > 25:
            continue
        args.append((line, gBand, rBand, catDir))

    results = pool.map(doit, args)

    grs = []
    for r in results:
        grs += r[2]*[r[1],]
    import pdb; pdb.set_trace()
    caw = np.array(grs[::1000])
    plt.hist(caw, bins=100, log=True)
    plt.xlabel("g-r")
    plt.ylabel("N")
    plt.show()
    # Median of all data is 0.51
    # Median of g-r<1.3 is 0.45
    
    # The most stars with this color are 
    # km15_5250.fits_g05_5470 1023230
    # km10_5500.fits_g10_5520 1062737
    # km10_5500.fits_g10_5520 1078438
    # km10_5500.fits_g10_5520 1100369
    # km15_5250.fits_g05_5470 1109509
    # km15_5250.fits_g05_5470 1157465
    # Lets assume that km15_5250.fits_g05_5470 defines our reference grid


