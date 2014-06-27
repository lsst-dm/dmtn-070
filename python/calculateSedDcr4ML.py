import sys, os
import numpy as np
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import cPickle

from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor, ExtraTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, AdaBoostRegressor

majorLocatorx   = MultipleLocator(3)
majorLocatory   = MultipleLocator(2)

# This one plots DCR (not refraction) vs. color: In partciular, DCR w.rt. 5000
# angstromgs.

contLevel  = 0.01
dcrLevel   = 5e-3
zds        = np.arange(0, 50, 1)

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

r5000      = {}
for zd in zds:
    r5000[int(zd)] = refract(500*10**-3, zd*np.pi/180.)*180./np.pi*3600. # arcsec 

def getresids(A, Y):
    cov  = np.linalg.inv(np.dot(A.T, A))
    soln = np.dot(cov, np.dot(A.T, Y))
    pred = np.dot(soln, A.T)
    return soln, Y-pred

def summstats(z, res, n):
    bins = list(set(z))
    med  = []
    std  = []
    fbad = [] # Lets just say the number of SEDs with residuals larger than dcrLevel
    for bin in bins:
        idx = np.where(z == bin)
        med.append(np.median(res[idx]))
        std.append(0.741 * (np.percentile(res[idx], 75) - np.percentile(res[idx], 25)))
        idx2 = np.where(np.abs(res[idx])>dcrLevel)
        fbad.append(np.sum(n[idx][idx2])/np.sum(n[idx]))
    return np.array(bins), np.array(med), np.array(std), np.array(fbad)

def doplot(regressor, X_test, y_test, z_test, n_test, fnames, rname, outfile):
    res  = regressor.predict(X_test) - y_test
    z, med, std, fbad = summstats(z_test, res, n_test)
    imports = regressor.feature_importances_
    idx     = np.argsort(imports)[::-1]
    cumh    = np.cumsum(imports[idx])
    y       = np.arange(len(fnames))
    fig     = plt.figure()
    plt.barh(y, imports[idx], align="center", log=True)
    plt.yticks(y, fnames[idx], fontsize=8, weight="bold")
    plt.xlabel("Importance")
    plt.ylim(0, len(fnames))
    plt.title("%s feature importance" % (rname))
    plt.ylim(-1.5, len(fnames)+0.5)
    xmin, xmax = plt.xlim()
    plt.xlim(xmin, 2)
    plt.savefig(outfile)
    return z, med, std, fbad

def doregress(X_train, y_train, n_train, X_test, y_test, n_test, band, fnames):
    lin  = LinearRegression()
    lin.fit(X_train, y_train)
    lres  = lin.predict(X_test) - y_test
    zl, ml, sl, fl = summstats(z_test, lres, n_test)
    
    #gbr1 = GradientBoostingRegressor(loss="ls")
    #gbr2 = GradientBoostingRegressor(loss="lad")
    #gbr1.fit(X_train, y_train)
    #gbr2.fit(X_train, y_train)
    #g1res = gbr1.predict(X_test) - y_test
    #g2res = gbr2.predict(X_test) - y_test
    #g1z, g1med, g1std = summstats(z_test, g1res)
    #g2z, g2med, g2std = summstats(z_test, g2res)
    
    #ada   = AdaBoostRegressor()
    #ada.fit(X_train, y_train)
    #ares  = ada.predict(X_test) - y_test
    #az, amed, astd = summstats(z_test, ares)

    # Some of these appear to be unstable
    # I.e. feature importance changes
    #for extension in ("A", "B", "C", "D", "E"):
    for extension in ("A", ):
        print "# Regressing", extension

        xtr   = ExtraTreeRegressor()
        xtr.fit(X_train, y_train)
        zx, mx, sx, fx = doplot(xtr, X_test, y_test, z_test, n_test, fnames, "%s-band ExtraTreeRegressor"%(band), "DCR_%s_%s_ext.png"%(band, extension))

        xtrw  = ExtraTreeRegressor()
        xtrw.fit(X_train, y_train, sample_weight=n_train)
        zxw, mxw, sxw, fxw = doplot(xtrw, X_test, y_test, z_test, n_test, fnames, "%s-band weighted ExtraTreeRegressor"%(band), "DCR_%s_%s_ext_weight.png"%(band, extension))

        ####

        tree = DecisionTreeRegressor()
        tree.fit(X_train, y_train)
        zt, mt, st, ft = doplot(tree, X_test, y_test, z_test, n_test, fnames, "%s-band DecisionTreeRegressor"%(band), "DCR_%s_%s_tree.png"%(band, extension))

        treew = DecisionTreeRegressor()
        treew.fit(X_train, y_train, sample_weight=n_train)
        ztw, mtw, stw, ftw = doplot(treew, X_test, y_test, z_test, n_test, fnames, "%s-band weighted DecisionTreeRegressor"%(band), "DCR_%s_%s_tree_weight.png"%(band, extension))

        ####
        weights = n_train
        nt      = 50

        rfr  = RandomForestRegressor(n_estimators=nt)
        rfr.fit(X_train, y_train)
        zr, mr, sr, fr = doplot(rfr, X_test, y_test, z_test, n_test, fnames, "%s-band RandomForestRegressor"%(band), "DCR_%s_%s_%d_rfr.png"%(band, extension, nt))
                
        rfrw  = RandomForestRegressor(n_estimators=nt)
        rfrw.fit(X_train, y_train, sample_weight=weights)
        zrw, mrw, srw, frw = doplot(rfrw, X_test, y_test, z_test, n_test, fnames, "%s-band weighted RandomForestRegressor"%(band), "DCR_%s_%s_%d_rfr_weight.png"%(band, extension, nt))
        print "RF %d : %.5e +/- %.5e vs weighted %.5e +/- %.5e" % (nt, 
                                                                   np.median(fr), 0.741 * (np.percentile(fr, 75) - np.percentile(fr, 25)),
                                                                   np.median(frw), 0.741 * (np.percentile(frw, 75) - np.percentile(frw, 25)))

     
        ####

        # Compare all models
        fig, (sp1, sp2, sp3) = plt.subplots(3, 1, sharex=True, figsize=(16,12))
        sp1.plot(zl, ml, "r-", label="LinearRegression")
        sp1.plot(zt, mt, "b-", label="DecisionTreeRegressor")
        sp1.plot(zr, mr, "g-", label="RandomForestRegressor")
        sp1.plot(zx, mx, "m-", label="ExtraTreeRegressor")

        sp2.plot(zl[np.where(sl>0.)], sl[np.where(sl>0.)], "r-")
        sp2.plot(zt[np.where(st>0.)], st[np.where(st>0.)], "b-")
        sp2.plot(zr[np.where(sr>0.)], sr[np.where(sr>0.)], "g-")
        sp2.plot(zx[np.where(sx>0.)], sx[np.where(sx>0.)], "m-")
        ymin, ymax = sp2.get_ylim()
        sp2.set_ylim(max(1e-7,ymin), 1e-1)

        sp3.plot(zl[np.where(fl>0.)], fl[np.where(fl>0.)], "r-")
        sp3.plot(zt[np.where(ft>0.)], ft[np.where(ft>0.)], "b-")
        sp3.plot(zr[np.where(fr>0.)], fr[np.where(fr>0.)], "g-")
        sp3.plot(zx[np.where(fx>0.)], fx[np.where(fx>0.)], "m-")
        ymin, ymax = sp3.get_ylim()
        sp3.set_ylim(max(1e-7,ymin), 1.1)

        sp1.legend(loc=2, fancybox=True)
        sp1.set_title("Mean DCR residual w.r.t. 500nm (arcsec)", weight="bold")
        sp2.set_ylabel("RMS residual (arcsec)", weight="bold")
        sp3.set_ylabel("f_tot with dDCR>%.3f"%(dcrLevel), weight="bold")
        sp3.set_xlabel("Zenith distance (deg)", weight="bold")
        sp1.axhline(y=0, c='k', linestyle='--', alpha=0.5)
        sp2.axhline(y=dcrLevel, c='k', linestyle='--', alpha=0.5)
        sp3.axhline(y=0.01, c='k', linestyle='--', alpha=0.5)
        sp2.semilogy()
        sp3.semilogy()
        plt.savefig("DCR_%s_%s.png" % (band, extension))

        ###

        fig, (sp1, sp2, sp3) = plt.subplots(3, 1, sharex=True, figsize=(16,12))
        sp1.plot(zl,  ml,  "r-", label="LinearRegression")
        sp1.plot(ztw, mtw, "b-", label="DecisionTreeRegressor weighted")
        sp1.plot(zrw, mrw, "g-", label="RandomForestRegressor weighted")
        sp1.plot(zxw, mxw, "m-", label="ExtraTreeRegressor weighted")

        sp2.plot(zl[np.where(sl>0.)],  sl[np.where(sl>0.)],  "r-")
        sp2.plot(ztw[np.where(stw>0.)], stw[np.where(stw>0.)], "b-")
        sp2.plot(zrw[np.where(srw>0.)], srw[np.where(srw>0.)], "g-")
        sp2.plot(zxw[np.where(sxw>0.)], sxw[np.where(sxw>0.)], "m-")
        ymin, ymax = sp2.get_ylim()
        sp2.set_ylim(max(1e-7,ymin), 1e-1)

        sp3.plot(zl[np.where(fl>0.)],  fl[np.where(fl>0.)],  "r-")
        sp3.plot(ztw[np.where(ftw>0.)], ftw[np.where(ftw>0.)], "b-")
        sp3.plot(zrw[np.where(frw>0.)], frw[np.where(frw>0.)], "g-")
        sp3.plot(zxw[np.where(fxw>0.)], fxw[np.where(fxw>0.)], "m-")
        ymin, ymax = sp3.get_ylim()
        sp3.set_ylim(max(1e-7,ymin), 1.1)

        sp1.legend(loc=2, fancybox=True)
        sp1.set_title("Mean DCR residual w.r.t. 500nm (arcsec)", weight="bold")
        sp2.set_ylabel("RMS residual (arcsec)", weight="bold")
        sp3.set_ylabel("f_tot with dDCR>%.3f"%(dcrLevel), weight="bold")
        sp3.set_xlabel("Zenith distance (deg)", weight="bold")
        sp1.axhline(y=0, c='k', linestyle='--', alpha=0.5)
        sp2.axhline(y=dcrLevel, c='k', linestyle='--', alpha=0.5)
        sp3.axhline(y=0.01, c='k', linestyle='--', alpha=0.5)
        sp2.semilogy()
        sp3.semilogy()
        plt.savefig("DCR_%s_%s_weight.png" % (band, extension))


def doit(args):
    sed = args
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
    
    # Resample array elements
    wavelen        = star.wavelen
    fnu            = star.fnu
    wavelenu, fnuu = star.resampleSED(wavelen, fnu, wavelen_match=uBand.wavelen)
    waveleng, fnug = star.resampleSED(wavelen, fnu, wavelen_match=gBand.wavelen)
    wavelenr, fnur = star.resampleSED(wavelen, fnu, wavelen_match=rBand.wavelen)
    waveleni, fnui = star.resampleSED(wavelen, fnu, wavelen_match=iBand.wavelen)
    wavelenz, fnuz = star.resampleSED(wavelen, fnu, wavelen_match=zBand.wavelen)
    fluxu = fnuu * uBand.phi
    fluxg = fnug * gBand.phi
    fluxr = fnur * rBand.phi
    fluxi = fnui * iBand.phi
    fluxz = fnuz * zBand.phi

    fluxes = [fluxu, fluxg, fluxr, fluxi, fluxz]
    nflux  = 5
    colors = []
    for i in range(5):
        for j in range(i+1, 5):
            color = -2.5 * np.log10( np.sum(fluxes[i]) / np.sum(fluxes[j]) )
            colors.append(color)

    dcr = []
    for zd in zds:
        # Find R 
        ru = getOffset(wavelenu, fluxu, zd) - r5000[int(zd)]
        rg = getOffset(waveleng, fluxg, zd) - r5000[int(zd)]
        rr = getOffset(wavelenr, fluxr, zd) - r5000[int(zd)]
        ri = getOffset(waveleni, fluxi, zd) - r5000[int(zd)]
        rz = getOffset(wavelenz, fluxz, zd) - r5000[int(zd)]
        dcr.append((zd, ru, rg, rr, ri, rz))

    return sed, colors, dcr
    
def integrate(sed, bp):
    wavelen = sed.wavelen
    fnu = sed.fnu
    if sed.needResample(wavelen=wavelen, wavelen_match=bp.wavelen):
        wavelen, fnu = sed.resampleSED(wavelen, fnu, wavelen_match=bp.wavelen)
    return np.sum(fnu * bp.phi)

if __name__ == "__main__":
    pickleFile = "calculateSedDcr4.pickle"
    if not os.path.isfile(pickleFile):
        catDir = os.environ["CAT_SHARE_DATA"]
        filtDir = os.environ["LSST_THROUGHPUTS_BASELINE"]
        uBand  = Bandpass()
        gBand  = Bandpass()
        rBand  = Bandpass()
        iBand  = Bandpass()
        zBand  = Bandpass()
        uBand.readThroughput(os.path.join(filtDir, "total_u.dat"))
        gBand.readThroughput(os.path.join(filtDir, "total_g.dat"))
        rBand.readThroughput(os.path.join(filtDir, "total_r.dat"))
        iBand.readThroughput(os.path.join(filtDir, "total_i.dat"))
        zBand.readThroughput(os.path.join(filtDir, "total_z.dat"))
        uBand.sbTophi()
        gBand.sbTophi()
        rBand.sbTophi()
        iBand.sbTophi()
        zBand.sbTophi()
    
        seen = []
        for line in open(sys.argv[1]).readlines()[1:]:
            sed, mag, nsed = line.split()

            # These seem to have funky colors
            if sed.startswith("burrows") or sed.startswith("L"):
                continue

            if not sed in seen:
                seen.append(sed)
                
        doit(seen[0])
        import multiprocessing
        pool = multiprocessing.Pool(multiprocessing.cpu_count()-2)
        results = pool.map(doit, seen)
    
        cPickle.dump(results, open(pickleFile, "wb"))
    else:
        results = cPickle.load(open(pickleFile, "rb"))


    # A little bit of processing I should have done:
    # How many objects, per passband, above the detection limit
    labels = [r[0] for r in results]
    bands  = ("u", "g", "r", "i", "z")
    cnames = []
    for i in range(5):
        for j in range(i+1,5):
            cnames.append("%s-%s" % (bands[i], bands[j]))

    print "# SORTING..."
    count  = {}
    for line in open(sys.argv[1]).readlines()[1:]:
        sed, mag, nsed = line.split()
        if sed.startswith("burrows") or sed.startswith("L"):
            continue

        if not sed in count.keys():
            count[sed] = {}
            for b in bands:
                count[sed][b] = 0

        idx    = labels.index(sed)
        colors = results[idx][1]
        ug     = colors[0]
        gr     = colors[4]
        ri     = colors[7]
        iz     = colors[9]
        gmag   = float(mag) + 1.0 # THERE WAS A BUG IN THE MAGS, 1 MAG TOO BRIGHT
        umag   = gmag + ug
        rmag   = gmag - gr
        imag   = rmag - ri
        zmag   = imag - iz
        if (umag>16) and (umag<23.9): count[sed]["u"] += int(nsed)
        if (gmag>16) and (gmag<25.0): count[sed]["g"] += int(nsed)
        if (rmag>16) and (rmag<24.7): count[sed]["r"] += int(nsed)
        if (imag>16) and (imag<24.0): count[sed]["i"] += int(nsed)
        if (zmag>16) and (zmag<23.3): count[sed]["z"] += int(nsed)

    # DCR in a given band
    bestResults = {}
    ncolor      = len(cnames)
    nzd         = len(zds)
    nresult     = len(results)
    fnames      = []
    for band, idx in zip(bands, (1,2,3,4,5)):
        print "DCR IN %s BAND" % (band)

        X = np.empty((nresult*nzd, 3*ncolor+2))
        y = np.empty((nresult*nzd))
        z = np.empty((nresult*nzd))
        n = np.empty((nresult*nzd))
        # Build a model with ZD, color, and color*ZD
        for r, result in enumerate(results):
            colors = result[1]
            zd     = np.array([x[0]   for x in result[2]])
            dcr    = np.array([x[idx] for x in result[2]])
            tanz   = np.tan(zd*np.pi/180.) 
            X[r*nzd:(r+1)*nzd,:ncolor]           = colors
            X[r*nzd:(r+1)*nzd,ncolor:2*ncolor]   = np.outer(tanz, colors)
            X[r*nzd:(r+1)*nzd,2*ncolor:3*ncolor] = np.outer(tanz**3, colors)
            X[r*nzd:(r+1)*nzd,3*ncolor]          = tanz * np.ones(nzd)
            X[r*nzd:(r+1)*nzd,3*ncolor+1]        = tanz**3 * np.ones(nzd)
            y[r*nzd:(r+1)*nzd]                   = dcr
            z[r*nzd:(r+1)*nzd]                   = zd
            n[r*nzd:(r+1)*nzd]                   = count[result[0]][band]

            if len(fnames) == 0:
                [fnames.append(x) for x in cnames]
                [fnames.append("tanz "+x) for x in cnames]
                [fnames.append("tanz^3 "+x) for x in cnames]
                fnames.append("tanz")
                fnames.append("tanz^3")
                fnames = np.array(fnames)
                
        X_train = np.vstack((X[0::3,:], X[1::3,:]))
        y_train = np.hstack((y[0::3],   y[1::3]))
        n_train = np.hstack((n[0::3],   n[1::3]))
        X_test  = X[2::3,:]
        y_test  = y[2::3]
        z_test  = z[2::3]
        n_test  = n[2::3]

        # get rid of bins with no stars
        idx     = np.where(n_train > 0)[0]
        X_train = X_train[idx,:]
        y_train = y_train[idx]
        n_train = n_train[idx]

        idx     = np.where(n_test > 0)[0]
        X_test  = X_test[idx,:]
        y_test  = y_test[idx]
        z_test  = z_test[idx]
        n_test  = n_test[idx]
        
        doregress(X_train, y_train, n_train, X_test, y_test, n_test, band, fnames)
