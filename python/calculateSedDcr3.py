import sys, os
import numpy as np
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import cPickle

majorLocatorx   = MultipleLocator(3)
majorLocatory   = MultipleLocator(2)

# This one plots refraction vs. color: What do I need to measure to
# estimate refraction?

contLevel  = 0.01
dcrLevel   = 5e-3
sedName    = "km15_5250.fits_g05_5470.gz"
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

def getresids(A, Y):
    cov  = np.linalg.inv(np.dot(A.T, A))
    soln = np.dot(cov, np.dot(A.T, Y))
    pred = np.dot(soln, A.T)
    return soln, Y-pred

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

    refrac = []
    for zd in zds:
        # Find R 
        ru = getOffset(waveleng, fluxu, zd)
        rg = getOffset(waveleng, fluxg, zd)
        rr = getOffset(wavelenr, fluxr, zd)
        ri = getOffset(waveleni, fluxi, zd)
        rz = getOffset(wavelenz, fluxz, zd)
        refrac.append((zd, ru, rg, rr, ri, rz))

    return sed, colors, refrac
    
def integrate(sed, bp):
    wavelen = sed.wavelen
    fnu = sed.fnu
    if sed.needResample(wavelen=wavelen, wavelen_match=bp.wavelen):
        wavelen, fnu = sed.resampleSED(wavelen, fnu, wavelen_match=bp.wavelen)
    return np.sum(fnu * bp.phi)

if __name__ == "__main__":
    pickleFile = "calculateSedDcr3.pickle"
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
            # Leading to u-band refraction RMS residuals of 0.06-0.07" 
            # Without it is 0.01"
            if sed.startswith("burrows") or sed.startswith("L"):
                continue

            if not sed in seen:
                seen.append(sed)
                
        import multiprocessing
        pool = multiprocessing.Pool(multiprocessing.cpu_count()-2)
        results = pool.map(doit, seen)
    
        cPickle.dump(results, open(pickleFile, "wb"))
    else:
        results = cPickle.load(open(pickleFile, "rb"))

    bands  = ("u", "g", "r", "i", "z")
    bcolor = ("b", "g", "r", "m", "k")
    colors = []
    for i in range(5):
        for j in range(i+1,5):
            colors.append("%s-%s" % (bands[i], bands[j]))

    # Refraction in a given band
    bestResults = {}
    for band, idx in zip(bands, (1,2,3,4,5)):
        print "REFRACTION IN %s BAND" % (band)

        # Find out which color works best in predicting the refraction for each passband
        oresults = {}

        #for cidx in range(10):
        # None do, composite colors basically win
        for cidx in range(0):

            # NOTE: the terms with pure color do not seem to work, the
            # solution is biased with mean residuals != 0.0

            # Model12 is refraction = A * tanz + B * tanz**2
            A12 = np.empty((0,2))
            Y12 = np.empty((0,))
            # Model13 is refraction = A * tanz + B * tanz**3
            A13 = np.empty((0,2))
            Y13 = np.empty((0,))

            # Model2 is refraction = A * tanz + C * tanz * color
            A2  = np.empty((0,2))
            Y2  = np.empty((0,))

            # Model32 is refraction = A * tanz + B * tanz**2 + C * tanz * color
            A32 = np.empty((0,3))
            Y32 = np.empty((0,))
            # Model33 is refraction = A * tanz + B * tanz**3 + C * tanz * color
            A33 = np.empty((0,3))
            Y33 = np.empty((0,))
            # Model33 is refraction = A * tanz + B * tanz**3 + C * tanz**5 + D * tanz * color
            A35 = np.empty((0,4))
            Y35 = np.empty((0,))

            # Model32 is refraction = A * tanz + B * tanz**2 + C * tanz * color + D * tanz**2 * color
            A42 = np.empty((0,4))
            Y42 = np.empty((0,))
            # Model33 is refraction = A * tanz + B * tanz**3 + C * tanz * color + D * tanz**3 * color
            A43 = np.empty((0,4))
            Y43 = np.empty((0,))
            # Model33 is refraction = A * tanz + B * tanz**3 + C * tanz**3 + D * tanz * color + E * tanz**3 * color + F * tanz**3 * color
            A45 = np.empty((0,6))
            Y45 = np.empty((0,))


            for result in results:
                color = result[1][cidx]
                zd    = np.array([x[0]   for x in result[2]])
                ref   = np.array([x[idx] for x in result[2]])
                tanz  = np.tan(zd*np.pi/180.) 

                if ref[-1] < 44.5 and band == "u":
                    print result[0]

                A12  = np.vstack((A12, np.vstack((tanz, tanz**2)).T))
                Y12  = np.append(Y12, ref)
                A13  = np.vstack((A13, np.vstack((tanz, tanz**3)).T))
                Y13  = np.append(Y13, ref)

                A2   = np.vstack((A2, np.vstack((tanz, tanz*color)).T))
                Y2   = np.append(Y2, ref)

                A32  = np.vstack((A32, np.vstack((tanz, tanz**2, tanz*color)).T))
                Y32  = np.append(Y32, ref)
                A33  = np.vstack((A33, np.vstack((tanz, tanz**3, tanz*color)).T))
                Y33  = np.append(Y33, ref)
                A35  = np.vstack((A35, np.vstack((tanz, tanz**3, tanz**5, tanz*color)).T))
                Y35  = np.append(Y35, ref)

                A42  = np.vstack((A42, np.vstack((tanz, tanz**2, tanz*color, tanz**2*color)).T))
                Y42  = np.append(Y42, ref)
                A43  = np.vstack((A43, np.vstack((tanz, tanz**3, tanz*color, tanz**3*color)).T))
                Y43  = np.append(Y43, ref)
                A45  = np.vstack((A45, np.vstack((tanz, tanz**3, tanz**5, tanz*color, tanz**3*color, tanz**5*color)).T))
                Y45  = np.append(Y45, ref)

            soln12, resid12 = getresids(A12, Y12)
            soln13, resid13 = getresids(A13, Y13)
            soln2,  resid2  = getresids(A2, Y2)
            soln32, resid32 = getresids(A32, Y32)
            soln33, resid33 = getresids(A33, Y33)
            soln35, resid35 = getresids(A35, Y35)
            soln42, resid42 = getresids(A42, Y42)
            soln43, resid43 = getresids(A43, Y43)
            soln45, resid45 = getresids(A45, Y45)

            #print " %s model 1 : resids %7.4f +/- %5.4f as" % (colors[cidx], np.mean(resid1), np.std(resid1)), soln1
            #print " %s model 2 : resids %7.4f +/- %5.4f as" % (colors[cidx], np.mean(resid2), np.std(resid2)), soln2
            #print " %s model 3 : resids %7.4f +/- %5.4f as" % (colors[cidx], np.mean(resid3), np.std(resid3)), soln3
            #print " %s model 4 : resids %7.4f +/- %5.4f as" % (colors[cidx], np.mean(resid4), np.std(resid4)), soln4
            #print " %s model 5 : resids %7.4f +/- %5.4f as" % (colors[cidx], np.mean(resid5), np.std(resid5)), soln5

            oresults[np.std(resid12)] =  (" %s model 1_2 : resids %8.5f +/- %6.5f as" % (colors[cidx], np.mean(resid12), np.std(resid12)), resid12)
            oresults[np.std(resid13)] =  (" %s model 1_3 : resids %8.5f +/- %6.5f as" % (colors[cidx], np.mean(resid13), np.std(resid13)), resid13)
            oresults[np.std(resid2)]  =  (" %s model 2   : resids %8.5f +/- %6.5f as" % (colors[cidx], np.mean(resid2),  np.std(resid2)),  resid2)
            oresults[np.std(resid32)] =  (" %s model 3_2 : resids %8.5f +/- %6.5f as" % (colors[cidx], np.mean(resid32), np.std(resid32)), resid32)
            oresults[np.std(resid33)] =  (" %s model 3_3 : resids %8.5f +/- %6.5f as" % (colors[cidx], np.mean(resid33), np.std(resid33)), resid33)
            oresults[np.std(resid35)] =  (" %s model 3_5 : resids %8.5f +/- %6.5f as" % (colors[cidx], np.mean(resid35), np.std(resid35)), resid35)
            oresults[np.std(resid42)] =  (" %s model 4_2 : resids %8.5f +/- %6.5f as" % (colors[cidx], np.mean(resid42), np.std(resid42)), resid42)
            oresults[np.std(resid43)] =  (" %s model 4_3 : resids %8.5f +/- %6.5f as" % (colors[cidx], np.mean(resid43), np.std(resid43)), resid43)
            oresults[np.std(resid45)] =  (" %s model 4_5 : resids %8.5f +/- %6.5f as" % (colors[cidx], np.mean(resid45), np.std(resid45)), resid45)

        # Finally, use all the colors to see which are useful
        # Model6 is refraction = A * tanz + Sum_i ( B_i tanz * color_i )
        # This does not look good; ignore
        #A6 = np.empty((0,6))
        #Y6 = np.empty((0,))

        # Model 7: Only u-g, g-r
        A7  = np.empty((0,3))
        Y7  = np.empty((0,))
        A73 = np.empty((0,6))
        Y73 = np.empty((0,))

        # Model 8: Only g-r, r-i
        A8  = np.empty((0,3))
        Y8  = np.empty((0,))
        A83 = np.empty((0,6))
        Y83 = np.empty((0,))

        # Model 9: Only r-i, i-z
        A9  = np.empty((0,3))
        Y9  = np.empty((0,))
        A93 = np.empty((0,6))
        Y93 = np.empty((0,))

        for result in results:
            zd    = np.array([x[0]   for x in result[2]])
            ref   = np.array([x[idx] for x in result[2]])
            tanz  = np.tan(zd*np.pi/180.) 
            #A6    = np.vstack((A6, np.vstack((tanz, tanz*result[1][0], tanz*result[1][4], tanz*result[1][5], tanz*result[1][7], tanz*result[1][9])).T))
            #Y6    = np.append(Y6, ref)

            A7    = np.vstack((A7, np.vstack((tanz, tanz*result[1][0], tanz*result[1][4])).T))
            Y7    = np.append(Y7, ref)
            A73   = np.vstack((A73, np.vstack((tanz, tanz*result[1][0], tanz*result[1][4], tanz**3, tanz**3*result[1][0], tanz**3*result[1][4])).T))
            Y73   = np.append(Y73, ref)

            A8    = np.vstack((A8, np.vstack((tanz, tanz*result[1][4], tanz*result[1][7])).T))
            Y8    = np.append(Y8, ref)
            A83   = np.vstack((A83, np.vstack((tanz, tanz*result[1][4], tanz*result[1][7], tanz**3, tanz**3*result[1][4], tanz**3*result[1][7])).T))
            Y83   = np.append(Y83, ref)

            A9    = np.vstack((A9, np.vstack((tanz, tanz*result[1][7], tanz*result[1][9])).T))
            Y9    = np.append(Y9, ref)
            A93   = np.vstack((A93, np.vstack((tanz, tanz*result[1][7], tanz*result[1][9], tanz**3, tanz**3*result[1][7], tanz**3*result[1][9])).T))
            Y93   = np.append(Y93, ref)

        #soln6, resid6 = getresids(A6, Y6)
        soln7,  resid7  = getresids(A7, Y7)
        soln8,  resid8  = getresids(A8, Y8)
        soln9,  resid9  = getresids(A9, Y9)
        soln73, resid73 = getresids(A73, Y73)
        soln83, resid83 = getresids(A83, Y83)
        soln93, resid93 = getresids(A93, Y93)
        #oresults[np.std(resid6)] =  (" %s model 6   : resids %8.5f +/- %6.5f as" % ("   ", np.mean(resid6),  np.std(resid6)),  resid6)
        oresults[np.std(resid7)]  =  (" %s model 7   : resids %8.5f +/- %6.5f as" % ("   ", np.mean(resid7),  np.std(resid7)),  resid7)
        oresults[np.std(resid73)] =  (" %s model 7_3 : resids %8.5f +/- %6.5f as" % ("   ", np.mean(resid73), np.std(resid73)), resid73)
        oresults[np.std(resid8)]  =  (" %s model 8   : resids %8.5f +/- %6.5f as" % ("   ", np.mean(resid8),  np.std(resid8)),  resid8)
        oresults[np.std(resid83)] =  (" %s model 8_3 : resids %8.5f +/- %6.5f as" % ("   ", np.mean(resid83), np.std(resid83)), resid83)
        oresults[np.std(resid9)]  =  (" %s model 9   : resids %8.5f +/- %6.5f as" % ("   ", np.mean(resid9),  np.std(resid9)),  resid9)
        oresults[np.std(resid93)] =  (" %s model 9_3 : resids %8.5f +/- %6.5f as" % ("   ", np.mean(resid93), np.std(resid93)), resid93)


        ### SORT

        keys = oresults.keys()
        keys.sort()
        for key in keys:
            print oresults[key][0]

        bestResults[band] = oresults[keys[0]][1]

    fig = plt.figure()
    # resids are in blocks of 50
    for bidx, band in enumerate(bands):
        ys  = []
        dys = []
        for i in range(len(zds)):
            data = bestResults[band][i::len(zds)]
            med  = np.median(data)
            std  = 0.741 * (np.percentile(data, 75) - np.percentile(data, 25))
            ys.append(med)
            dys.append(std)
        plt.plot(zds, dys, "%s-" % (bcolor[bidx]), label=band)
    plt.legend()
    plt.xlabel("Zenith Distance", weight="bold", fontsize=14)
    plt.ylabel("RMS around Refraction fit", weight="bold", fontsize=14)
    plt.axhline(y=0.005, c='k', linestyle=':')
    plt.setp(plt.gca().get_xticklabels()+plt.gca().get_yticklabels(), weight="bold", fontsize=12)
    plt.semilogy()
        
    import pdb; pdb.set_trace()

        #plt.figure()
        #plt.plot(A5[:,0], resid5, "r.")
        #plt.figure()
        #plt.plot(A5[:,0], Y5, "r.")
        #plt.show()

        
        
