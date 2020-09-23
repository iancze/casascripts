import casatasks


def fit_gaussian(imagename, region, dooff=False):
    """
    Wrapper for imfit in CASA to fit a single Gaussian component to a selected
    region of the image.

    Args:
        imagename (string): Name of CASA image (ending in .image)
        region (string): CASA region format, e.g.,
            'circle[[200pix, 200pix], 3arcsec]'
        dooff (bool): allow for fitting a zero-level offset
    """
    imfitdict = imfit(imagename=imagename, region=region, dooff=dooff)
    print(imfitdict)
    # Check if the source was resolved
    was_resolved = not imfitdict["deconvolved"]["component0"]["ispoint"]
    # Get the coordinate system
    coordsystem = imfitdict["deconvolved"]["component0"]["shape"]["direction"]["refer"]
    # Get the parameters
    headerlist = imhead(imagename)
    phasecenter_ra, phasecenter_dec = headerlist["refval"][:2]
    peak_ra = imfitdict["deconvolved"]["component0"]["shape"]["direction"]["m0"][
        "value"
    ]
    peak_dec = imfitdict["deconvolved"]["component0"]["shape"]["direction"]["m1"][
        "value"
    ]
    xcen, ycen = headerlist["refpix"][:2]
    deltax, deltay = headerlist["incr"][:2]
    peak_x = xcen + np.unwrap(np.array([0, peak_ra - phasecenter_ra]))[
        1
    ] / deltax * np.cos(phasecenter_dec)
    peak_y = ycen + (peak_dec - phasecenter_dec) / deltay

    # get the total flux

    # Print
    if coordsystem == "J2000":
        print(
            "#Peak of Gaussian component identified with imfit: J2000 %s"
            % au.rad2radec(imfitdict=imfitdict, hmsdms=True, delimiter=" ")
        )
    elif coordsystem == "ICRS":
        print(
            "#Peak of Gaussian component identified with imfit: ICRS %s"
            % au.rad2radec(imfitdict=imfitdict, hmsdms=True, delimiter=" ")
        )
        J2000coords = au.ICRSToJ2000(au.rad2radec(imfitdict=imfitdict, delimiter=" "))
        print("#Peak in J2000 coordinates: %s" % J2000coords)
    else:
        print(
            "#If the coordinates aren't in ICRS or J2000, then something \
        weird is going on"
        )
    # If the object was resolved, print the inclination, PA, major and minor axis
    if was_resolved:
        PA = imfitdict["deconvolved"]["component0"]["shape"]["positionangle"]["value"]
        majoraxis = imfitdict["deconvolved"]["component0"]["shape"]["majoraxis"][
            "value"
        ]
        minoraxis = imfitdict["deconvolved"]["component0"]["shape"]["minoraxis"][
            "value"
        ]

        print("#PA of Gaussian component: %.2f deg" % PA)
        print(
            "#Inclination of Gaussian component: %.2f deg"
            % (np.arccos(minoraxis / majoraxis) * 180 / np.pi,)
        )
    print("#Pixel coordinates of peak: x = %.3f y = %.3f" % (peak_x, peak_y))


def plot_deprojected(
    filelist, incl=0, PA=0, offx=0, offy=0, fluxscale=None, uvbins=None, show_err=True
):
    """
    Plots real and imaginary deprojected visibilities from a list of .npz files

    Args:
        filelist (list): names of .npz files storing visibility data
        incl (float): Inclination of disk [degrees]
        PA (float): Position angle of disk [degrees]
        offx (float): Horizontal offset of disk center from phase center
            [arcseconds]
        offy (float): Vertical offset of disk center from phase center
            [arcseconds]
        fluxscale (list): scaling factors to multiply the visibility values
            by before plotting. Default value is set to all ones.
        uvbins: Array of bins at which to plot the visibility values, in lambda.
            By default, the range plotted will be from 10 to 1000 kilolambda
        show_err (bool): plot error bars.

    More information about the deprojection is in A. M. Hughes et al. 2007 and
    Lay et al. 1997, Fig 3e.
    """
    if fluxscale is None:
        fluxscale = np.ones(len(filelist))
    assert len(filelist) == len(fluxscale)

    if uvbins is None:
        uvbins = 10.0 + 10.0 * np.arange(100)

    minvis = np.zeros(len(filelist))
    maxvis = np.zeros(len(filelist))
    fig, ax = plt.subplots(2, 1, sharex=True)
    for i, filename in enumerate(filelist):

        # read in the data
        inpf = np.load(filename)
        u = inpf["u"]
        v = inpf["v"]
        vis = fluxscale[i] * inpf["Vis"]
        wgt = inpf["Wgt"]

        # deproject the visibilities and do the annular averaging
        vp = deproject_vis(
            [u, v, vis, wgt], bins=uvbins, incl=incl, PA=PA, offx=offx, offy=offy
        )
        vp_rho, vp_vis, vp_sig = vp

        # calculate min, max of deprojected, averaged reals (for visualization)
        minvis[i] = np.min(vp_vis.real)
        maxvis[i] = np.max(vp_vis.real)

        # plot the profile
        if show_err:
            ax[0].errorbar(
                1e-3 * vp_rho, vp_vis.real, yerr=vp_sig.real, label=filename, fmt="."
            )
            ax[1].errorbar(
                1e-3 * vp_rho, vp_vis.imag, yerr=vp_sig.imag, label=filename, fmt="."
            )
        else:
            ax[0].plot(1e-3 * vp_rho, vp_vis.real, "o", markersize=2.8, label=filename)
            ax[1].plot(1e-3 * vp_rho, vp_vis.imag, "o", markersize=2.8, label=filename)

    allmaxvis = np.max(maxvis)
    allminvis = np.min(minvis)
    if (allminvis < 0) or (allminvis - 0.1 * allmaxvis < 0):
        ax[0].axis([0, np.max(uvbins), allminvis - 0.1 * allmaxvis, 1.1 * allmaxvis])
        ax[1].axis([0, np.max(uvbins), allminvis - 0.1 * allmaxvis, 1.1 * allmaxvis])
    else:
        ax[0].axis([0, np.max(uvbins), 0.0, 1.1 * allmaxvis])
        ax[1].axis([0, np.max(uvbins), 0.0, 1.1 * allmaxvis])

    ax[0].plot([0, np.max(uvbins)], [0, 0], "--k")
    ax[1].plot([0, np.max(uvbins)], [0, 0], "--k")
    plt.xlabel("deprojected baseline length [kilo$\lambda$]")
    ax[0].set_ylabel("average real [Jy]")
    ax[1].set_ylabel("average imag [Jy]")
    ax[0].legend()
    plt.savefig("deprojected_PA{:.1f}_i{:.1f}.png".format(PA, incl))
    plt.show(block=False)


def estimate_flux_scale(
    reference, comparison, incl=0, PA=0, uvbins=None, offx=0, offy=0
):
    """
    Calculates the weighted average of the flux ratio between two observations
    of a source. The minimum baseline compared is the longer of the minimum
    baselines in the individual datasets. The longest baseline compared is
    either the shorter of the longest baselines in the individual datasets,
    or 800 kilolambda.

    Args:
        reference (str): name of .npz file holding the reference dataset
            (with the "correct" flux")
        comparison (str): name of .npz file holding the comparison dataset
            (with the flux ratio being checked)
        filelist (list): names of .npz files storing visibility data
        incl (float): Inclination of disk [degrees]
        PA (float): Position angle of disk [degrees]
        offx (float): Horizontal offset of disk center from phase center
            [arcseconds]
        offy (float): Vertical offset of disk center from phase center
            [arcseconds]
        uvbins (array): bins at which to compare the visibility values [lambda]
            By default, the minimum baseline compared is the longer of the
            minimum baselines in the individual datasets. The longest
            baseline compared is either the shorter of the longest baselines
            in the individual datasets, or 800 kilolambda, whichever
            comes first.
    """

    inpf = np.load(reference)
    u_ref = inpf["u"]
    v_ref = inpf["v"]
    vis_ref = inpf["Vis"]
    wgt_ref = inpf["Wgt"]

    inpf = np.load(comparison)
    u_comp = inpf["u"]
    v_comp = inpf["v"]
    vis_comp = inpf["Vis"]
    wgt_comp = inpf["Wgt"]

    uvdist_ref = np.sqrt(u_ref ** 2 + v_ref ** 2)
    uvdist_comp = np.sqrt(u_comp ** 2 + v_comp ** 2)

    mindist = np.max(np.array([np.min(uvdist_ref), np.min(uvdist_comp)]))
    # the maximum baseline we want to compare is the longest shared baseline or
    # 800 kilolambda, whichever comes first (we don't want to go out to a
    # baseline that's too long because phase decorrelation becomes a bigger
    # issue at longer baselines.
    maxdist = np.min(np.array([np.max(uvdist_ref), np.max(uvdist_ref), 8e5]))
    if uvbins is None:
        uvbins = mindist / 1.0e3 + 10.0 * np.arange(
            np.floor((maxdist - mindist) / 1.0e4)
        )

    # deproject the visibilities and do the annular averaging
    vp = deproject_vis(
        [u_ref, v_ref, vis_ref, wgt_ref],
        bins=uvbins,
        incl=incl,
        PA=PA,
        offx=offx,
        offy=offy,
    )
    ref_rho, ref_vis, ref_sig = vp

    # deproject the visibilities and do the annular averaging
    vp = deproject_vis(
        [u_comp, v_comp, vis_comp, wgt_comp],
        bins=uvbins,
        incl=incl,
        PA=PA,
        offx=offx,
        offy=offy,
    )
    comp_rho, comp_vis, comp_sig = vp

    maxlen = np.min(np.array([len(comp_rho), len(ref_rho)]))

    # we only want to compare overlapping baseline intervals
    rho_intersection = np.intersect1d(ref_rho, comp_rho)

    # they're the same for the real and imaginary components
    comp_sig_intersection = comp_sig[np.where(np.in1d(comp_rho, rho_intersection))].real
    comp_vis_intersection = comp_vis[np.where(np.in1d(comp_rho, rho_intersection))]
    ref_sig_intersection = ref_sig[np.where(np.in1d(ref_rho, rho_intersection))].real
    ref_vis_intersection = ref_vis[np.where(np.in1d(ref_rho, rho_intersection))]

    ratio = np.abs(comp_vis_intersection) / np.abs(ref_vis_intersection)
    err = ratio * np.sqrt(
        (comp_sig_intersection / np.abs(comp_vis_intersection)) ** 2
        + (ref_sig_intersection / np.abs(ref_vis_intersection)) ** 2
    )

    w = 1 / err ** 2
    ratio_avg = np.sum(w * ratio) / np.sum(w)
    print(
        "#The ratio of the fluxes of %s to %s is %.5f"
        % (comparison, reference, ratio_avg)
    )
    print(
        "#The scaling factor for gencal is %.3f for your comparison \
        measurement"
        % (sqrt(ratio_avg))
    )
    print(
        "#The error on the weighted mean ratio is %.3e, although it's \
        likely that the weights in the measurement sets are off by some \
        constant factor"
        % (1 / np.sqrt(np.sum(w)),)
    )
    plt.figure()
    plt.errorbar(
        1e-3 * rho_intersection, ratio, yerr=err, fmt=".", label="Binned ratios"
    )
    plt.plot(
        1e-3 * rho_intersection,
        np.ones_like(ratio) * ratio_avg,
        label="weighted average",
    )
    plt.ylabel("Visibility amplitude ratios")
    plt.xlabel("UV distance (kilolambda)")
    plt.legend()
    plt.savefig("ratio_PA{:.1f}_i{:.1f}.png".format(PA, incl))
    plt.show(block=False)


def estimate_SNR(imagename, disk_mask, noise_mask, chans=None):
    """
    Estimate peak SNR of source

    Args:
        imagename (string): Image name ending in '.image'
        disk_mask: in the CASA region format
        noise_mask (string): Annulus to measure image rms, in the CASA region
            format, e.g. 'annulus[[500pix, 500pix],["1arcsec", "2arcsec"]]'
    """
    headerlist = casatasks.imhead(imagename, mode="list")
    beammajor = headerlist["beammajor"]["value"]
    beamminor = headerlist["beamminor"]["value"]
    beampa = headerlist["beampa"]["value"]
    print("# %s" % imagename)
    print(
        "# Beam %.3f arcsec x %.3f arcsec (%.2f deg)" % (beammajor, beamminor, beampa)
    )
    if chans is not None:
        disk_stats = casatasks.imstat(
            imagename=imagename, region=disk_mask, chans=chans
        )
    else:
        disk_stats = casatasks.imstat(imagename=imagename, region=disk_mask)
    disk_flux = disk_stats["flux"][0]
    print("# Flux inside disk mask: %.3f mJy" % (disk_flux * 1000,))
    peak_intensity = disk_stats["max"][0]
    print("# Peak intensity of source: %.3f mJy/beam" % (peak_intensity * 1000,))
    # print("npoints", disk_stats["npts"])
    # print("sum", disk_stats["sum"])
    # print("mean", disk_stats["mean"])
    # print("sigma", disk_stats["sigma"])
    if chans is not None:
        rms = casatasks.imstat(imagename=imagename, region=noise_mask, chans=chans)[
            "rms"
        ][0]
    else:
        rms = casatasks.imstat(imagename=imagename, region=noise_mask)["rms"][0]
    print("# rms: %.3e mJy/beam" % (rms * 1000,))
    SNR = peak_intensity / rms
    print("# Peak SNR: %.3f" % (SNR,))
    print()


def get_station_numbers(msfile, antenna_name):
    """
    Get the station numbers for all observations in which the given
    antenna appears.

    Args:
        msfile (string): name of measurement set
        antenna_name (string): name of antenna (e.g. "DA48")
    """
    tb.open(msfile + "/ANTENNA")
    ant_names = tb.getcol("NAME")
    ant_stations = tb.getcol("STATION")
    tb.close()

    ant_numbers = np.where(ant_names == antenna_name)[0]

    tb.open(msfile)
    antenna1 = tb.getcol("ANTENNA1")
    obsid = tb.getcol("OBSERVATION_ID")
    tb.close()

    for i in ant_numbers:
        matching_obs = np.unique(obsid[np.where(antenna1 == i)])
        for j in matching_obs:
            print("#Observation ID %d: %s@%s" % (j, antenna_name, ant_stations[i]))
