def split_all_obs(msfile, nametemplate):
    """
    Split out individual observations in a measurement set 

    Args:
        msfile (string): name of measurement set, ending in '.ms'   
        nametemplate (string): template name of output measurement sets for 
            individual observations 
    """
    tb.open(msfile)
    spw_col = tb.getcol("DATA_DESC_ID")
    obs_col = tb.getcol("OBSERVATION_ID")
    field_col = tb.getcol("FIELD_ID")
    tb.close()

    obs_ids = np.unique(obs_col)

    # yes, it would be more logical to split out by observation id,
    # but splitting out by observation id in practice leads to some issues with
    # the metadata
    for i in obs_ids:
        spws = np.unique(spw_col[np.where(obs_col == i)])
        # sometimes the MS secretly has multiple field IDs lurking even if
        # listobs only shows one field
        fields = np.unique(field_col[np.where(obs_col == i)])
        if len(spws) == 1:
            spw = str(spws[0])
        else:
            spw = "%d~%d" % (spws[0], spws[-1])

        if len(fields) == 1:
            field = str(fields[0])
        else:
            field = "%d~%d" % (fields[0], fields[-1])
        # start of CASA commands
        outputvis = nametemplate + "%d.ms" % i
        os.system("rm -rf " + outputvis)
        print("#Saving observation %d of %s to %s" % (i, msfile, outputvis))
        split(vis=msfile, spw=spw, field=field, outputvis=outputvis, datacolumn="data")


def export_MS(msfile):
    """
    Spectrally averages visibilities to a single channel per SPW and exports 
        to .npz file.

    Args:
        msfile (string): Name of CASA measurement set, ending in '.ms' 

    Returns: 
        None
    """
    filename = msfile
    if filename[-3:] != ".ms":
        print("MS name must end in '.ms'")
        return
    # strip off the '.ms'
    MS_filename = filename.replace(".ms", "")

    # get information about spectral windows

    tb.open(MS_filename + ".ms/SPECTRAL_WINDOW")
    num_chan = tb.getcol("NUM_CHAN").tolist()
    tb.close()

    # spectral averaging (1 channel per SPW)

    os.system("rm -rf %s" % MS_filename + "_spavg.ms")
    split(
        vis=MS_filename + ".ms",
        width=num_chan,
        datacolumn="data",
        outputvis=MS_filename + "_spavg.ms",
    )

    # get the data tables
    tb.open(MS_filename + "_spavg.ms")
    data = np.squeeze(tb.getcol("DATA"))
    flag = np.squeeze(tb.getcol("FLAG"))
    uvw = tb.getcol("UVW")
    weight = tb.getcol("WEIGHT")
    spwid = tb.getcol("DATA_DESC_ID")
    tb.close()

    # get frequency information
    tb.open(MS_filename + "_spavg.ms/SPECTRAL_WINDOW")
    freqlist = np.squeeze(tb.getcol("CHAN_FREQ"))
    tb.close()

    # get rid of any flagged columns
    good = np.squeeze(np.any(flag, axis=0) == False)
    data = data[:, good]
    weight = weight[:, good]
    uvw = uvw[:, good]
    spwid = spwid[good]

    # compute spatial frequencies in lambda units
    get_freq = lambda ispw: freqlist[ispw]
    # get spectral frequency corresponding to each datapoint
    freqs = get_freq(spwid)
    u = uvw[0, :] * freqs / 2.9979e8
    v = uvw[1, :] * freqs / 2.9979e8

    # average the polarizations
    Re = np.sum(data.real * weight, axis=0) / np.sum(weight, axis=0)
    Im = np.sum(data.imag * weight, axis=0) / np.sum(weight, axis=0)
    Vis = Re + 1j * Im
    Wgt = np.sum(weight, axis=0)

    # output to npz file and delete intermediate measurement set
    os.system("rm -rf %s" % MS_filename + "_spavg.ms")
    os.system("rm -rf " + MS_filename + ".vis.npz")
    np.savez(MS_filename + ".vis", u=u, v=v, Vis=Vis, Wgt=Wgt)
    print("#Measurement set exported to %s" % (MS_filename + ".vis.npz",))


def deproject_vis(
    data, bins=np.array([0.0]), incl=0.0, PA=0.0, offx=0.0, offy=0.0, errtype="mean"
):
    """
    Deproject and azimuthally average visibilities 

    Args:
    data (4-tuple): u,v, visibilities, and weight arrays 
    bins (1-D array):  uv distance bins [kilolambda]
    incl (float): Inclination of disk [degrees]
    PA (float): Position angle of disk [degrees]
    offx (float): Horizontal offset of disk center from phase center 
        [arcseconds]
    offy (float): Vertical offset of disk center from phase center [arcseconds]

    Returns:
        (3-tuple) uv distance bins (1D array), visibilities (1D array), 
            errors on averaged visibilities (1D array) 
    """

    # - read in, parse data
    u, v, vis, wgt = data
    # - convert keywords into relevant units
    inclr = np.radians(incl)
    PAr = 0.5 * np.pi - np.radians(PA)
    offx *= -np.pi / (180.0 * 3600.0)
    offy *= -np.pi / (180.0 * 3600.0)

    # - change to a deprojected, rotated coordinate system
    uprime = u * np.cos(PAr) + v * np.sin(PAr)
    vprime = (-u * np.sin(PAr) + v * np.cos(PAr)) * np.cos(inclr)
    rhop = np.sqrt(uprime ** 2 + vprime ** 2)

    # - phase shifts to account for offsets
    shifts = np.exp(-2.0 * np.pi * 1.0j * (u * -offx + v * -offy))
    visp = vis * shifts
    realp = visp.real
    imagp = visp.imag

    # - if requested, return a binned (averaged) representation
    if bins.size > 1.0:
        avbins = 1e3 * bins  # scale to lambda units (input in klambda)
        bwid = 0.5 * (avbins[1] - avbins[0])
        bvis = np.zeros_like(avbins, dtype="complex")
        berr = np.zeros_like(avbins, dtype="complex")
        for ib in np.arange(len(avbins)):
            inb = np.where((rhop >= avbins[ib] - bwid) & (rhop < avbins[ib] + bwid))
            if len(inb[0]) >= 5:
                bRe, eRemu = np.average(realp[inb], weights=wgt[inb], returned=True)
                eRese = np.std(realp[inb])
                bIm, eImmu = np.average(imagp[inb], weights=wgt[inb], returned=True)
                eImse = np.std(imagp[inb])
                bvis[ib] = bRe + 1j * bIm
                if errtype == "scat":
                    berr[ib] = eRese + 1j * eImse
                else:
                    berr[ib] = 1.0 / np.sqrt(eRemu) + 1j / np.sqrt(eImmu)
            else:
                bvis[ib] = 0 + 1j * 0
                berr[ib] = 0 + 1j * 0
        parser = np.where(berr.real != 0)
        output = avbins[parser], bvis[parser], berr[parser]
        return output

    # - if not, returned the unbinned representation
    output = rhop, realp + 1j * imagp, 1.0 / np.sqrt(wgt)

    return output


def rescale_flux(vis, gencalparameter):
    """
    Rescale visibility fluxes using gencal, then split into a new 
    measurement set.
 
    Args:
        vis (string): Measurement set name, ending in ms 
        gencalparameter (list): flux rescaling parameters to be passed to 
            `parameter` for gencal task
    """
    caltable = "scale_" + vis.replace(".ms", ".gencal")
    os.system("rm -rf " + caltable)
    gencal(vis=vis, caltable=caltable, caltype="amp", parameter=gencalparameter)
    applycal(vis=vis, gaintable=caltable, calwt=True, flagbackup=True)
    vis_rescaled = vis.replace(".ms", "_rescaled.ms")
    print("#Splitting out rescaled values into new MS: %s" % (vis_rescaled,))
    os.system("rm -rf " + vis_rescaled + "*")
    split(vis=vis, outputvis=vis_rescaled, datacolumn="corrected")
