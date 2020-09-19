def LSRKvel_to_chan(msfile, field, obsid, spw, restfreq, LSRKvelocity):
    """
    Identifies the channel(s) corresponding to input LSRK velocities. 
    Useful for choosing which channels to split out or flag if a line is 
    expected to be present

    Args:
        msfile (string): name of measurement set 
        field (str): field name
        spw (int): Spectral window number 
        obsid (int): Observation ID corresponding to the selected 
            spectral window 
        restfreq (float): Rest frequency [Hz]
        LSRKvelocity (float or array of floats): input velocity in LSRK frame 
            km/s]
            

    Returns:
        (array) or (int) channel number most closely corresponding to 
            input LSRK velocity 
    """
    cc = 299792458.0  # speed of light in m/s

    print("In LSRKvel", field, spw, restfreq, LSRKvelocity)

    # IC: the problem we have here is that for the second execution,
    # "DATA_DESC_ID" does not map to the same spectral windows that `listobs`
    # tells us about it seems like these things are only used to find obsid,
    # so maybe we can find it a different way
    # just specifying it works fine

    # tb.open(msfile)
    # spw_col = tb.getcol('DATA_DESC_ID')
    # obs_col = tb.getcol('OBSERVATION_ID')
    # print("spw_col", spw_col)
    # print("spw_col_unique", np.unique(spw_col))
    # print("obs_col", obs_col)

    # tb.close()
    # obsid = np.unique(obs_col[np.where(spw_col==spw)])
    # print("obsid", obsid)

    tb.open(msfile + "/SPECTRAL_WINDOW")
    chanfreqs = tb.getcol("CHAN_FREQ", startrow=spw, nrow=1)
    tb.close()
    tb.open(msfile + "/FIELD")
    fieldnames = tb.getcol("NAME")
    tb.close()
    tb.open(msfile + "/OBSERVATION")
    obstime = np.squeeze(tb.getcol("TIME_RANGE", startrow=obsid, nrow=1))[0]
    tb.close()
    nchan = len(chanfreqs)
    ms.open(msfile)
    lsrkfreqs = ms.cvelfreqs(
        spwids=[spw],
        fieldids=np.where(fieldnames == field)[0][0],
        mode="channel",
        nchan=nchan,
        obstime=str(obstime) + "s",
        start=0,
        outframe="LSRK",
    )
    # converted to LSRK velocities in km/s
    chanvelocities = (restfreq - lsrkfreqs) / restfreq * cc / 1.0e3
    ms.close()
    if type(LSRKvelocity) == np.ndarray:
        outchans = np.zeros_like(LSRKvelocity)
        for i in range(len(LSRKvelocity)):
            outchans[i] = np.argmin(np.abs(chanvelocities - LSRKvelocity[i]))
        return outchans
    else:
        return np.argmin(np.abs(chanvelocities - LSRKvelocity))


def get_flagchannels(ms_dict, output_prefix, velocity_range=np.array([-20, 20])):
    """
    Identify channels to flag based on provided velocity range of the 
        line emission

    Args:
        ms_dict (dictionary): Dictionary of information about measurement set
        output_prefix (string): Prefix for all output file names 
        velocity_range np.array([min_velocity, max_velocity]): Velocity range 
            (in km/s) over which line emission has been identified

    Returns:
        String of channels to be flagged, in a format that can be passed to 
            the spw parameter in CASA's flagdata task. 
    """
    flagchannels_string = ""
    for j, spw in enumerate(ms_dict["line_spws"]):
        chans = LSRKvel_to_chan(
            ms_dict["vis"],
            ms_dict["field"],
            ms_dict["obsids"][j],
            spw,
            ms_dict["line_freqs"][j],
            velocity_range,
        )
        if j == 0:
            flagchannels_string += "%d:%d~%d" % (
                spw,
                np.min([chans[0], chans[1]]),
                np.max([chans[0], chans[1]]),
            )
        else:
            flagchannels_string += ", %d:%d~%d" % (
                spw,
                np.min([chans[0], chans[1]]),
                np.max([chans[0], chans[1]]),
            )
    print(
        "# Flagchannels input string for %s: '%s'"
        % (ms_dict["name"], flagchannels_string)
    )

    return flagchannels_string

