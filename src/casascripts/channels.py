import table

from .constants import cc


def vel_to_chan(msfile, field, obsid, spw, restfreq, vel):
    """
    Identifies the channel(s) corresponding to input LSRK velocities. 
    Useful for choosing which channels to split out or flag if a line is 
    expected to be present

    Args:
        msfile (string): name of measurement set 
        field (string): field name
        spw (int): Spectral window number 
        obsid (int): Observation ID corresponding to the selected 
            spectral window 
        restfreq (float): Rest frequency [Hz]
        vel (float or array of floats): input velocity in LSRK frame 
            km/s]
            

    Returns:
        (array) or (int) channel number most closely corresponding to 
            input LSRK velocity 
    """

    table.open(msfile + "/SPECTRAL_WINDOW")
    chanfreqs = table.getcol("CHAN_FREQ", startrow=spw, nrow=1)
    table.close()
    table.open(msfile + "/FIELD")
    fieldnames = table.getcol("NAME")
    table.close()
    table.open(msfile + "/OBSERVATION")
    obstime = np.squeeze(table.getcol("TIME_RANGE", startrow=obsid, nrow=1))[0]
    table.close()
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
    # convert to LSRK velocities [km/s]
    chanvelocities = (restfreq - lsrkfreqs) / restfreq * cc_kms
    ms.close()
    if type(vel) == np.ndarray:
        outchans = np.zeros_like(vel)
        for i in range(len(vel)):
            outchans[i] = np.argmin(np.abs(chanvelocities - vel[i]))
        return outchans
    else:
        return np.argmin(np.abs(chanvelocities - vel))


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
        chans = vel_to_chan(
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

