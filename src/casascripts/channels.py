from casatools import table, ms
import numpy as np

from .constants import cc_kms


def vel_to_chan(msfile, field, obsid, spw, restfreq, vel):
    """
    Identifies the channel(s) corresponding to input LSRK velocities.
    Useful for choosing which channels to split out or flag if a line is
    expected to be present.

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

    tb = table()
    mstool = ms()
    # open the file
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
    mstool.open(msfile)
    lsrkfreqs = mstool.cvelfreqs(
        spwids=[spw],
        mode="channel",
        nchan=nchan,
        obstime=str(obstime) + "s",
        start=0,
        outframe="LSRK",
    )
    # convert to LSRK velocities [km/s]
    chanvelocities = (restfreq - lsrkfreqs) / restfreq * cc_kms
    mstool.close()
    if type(vel) == np.ndarray:
        outchans = np.zeros_like(vel)
        for i in range(len(vel)):
            outchans[i] = np.argmin(np.abs(chanvelocities - vel[i]))
        return outchans
    else:
        return np.argmin(np.abs(chanvelocities - vel))


def get_flagchannels(
    msfile, field, obsids, spws, restfreqs, velocity_range=np.array([-20, 20])
):
    """
    Identify channels to flag based on provided velocity range of the
        line emission

    Args:
        ms_dict (dictionary): Dictionary of information about measurement set
        velocity_range np.array([min_velocity, max_velocity]): Velocity range
            (in km/s) over which line emission has been identified

    Returns:
        String of channels to be flagged, in a format that can be passed to
            the spw parameter in CASA's flagdata task.
    """
    flagchannels_string = ""
    for j, (obsid, spw, restfreq) in enumerate(zip(obsids, spws, restfreqs)):
        chans = vel_to_chan(
            msfile,
            field,
            obsid,
            spw,
            restfreq,
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
    return flagchannels_string
