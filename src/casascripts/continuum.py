from casatools import table
from casatasks import flagmanager, flagdata, split
import os


def avg_cont(
    msfile,
    outputvis,
    flagchannels="",
    maxchanwidth=125,
    datacolumn="data",
    spws=None,
    width_array=None,
):
    """
    Produce spectrally averaged continuum measurement sets

    Args:
        msfile (string): input file name
        outputvis (string): output file name
        flagchannels: Argument to be passed for flagchannels parameter in
            flagdata task
        maxchanwidth: Maximum width of channel (MHz). This is the value
            recommended by ALMA for Band 6 to avoid bandwidth smearing
        datacolumn: Column to pull from for continuum averaging (usually will
            be 'data', but may sometimes be 'corrected' if there was flux
            rescaling applied)
        width_array (array): Argument to be passed to CASA for the width
            parameter in split. If not set, all SPWs will be selected by default.

    Returns:
        None
    """
    tb = table()

    tb.open(msfile + "/SPECTRAL_WINDOW")
    total_bw = tb.getcol("TOTAL_BANDWIDTH")
    num_chan = tb.getcol("NUM_CHAN")
    tb.close()

    timebin = "0s"  # default in CASA

    if os.path.isdir(msfile + ".flagversions/flags.before_cont_flags"):
        # clear out old versions of the flags
        flagmanager(vis=msfile, mode="delete", versionname="before_cont_flags")

    # save flag state before flagging spectral lines
    flagmanager(
        vis=msfile,
        mode="save",
        versionname="before_cont_flags",
        comment="Flag states before spectral lines are flagged",
    )
    # flag spectral lines
    flagdata(
        vis=msfile,
        mode="manual",
        spw=flagchannels,
        flagbackup=False,
    )

    os.system("rm -rf " + outputvis)
    split(
        vis=msfile,
        outputvis=outputvis,
        spw=spws,
        width=width_array,
        timebin=timebin,
        datacolumn=datacolumn,
        intent="OBSERVE_TARGET#ON_SOURCE",
        keepflags=False,
    )

    # restore flagged spectral line channels
    flagmanager(vis=msfile, mode="restore", versionname="before_cont_flags")

    print("#Averaged continuum dataset saved to %s" % outputvis)
