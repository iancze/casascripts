def avg_cont(
    ms_dict,
    output_prefix,
    flagchannels="",
    maxchanwidth=125,
    datacolumn="data",
    contspws=None,
    width_array=None,
):
    """
    Produce spectrally averaged continuum measurement sets 

    Args:
        ms_dict: Dictionary of information about measurement set
        output_prefix (string): Prefix for all output file names 
        flagchannels: Argument to be passed for flagchannels parameter in 
            flagdata task
        maxchanwidth: Maximum width of channel (MHz). This is the value 
            recommended by ALMA for Band 6 to avoid bandwidth smearing
        datacolumn: Column to pull from for continuum averaging (usually will 
            be 'data', but may sometimes be 'corrected' if there was flux 
            rescaling applied)
        contspws (string): Argument to be passed to CASA for the spw parameter 
            in split. If not set, all SPWs will be selected by default. 
        width_array (array): Argument to be passed to CASA for the width 
            parameter in split. If not set, all SPWs will be selected by default. 

    Returns:
        None
    """
    msfile = ms_dict["vis"]
    tb.open(msfile + "/SPECTRAL_WINDOW")
    total_bw = tb.getcol("TOTAL_BANDWIDTH")
    num_chan = tb.getcol("NUM_CHAN")
    tb.close()
    if width_array is None and contspws is None:
        # array of number of channels to average to form an output channel
        # (to be passed to mstransform)
        width_array = (
            (num_chan / np.ceil(total_bw / (1.0e6 * maxchanwidth)))
            .astype("int")
            .tolist()
        )
        contspws = "%d~%d" % (0, len(total_bw) - 1)  # by default select all SPWs
    elif (width_array is not None and contspws is None) or (
        width_array is None and contspws is not None
    ):
        raise ValueError(
            "If either contspws or width_array is set to a value,\
            the other parameter has to be manually set as well"
        )

    if ms_dict["name"] == "LB1":
        timebin = "6s"
    else:
        timebin = "0s"  # default in CASA

    # start of CASA commands
    if len(flagchannels) == 0:
        outputvis = output_prefix + "_" + ms_dict["name"] + "_initcont.ms"
        os.system("rm -rf " + outputvis)
        split(
            vis=msfile,
            field=ms_dict["field"],
            spw=contspws,
            outputvis=outputvis,
            width=width_array,
            timebin=timebin,
            datacolumn=datacolumn,
            intent="OBSERVE_TARGET#ON_SOURCE",
            keepflags=False,
        )
    else:
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
            field=ms_dict["field"],
        )

        outputvis = output_prefix + "_" + ms_dict["name"] + "_initcont.ms"
        os.system("rm -rf " + outputvis)
        split(
            vis=msfile,
            field=ms_dict["field"],
            spw=contspws,
            outputvis=outputvis,
            width=width_array,
            timebin=timebin,
            datacolumn=datacolumn,
            intent="OBSERVE_TARGET#ON_SOURCE",
            keepflags=False,
        )

        # restore flagged spectral line channels
        flagmanager(vis=msfile, mode="restore", versionname="before_cont_flags")

    print("#Averaged continuum dataset saved to %s" % outputvis)

