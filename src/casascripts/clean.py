def tclean_wrapper(
    vis,
    imagename,
    scales,
    smallscalebias=0.6,
    mask="",
    threshold="0.2mJy",
    imsize=None,
    cellsize=None,
    interactive=False,
    robust=0.5,
    gain=0.3,
    niter=50000,
    cycleniter=300,
    uvtaper=[],
    savemodel="none",
    pbcor=False,
):
    """
    Wrapper for tclean with keywords set to values desired for the Large 
    Program imaging.  See the CASA 5.1.1 documentation for tclean to get the 
    definitions of all the parameters.
    """
    if imsize is None:
        if "LB" in vis or "combined" in vis:
            imsize = 3000
        elif "SB" in vis:
            imsize = 900
        else:
            print("Error: need to set imsize manually")

    if cellsize is None:
        if "LB" in vis or "combined" in vis:
            cellsize = ".003arcsec"
        elif "SB" in vis:
            cellsize = ".03arcsec"
        else:
            print("Error: need to set cellsize manually")

    for ext in [
        ".image",
        ".mask",
        ".model",
        ".pb",
        ".psf",
        ".residual",
        ".sumwt",
        ".image.pbcor",
    ]:
        os.system("rm -rf " + imagename + ext)
    tclean(
        vis=vis,
        imagename=imagename,
        specmode="mfs",
        deconvolver="multiscale",
        scales=scales,
        weighting="briggs",
        robust=robust,
        gain=gain,
        imsize=imsize,
        cell=cellsize,
        # set smallscalebias to CASA's default of 0.6 unless manually changed
        smallscalebias=smallscalebias,
        niter=niter,  # we want to end on the threshold
        interactive=interactive,
        threshold=threshold,
        cycleniter=cycleniter,
        cyclefactor=1,
        uvtaper=uvtaper,
        mask=mask,
        savemodel=savemodel,
        nterms=1,
        pbcor=pbcor,
    )

    # this step is a workaround a bug in tclean that doesn't always save the
    # model during multiscale clean. See the "Known Issues" section for
    # CASA 5.1.1 on NRAO's website
    # since we're using CASA 5.4, it seems like it's working
    # actually I take that back
    if savemodel == "modelcolumn":
        print()
        print("Running tclean a second time to save the model...")
        tclean(
            vis=vis,
            imagename=imagename,
            specmode="mfs",
            deconvolver="multiscale",
            scales=scales,
            weighting="briggs",
            robust=robust,
            gain=gain,
            imsize=imsize,
            cell=cellsize,
            smallscalebias=smallscalebias,
            niter=0,
            interactive=False,
            threshold=threshold,
            cycleniter=cycleniter,
            cyclefactor=1,
            uvtaper=uvtaper,
            mask="",
            savemodel=savemodel,
            calcres=False,
            calcpsf=False,
            nterms=1,
            pbcor=pbcor,
        )


def image_each_obs(
    ms_dict,
    prefix,
    scales,
    smallscalebias=0.6,
    mask="",
    threshold="0.2mJy",
    imsize=None,
    cellsize=None,
    interactive=False,
    robust=0.5,
    gain=0.3,
    niter=50000,
    cycleniter=300,
):
    """
    Wrapper for tclean that will loop through all the observations in a 
    measurement set and image them individually.

    Args:
        ms_dict (dict): information about measurement set
        prefix (string): prefix for all output file names 
    
    See the CASA 5.1.1 documentation for tclean to get the definitions of all 
    other parameters.
    """
    msfile = prefix + "_" + ms_dict["name"] + "_initcont.ms"
    tb.open(msfile + "/OBSERVATION")
    # picked an arbitrary column to count the number of observations
    num_observations = (tb.getcol("TIME_RANGE")).shape[1]
    tb.close()

    if imsize is None:
        if ms_dict["name"] == "LB1":
            imsize = 3000
        else:
            imsize = 900

    if cellsize is None:
        if ms_dict["name"] == "LB1":
            cellsize = ".003arcsec"
        else:
            imsize = 900
            cellsize = ".03arcsec"

    # start of CASA commands
    for i in range(num_observations):
        observation = "%d" % i
        imagename = prefix + "_" + ms_dict["name"] + "_initcont_exec%s" % observation
        for ext in [".image", ".mask", ".model", ".pb", ".psf", ".residual", ".sumwt"]:
            os.system("rm -rf " + imagename + ext)
        tclean(
            vis=msfile,
            imagename=imagename,
            observation=observation,
            specmode="mfs",
            deconvolver="multiscale",
            scales=scales,
            weighting="briggs",
            robust=robust,
            gain=gain,
            imsize=imsize,
            cell=cellsize,
            smallscalebias=smallscalebias,
            niter=niter,
            interactive=interactive,
            threshold=threshold,
            cycleniter=cycleniter,
            cyclefactor=1,
            mask=mask,
            nterms=1,
        )

    # to delete the model, use
    # delmod(vis=msfile, otf=True, scr=True)

    print(
        "Each observation saved in the format %sOBSERVATIONNUMBER.image"
        % (prefix + "_" + ms_dict["name"] + "_initcont_exec",)
    )

