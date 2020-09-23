from casatasks import tclean

import os
import shutil


def tclean_wrapper(
    vis,
    imagename,
    imsize=None,
    cellsize=None,
    interactive=False,
    robust=0.5,
    niter=50000,
    pbcor=True,
    **kwargs
):

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
        fname = imagename + ext
        if os.path.exists(fname):
            shutil.rmtree(fname)

    tclean(
        vis=vis,
        imagename=imagename,
        weighting="briggs",
        robust=robust,
        imsize=imsize,
        cell=cellsize,
        niter=niter,  # we want to end on the threshold
        interactive=interactive,
        pbcor=pbcor,
        **kwargs
    )

    # # this step is a workaround a bug in tclean that doesn't always save the
    # # model during multiscale clean. See the "Known Issues" section for
    # # CASA 5.1.1 on NRAO's website
    # # since we're using CASA 5.4, it seems like it's working
    # # actually I take that back
    # if savemodel == "modelcolumn":
    #     print()
    #     print("Running tclean a second time to save the model...")
    #     tclean(
    #         vis=vis,
    #         imagename=imagename,
    #         specmode="mfs",
    #         deconvolver="multiscale",
    #         scales=scales,
    #         weighting="briggs",
    #         robust=robust,
    #         gain=gain,
    #         imsize=imsize,
    #         cell=cellsize,
    #         smallscalebias=smallscalebias,
    #         niter=0,
    #         interactive=False,
    #         threshold=threshold,
    #         cycleniter=cycleniter,
    #         cyclefactor=1,
    #         uvtaper=uvtaper,
    #         mask="",
    #         savemodel=savemodel,
    #         calcres=False,
    #         calcpsf=False,
    #         nterms=1,
    #         pbcor=pbcor,
    #     )