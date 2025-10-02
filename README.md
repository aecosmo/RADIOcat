The files contain a set of codes for making a mosaic from radio images. Here, I have used uGMRT band 2 images centered on a single frequency of 147.4 MHz. 
The following packages are required: Astropy, pybdsf, CASA, and aegean (for completeness correction only). 
You can use these codes to make a mosaic, build a source catalogue, classify point-like and extend sources, and compute source counts with necessary corrections such as false detection rate, completeness, and visibility area. 
Necessary images can be requested from the developer (me) and the collaborators. 
The TGSS, GLEAM, and GLEAM-X are available in the Vizier Astronomical database. The uGMRT catalogue is present in the directory, and will soon be moved to Vizier. 


To use these codes you can follow these instructions:
    
    Start with the *.SP2B.PBCOR.FITS images
    PSF matching: casa -c imsmooth.py
    mosaicing.ipynb
        Reads PSF-matched images, runs pybdsf on these to generate the rms maps (which are used as weights in the next step) 
        Makes the mosaic using PSF-matched images and their weights, runs pybdsf on the mosaic to generate the rms map of the mosaic
        Calculates the median rms of the image. The number goes to the abstract of the paper. 
    
    plots.ipynb
        Make the plots: obs strategy, multiplePCs, mosaic, triple plot, rms map, completeness, 
    
    Pointlike or extended. pointextended.ipynb
    
    crossmatch.ipynb astrometry and flux comparison. NVSS-uGMRT cross-matched spectral index. 
    
    catalogue.ipynb to print a sample catalogue
    
    crossmatch-bright-spind.ipynb finds the bright sources and extracts their information from different catalogues. 
    
    image_summary.inynb gives a summary of the mosaic


    
