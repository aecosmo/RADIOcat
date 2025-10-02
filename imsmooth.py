input_images = ["AAAA.SP2B.PBCOR.FITS", "BBBB.SP2B.PBCOR.FITS", "CCCC.SP2B.PBCOR.FITS", "DDDD.SP2B.PBCOR.FITS", 
                "EEEE.SP2B.PBCOR.FITS", "FFFF.SP2B.PBCOR.FITS", "GGGG.SP2B.PBCOR.FITS"]
                

output_images = ["AAAA.SP2B.PBCOR.SMOOTH.FITS", "BBBB.SP2B.PBCOR.SMOOTH.FITS", "CCCC.SP2B.PBCOR.SMOOTH.FITS", "DDDD.SP2B.PBCOR.SMOOTH.FITS", 
                "EEEE.SP2B.PBCOR.SMOOTH.FITS", "FFFF.SP2B.PBCOR.SMOOTH.FITS", "GGGG.SP2B.PBCOR.SMOOTH.FITS"]
                
for ii in range(len(input_images)):
    # From within CASA
    importfits(fitsimage=input_images[ii], imagename=input_images[ii]+'.image', overwrite=True)
    imsmooth(imagename=input_images[ii]+'.image', outfile=output_images[ii]+'.image', kernel='gauss', targetres=True, major='40arcsec', minor='22arcsec', pa="23deg", overwrite=True)
    exportfits(imagename=output_images[ii]+'.image', fitsimage=output_images[ii], overwrite=True)
    
import os
os.system('rm -rf *.image')
