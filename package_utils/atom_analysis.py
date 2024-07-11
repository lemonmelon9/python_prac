# Drift correction of STM images obtained on monoclinic Ta2NiSe5.
def Ta2NiSe5_driftcorr(z, a_image_nm, c_image_nm, scansize_nm, origin_line=0, origin_pixel=0):
    import numpy as np
    lines, pixels = np.shape(z)
    a, c = 0.34916, 1.565 # lattice constants (nm)
    lines_corr = int(round(lines * (a_image_nm / a)))
    pixels_corr = int(round(pixels * (c_image_nm / c)))
    # print (lines_corr, pixels_corr)
    z_corr = np.copy(z)[origin_line:origin_line+lines_corr, origin_pixel:origin_pixel+pixels_corr]
    
    # length of scansize must be 2.
    # unitpx in \AA unit.
    a_unitpx = (10*scansize_nm[0])/lines_corr
    c_unitpx = (10*scansize_nm[1])/pixels_corr

    return z_corr, a_unitpx, c_unitpx