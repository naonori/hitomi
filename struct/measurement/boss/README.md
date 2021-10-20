
Here we compute the power spectrum (pk), the two-point correlation function (2pcf), the bispectrum (bk), the three-point correlation function (3pcf), and their BAO reconstructed versions from the BOSS galaxy data. In addition, we also calculate the two-point and three-point correlation functions of the window function for the BOSS survey using the random particle distribution.

----------------

Our algorithm is based on the FFT. Therefore, the results depend on the number of meshes separating the FFT boxes containing the particles. For a fair comparison with previous studies, we adopt the box size and the number of meshes used by [Beutler et al. 2016](https://academic.oup.com/mnras/article/466/2/2242/2712530). Download Beutler+2016's results from the link below and unzip them. 
 
    https://data.sdss.org/sas/dr12/boss/papers/clustering/Beutler_etal_DR12COMBINED_fullshape_powspec.tar.gz

Then, open the following four files contained in the extracted directory.

    Beutleretal_pk_monopole_DR12_NGC_z1_prerecon_120.dat
    Beutleretal_pk_monopole_DR12_NGC_z3_prerecon_120.dat
    Beutleretal_pk_monopole_DR12_SGC_z1_prerecon_120.dat
    Beutleretal_pk_monopole_DR12_SGC_z3_prerecon_120.dat

Users can find the following box sizes and mesh numbers:

for NGC_z1,
    Lx = 1350 Ly = 2450 Lz = 1400
    nx = 250 ny = 460 nz = 260

for NGC_z3,
    Lx = 1800 Ly = 3400 Lz = 1900
    nx = 340 ny = 650 nz = 360

for SGC_z1,
    Lx = 1000 Ly = 1900 Lz = 1100
    nx = 190 ny = 360 nz = 210

 for SGC_z3,
    Lx = 1000 Ly = 2600 Lz = 1500
    nx = 190 ny = 500 nz = 280.

Increasing the number of meshes allows for more accurate measurements at smaller scales.

----------------

Copy the three files `Makefile`, `main.cpp`, and `default_param.ini` in `$WORK/cosmo/hitomi_measurement`, and paste them into the user's preferred location. Users can edit these files to suit their convenience.


