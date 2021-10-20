
Download the SDSS DR12 BOSS galaxy data, the random particle data, and the Patchy-mock data from the following links and copy them into `$WORK/data/boss/`:

[galaxy_DR12v5_CMASSLOWZTOT_North.fits.gz](https://data.sdss3.org/sas/dr12/boss/lss/)
[galaxy_DR12v5_CMASSLOWZTOT_South.fits.gz](https://data.sdss3.org/sas/dr12/boss/lss/)
[random0_DR12v5_CMASSLOWZTOT_North.fits.gz](https://data.sdss3.org/sas/dr12/boss/lss/)
[random0_DR12v5_CMASSLOWZTOT_South.fits.gz](https://data.sdss3.org/sas/dr12/boss/lss/)
[random1_DR12v5_CMASSLOWZTOT_North.fits.gz](https://data.sdss3.org/sas/dr12/boss/lss/)
[random1_DR12v5_CMASSLOWZTOT_South.fits.gz](https://data.sdss3.org/sas/dr12/boss/lss/)

[Patchy-Mocks-DR12NGC-COMPSAM_V6C.tar.gz](https://data.sdss3.org/sas/dr12/boss/lss/dr12_multidark_patchy_mocks/ )
[Patchy-Mocks-DR12SGC-COMPSAM_V6C.tar.gz](https://data.sdss3.org/sas/dr12/boss/lss/dr12_multidark_patchy_mocks/ )
[Patchy-Mocks-Randoms-DR12NGC-COMPSAM_V6C_x100.tar.gz](https://data.sdss3.org/sas/dr12/boss/lss/dr12_multidark_patchy_mocks/ )
[Patchy-Mocks-Randoms-DR12SGC-COMPSAM_V6C_x100.tar.gz](https://data.sdss3.org/sas/dr12/boss/lss/dr12_multidark_patchy_mocks/ )

Then, unzip the .gz and .tar.gz files:

    gzip -d xxxx.gz (for .gz files)
    tar -zxvf yyyy.tar.gz (for .tar.gz files)

In order to rewrite the galaxy data as 3D Cartesian coordinates, run

    python fits2dat_galaxy.py

For the galaxy random data, run

    python fits2dat_galaxy_random.py

All the resulting files are placed in `$WORK/data/boss/galaxy_DR12v5_CMASSLOWZTOT`

----------

For the Multi-Dark Patchy mock data, run

    python fits2dat_mock.py

For the mock random data, run

    python fits2dat_mock_random.py

All the resulting files are placed in `$WORK/data/boss/Patchy-Mocks-DR12[NGC,SGC]-COMPSAM_V6C_ZBIN[1,2,3]`

