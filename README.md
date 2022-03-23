# Hitomi

| Main developer: | Naonori Sugiyama <<nao.s.sugiyama@gmail.com>>|
| ----- | ----- |
| **Author:** | **Naonori Sugiyama <<nao.s.sugiyama@gmail.com>>**|
| Licence: | MIT |

The purpose of the `Hitomi` project is to provide a comprehensive set of codes for cosmological analysis of anisotropic galaxy distributions using two- and three-point statistics: two-point correlation function (2PCF), power spectrum, three-point correlation function (3PCF), and bispectrum.

Specifically, this set of codes is divided into seven phases:

    [0] install
    [1] data
    [2] measurement
    [3] model
    [4] covariance
    [5] analysis
    [6] figure

`Hitomi` has `hitomi_measurement`, which measures cosmological statistics from galaxy samples, and `hitomi_theory`, which computes the corresponding theoretical model. Each code is written in c++, and `hitomi_theory` is further wrapped in python.

This set of codes is a collection of programs that were used to write the following papers.

1. [Naonori Sugiyama, Maresuke Shiraishi, and Teppei Okumura, 2018, MNRAS, 473, 2737](https://academic.oup.com/mnras/article/473/2/2737/4111164)
2. [Naonori Sugiyama, Shun Saito, Florian Beutler, and Hee-Jong Seo, 2019, MNRAS, 484, 364](https://academic.oup.com/mnras/article/484/1/364/5222671)
3. [Naonori Sugiyama, Shun Saito, Florian Beutler, and Hee-Jong Seo, 2020, MNRAS, 497, 1684](https://academic.oup.com/mnras/article/497/2/1684/5867793)
4. [Masato Shirasaki, Naonori Sugiyama, Ryuichi Takahashi, and Francisco-Shu Kitaura, 2021, Phys. Rev. D 103, 023506](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.023506)
5. [Naonori Sugiyama, Shun Saito, Florian Beutler, and Hee-Jong Seo, 2021, MNRAS, 501, 2862](https://academic.oup.com/mnras/article/501/2/2862/6019894)
6. Naonori Sugiyama and others, 2022 in preparation.

### Note

There is no detailed documentation for using our code, but all the steps are available as videos on **YOUTUBE**. Our code set is not wholly idealized, so users may need to do some manual work from time to time to complete the cosmological analysis. (e.g., creating some directories, changing input parameters manually, etc.). It is not easy to explain all those details in the documentation, but the video will make it clear.

# Demo

Here, we give a demonstration to complete the analysis that tests DHOST theories, a type of modified gravity theory, using the anisotropic 3PCF.

- [Constraints on modified gravity theories](DEMO_MG.md)

# Description

This project is still in progress. The ultimate goal of this project is to implement all of the program listed below.

The following list has "**implemented**" or "**not yet**" for each of the cases that have already been implemented in `Hitomi` and those that have not yet been implemented. In addition, cases that have already been calculated in the above-published papers but not yet implemented in `Hitomi` are appended with "**possible but not yet**". (as of February 08, 2022)

All `Hitomi` code demonstrations use the galaxy sample derived from the Baryon Oscillation Spectroscopic Survey ([BOSS](http://www.sdss3.org/science/boss_publications.php)), part of the Data Release 12 of the Sloan Digital Sky Survey III ([SDSS-III DR12](https://www.sdss.org/dr12/)), and the corresponding Multi-Dark Patchy mock catalogues([Patchy mock](http://www.skiesanduniverses.org/page/page-3/page-15/page-9/)). These sample codes would be easily applicable to galaxy samples from other galaxy redshift surveys.

### Note

"**implemented**" in the list below is a link to a video on **YOUTUBE** that explains how to use the code.

## [0] Install

Since `Hitomi` uses various external codes, users need to install them and then install the codes related to `Hitomi`.

- install all necessary libraries ([implemented](https://youtu.be/vlP7XIXZsUM))

The first step is to create a **WORK** directory where you will actually work. (In our demo case, `WORK=/mwork0/sugiymnn/WORK`.) 

    $> mkdir WORK

All the contents of the downloaded `Hitomi` directory will be properly copied into `WORK`. Therefore, users can put `Hitomi` anywhere they like. The demo case places it in the same directory as `WORK`.

    $> ls
    $> hitomi/ WORK/

Download the following source files and copy them into `hitomi/env`:
  
 - [Anaconda3-2021.05-Linux-x86_64.sh](https://www.anaconda.com/)
    - `hitomi_theory` and `class` require `python3`. `Anaconda3` is a convenient way to build `python3` and its surrounding environment.
 - [fftw-3.3.9.tar.gz](http://www.fftw.org/index.html)
    - Both `hitomi_measurement` and `hitomi_theory` require `fftw3`.
 - [gsl-2.7.tar.gz](https://www.gnu.org/software/gsl/)
    - Both `hitomi_measurement` and `hitomi_theory` require `gsl`.
 - [Cuba-4.2.1.tar.gz](http://www.feynarts.de/cuba/)
    - `hitomi_theory` uses `pycuba` to perform multiple integrations numerically. The original `cuba` library is required to install `pycuba`.
 - [MultiNest-master.zip](https://github.com/farhanferoz/MultiNest.git)
    - `montepython` requires `pymultinest`, and the original `multinest` is needed to install `pymultinest`.
 - [cmake-3.21.1.tar.gz](https://cmake.org/)
    - To install `pymultinest`, `cmake` is required.
 - [lapack-3.10.0.tar.gz](http://www.netlib.org/lapack/)
    - To install `pymultinest`, `lapack` is required.
 - [PyMultiNest-2.10.zip](https://github.com/JohannesBuchner/PyMultiNest/tree/v2.10)
    - `pymultinest` installs `pymultinest` and `pycuba` at the same time. `montepython` and `hitomi_theory` require `pymultinest` and `pycuba`, respectively.
 - [astropy-4.3.zip](https://github.com/astropy/astropy/tree/v4.3)
    - `astropy` is used to read fits files containing the BOSS galaxy distribution data.
 - [class_public-3.0.1.zip](https://github.com/lesgourg/class_public/tree/v3.0.1)
    - `class` is used to compute linear power spectra and various cosmological functions such as the linear growth factor and the linear growth rate. 
 - [montepython_public-3.4.zip](https://github.com/brinckmann/montepython_public/tree/3.4)
    - `montepython` is used to estimate cosmological parameters of interest.
 - [FFTLog-master.zip](https://github.com/slosar/FFTLog)
    - `hitomi_theory` uses `FFTLog` to compute Hankel transforms from the power spectrum and bispectrum to the 2PCF and 3PCF, respectively, and vice versa. `hitomi_theory` requires that the contents of `FFTLog` written in C be rewritten to be wrapped in python. Therefore, `FFTLog` must be installed anew by users.


After the above libraries are installed correctly, users can install `hitomi_measurement` and `hitomi_theory` in `hitomi/src`. All the libraries can be installed via `install.sh`.

    ./install.sh
  
This run will install all the necessary libraries in `$WORK/cosmo`. 

    $WORK/cosmo> ls
    $WORK/cosmo> anaconda3/  astropy/  class/  cmake/  cuba/  fftlog/  fftw3/  gsl/  hitomi_measurement/  hitomi_theory/  lapack/  montepython/  multinest/  pymultinest/

If `Hitomi` is installed correctly, the final directory structure will be:

    $WORK> ls
    $WORK> analysis/ cosmo/ covariance/ data/ figure/ measurement/ model/

### Note

1. We have assumed that users do not have an external network, so all libraries are installed from the source files. However, if users have access to an external network, it will be easier to build the environment. In that case, please refer to the text written in "install.sh" to install the above 14 libraries in the user's environment.

2. If users have already installed some of the above libraries (e.g., `fftw3`, `gsl`, `astropy`, etc.) in their environment, they do not need to re-install them in `$WORK/cosmo`.
    
3. Users need to specify where to install `fftw` and `gsl` in the user's environment in `hitomi/src/hitomi_measurement/Makefile` and `hitomi/src/hitomi_theory/setup.py`.

## [1] Data

The position data of the observed galaxies are given in terms of right ascension, declination, and redshift. In order to measure cosmological statistics such as the power spectrum, bispectrum, 2PCF, and 3PCF, these data need to be transformed into Cartesian coordinates (x,y,z) using fiducial cosmological parameters.

1. BOSS galaxy sample in the Cartesian coordinate system ([implemented](https://youtu.be/J-0u_gUwTqE))
2. BOSS random sample in the Cartesian coordinate system ([implemented](https://youtu.be/SBLs9TCbALk))
3. Patchy mock galaxy sample in the Cartesian coordinate system ([implemented](https://youtu.be/-aOCTbsruuM))
4. Patchy mock random sample in the Cartesian coordinate system ([implemented](https://youtu.be/-1z_yKcSegA))

## [2] Measurement

`Hitomi` can measure the Legendre-expanded 2PCF and power spectrum from an observed sample of galaxies (e.g., [SDSS DR12 BOSS](https://data.sdss.org/sas/dr12/boss/lss/)). It can also measure the 3PCF and bispectrum expanded using the Tripolar spherical harmonic (TripoSH) function proposed by [Sugiyama et al. 2019](https://academic.oup.com/mnras/article/484/1/364/5222671). All of these calculations are performed using the Fast Fourier Transform (FFT) algorithm ([fftw](http://www.fftw.org/index.html)). Therefore, the window 2PCF and window 3PCF also need to be measured to account for corrections to the survey geometry.

1. power spectrum ([implemented](https://youtu.be/Xe_Wvj1clck))
2. 2PCF ([implemented](https://youtu.be/yidW5qRTDY0))
3. window 2PCF ([implemented](https://youtu.be/gFz1VjWfMRg))
4. bispectrum ([implemented](https://youtu.be/9hVyLKhltv0))
5. 3PCF ([implemented](https://youtu.be/axCqkd2KNAU))
6. window 3PCF ([implemented](https://youtu.be/FJS9v9O03ls))

----------------------

`Hitomi` also implements a simple reconstruction of the galaxy distribution ([Shirasaki et al. 2021](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.023506)). Therefore, the reconstructed versions of all the statistics listed above can also be measured.

7. power spectrum after reconstruction ([implemented](https://youtu.be/c5K_kqq7kdc))
8. 2PCF after reconstruction ([implemented](https://youtu.be/SpKtpmg6Pzo))
9. window 2PCF after reconstruction ([implemented](https://youtu.be/Qsls-kc5xZE))
10. bispectrum after reconstruction ([implemented](https://youtu.be/w7v6Tc-Krg0))
11. 3PCF after reconstruction ([implemented](https://youtu.be/fpR53XZA_E4))
12. window 3PCF after reconstruction ([implemented](https://youtu.be/QZpPvtAeNGk))
 
----------------------

When assuming a global plane-parallel approximation with periodic boundary conditions, often considered for box simulation data, the algorithm for two- or three-point statistics using FFT is simplified. `Hitomi` implements programs that load the dark matter particle data produced by [Gadget2](https://wwwmpa.mpa-garching.mpg.de/gadget/) and the halo data produced by [rockstar](https://code.google.com/archive/p/rockstar/), and can measure two-point and three-point statistics for these dark matter particles and halos before and after reconstruction.
  
13. power spectra for dark matter and halos before and after reconstruction ([implemented](https://youtu.be/ufflObL0-bo))
14. 2PCFs for dark matter and halos before and after reconstruction ([implemented](https://youtu.be/I56s5QECpt8))
15. bispectra for dark matter and halos before and after reconstruction ([implemented](https://youtu.be/iR3aNJPhHvY))
16. 3PCFs for dark matter and halos before and after reconstruction ([implemented](https://youtu.be/h5oWVmSOGww))

----------------------

`Hitomi` creates new output directories corresponding to the number of the simulation data to be read when handling a large number of simulation data with the same cosmological parameters but generated by different initial conditions, such as mock simulations. One example is to compute multiple 2PCFs and 3PCFs from the Patchy mock data.

17. 2PCFs measured from the Patchy mock catalogues ([implemented](https://youtu.be/n6_ifs9bUeM))
18. 3PCFs measured from the Patchy mock catalogues ([implemented](https://youtu.be/JVV7ybXJUPk))

----------------------

`Hitomi` is basically a serial code, but it also implements MPI parallelization. `Hitomi` uses MPI to read multiple different input parameters simultaneously. In this case, no data is exchanged between MPIs, and parallelization can be achieved with 100% efficiency. For example, in the case of MPI parallelization with size=40, `Hitomi` can simultaneously load the mock simulation data from number 1 to 40 and compute the corresponding cosmological statistics (e.g., the 2PCF, power spectrum, 3PCF, and bispectrum) for these simulation data at once. As another example, by assigning one of the two scales on which the bispectrum and 3PCF depend on each MPI rank, it is possible to measure the bispectrum and 3PCF at high speed. With these examples as a guide, users can implement their own MPI parallelization for the various input parameters used by `Hitomi`.

19. 2PCFs measured from the Patchy mock catalogues using MPI ([implemented](https://youtu.be/eUgjqUwj3QY))
20. fast measurement of the 3PCF from a single galaxy sample using MPI ([implemented](https://youtu.be/DzeZf-nkkfk))

## [3] Model

`Hitomi` calculates the power spectrum and bispectrum using perturbation theory. By calculating the Hankel transforms of the power spectrum and bispectrum using [FFTLog](https://github.com/slosar/FFTLog), it can also calculate the 2PCF and 3PCF. The computed power spectrum and bispectrum are based on so-called tree-level solutions, which means that linear solutions are used for the power spectrum, and perturbation solutions up to the second-order are used for the bispectrum. Redshift space distortion (RSD) effects and bias effects are also included, as appropriate, in the respective theoretical models for the power spectrum and the bispectrum. For reference, the power spectrum, bispectrum, 2PCF, and 3PCF of the smoothed version without baryon acoustic oscillations (BAO), the so-called "no-wiggle" version, can also be calculated.

1. linear power spectrum and 2PCF ([implemented](https://youtu.be/WnWBEZSv-1Q))
2. linear power spectrum and 2PCF without BAO ([implemented](https://youtu.be/daajT-dgoUo))
3. tree-level bispectrum and 3PCF ([implemented](https://youtu.be/mZwCGbrDJFw))
4. tree-level bispectrum and 3PCF without BAO ([implemented](https://youtu.be/Kfivlm3XQqU))

----------------------

`Hitomi` can calculate the effect of nonlinear damping in BAO to incorporate higher-order perturbation effects beyond the tree level. A Gaussian function can describe this effect. In the case of the power spectrum, the model is given by P = (Plin-Pnw) * exp(-k^2 Sigma^2) + Pnw ([Eisenstein et al. 2007](https://iopscience.iop.org/article/10.1086/518755)), where Plin is the linear power spectrum, and Pnw is its no-wiggle version. The bispectrum case is given by [equation (1) in Sugiyama et al. 2021](https://academic.oup.com/mnras/article/501/2/2862/6019894), using Plin and Pnw.

5. linear power spectrum and 2PCF with nonlinear BAO damping ([implemented](https://youtu.be/RaUWRiBkVy0))
6. tree-level bispectrum and 3PCF with nonlinear BAO damping ([implemented](https://youtu.be/-ptwDhsHFs8))

----------------------

`Hitomi` can also calculate the power spectrum, bispectrum, 2PCF, and 3PCF in the reconstructed case. The same values of the linear bias parameter (b1) and the Gaussian smoothing parameter (R) entered when reconstructing the galaxy distribution can also be entered for the theoretical model to calculate a model that has a one-to-one correspondence with the measurement results. Since the reconstruction does not change the linear theory, only the effect of moderating the nonlinear damping effect of BAO is taken into account in the case of the power spectrum. In the reconstructed bispectrum case, the theoretical model considers both the impact of changing the second-order perturbation solution and of moderating the BAO damping.

7. linear power spectrum and 2PCF with reconstructed BAO damping ([implemented](https://youtu.be/yKy0o5WRLzc))
8. reconstructed tree-level bispectrum and 3PCF with reconstructed BAO damping ([implemented](https://youtu.be/uCRj_XyjXII))

----------------------

Theoretical models need to be computed fast in order to estimate the cosmological parameters. In particular, the bispectrum and 3PCF are more challenging to compute fast because they depend on two different scales and the number of data bins to be computed tends to increase compared to the case of the power spectrum and 2PCF. [Section 4.3 of Sugiyama et al. 2021](https://academic.oup.com/mnras/article/501/2/2862/6019894) proposed a method that decomposes the bispectrum in terms of each parameter of interest and that the parameter-independent terms are calculated in advance to speed up the calculation. 

9. parameter-decomposed power spectra and 2PCFs with nonlinear BAO damping ([implemented](https://youtu.be/W_nTnuDKcEU))
10. parameter-decomposed bispectra and 3PCFs with nonlinear BAO damping ([implemented](https://youtu.be/nuHTn9unSmI))

For the Alcock-Paczyn'ski (AP) effect, [Section 4.4 in Sugiyama et al. 2021](https://academic.oup.com/mnras/article/501/2/2862/6019894) introduced a fast method to calculate the bispectrum that includes the AP effect, using an approximation method that decomposes the bispectrum with the parameters that characterize the AP effect.．

11. approximate bispectrum and 3PCF with nonlinear BAO damping, including the AP effect (***possible but not yet***)

----------------------

`Hitomi` can compute the effect of three types of primordial non-Gaussianity (PNG) in the bispectrum and 3PCF: local, equilateral, and orthogonal types. However, it considers only the simplest case, where the coupling between PNGs and gravitational nonlinear effects is neglected. In this case, the PNG signals in the bispectrum and 3PCF are unchanged before and after reconstruction.

12. bispectra and 3PCFs generated from the local, equilateral, and orthogonal PNGs ([implemented](https://youtu.be/vb8Uxv1idg0))

----------------------

In any theoretical model that `Hitomi` calculates, higher-order nonlinear effects, the so-called loop corrections, which are essential at small scales, are not considered. (Except for the impact of nonlinear BAO damping). Therefore, the limit of applicability of the theoretical model given by `Hitomi` is limited to large scales. In fact, in the cosmological analysis of [Sugiyama et al. 2021](https://academic.oup.com/mnras/article/501/2/2862/6019894), only large scales above 80 Mpc/h were used for the 2PCF and 3PCF. That paper showed that the addition of the 3PCF can lead to improved results compared to the existing results using only the 2PCF, even considering only such large scales if focusing only on the parameters related to BAO, i.e., the Hubble parameter and the angular diameter distance. However, to better constrain parameters, such as fsigma8, for neutrino mass limits and testing modified gravity theories, it is essential to build models that are applicable up to small scales. Our goal is to develop such a theoretical model in the future.

13. power spectrum and 2PCF with loop corrections (**not yet**)
14. bispectrum and 3PCF with loop corrections (**not yet**)

The nonlinear loop correction also affects the reconstruction because the reconstruction method partially removes the effect of the nonlinear gravitational evolution.

15. reconstructed power spectrum and 2PCF with loop corrections (**not yet**)
16. reconstructed bispectrum and 3PCF with loop corrections (**not yet**)

Furthermore, the loop correction results in a gravitational nonlinear coupling with PNGs.

17. bispectra and 3PCFs generated from the local, equilateral, and orthogonal PNGs with loop corrections (**not yet**)
18. reconstructed bispectra and 3PCFs generated from the local, equilateral, and orthogonal PNGs with loop corrections (**not yet**)

## [4] Covariance Matrix

In the galaxy cosmological analysis, the currently dominant method for calculating the covariance matrix is to prepare a large number of mock samples that reproduce the observed galaxy distribution and calculate it directly using a large number of power spectra, bispectra, 2PCFs, and 3PCFs measured from these samples. In fact, in [Sugiyama et al. 2021](https://academic.oup.com/mnras/article/501/2/2862/6019894), the covariance matrices of the 2PCF and 3PCF were calculated using 2048 mock galaxy catalogues from the [Multi-Dark Patchy mock](https://data.sdss.org/sas/dr12/boss/lss/dr12_multidark_patchy_mocks/) prepared for the SDSS DR12 BOSS analysis.

0. covariance matrices computed from the 2048 Patchy mock catalogues ([implemented](DEMO_MG.md))

----------------------

However, the number of mock catalogues required is huge to obtain reliable inverse covariance matrices for bispectra and 3PCFs, which have far more data bins than power spectra and 2PCFs. In addition, for example, [DESI](https://www.desi.lbl.gov/) and [PFS](https://pfs.ipmu.jp/) are planning to observe a more comprehensive redshift range and a more extensive survey area than BOSS, and the computational cost of generating mock catalogues is increasing.

As an alternative method to solve this problem, the `Hitomi` project promotes the calculation of covariance matrices using perturbation theory. The covariance matrix computed by perturbation theory corresponds to the one computed by an infinite number of simulations, and the problem of the number of mock catalogues no longer arises. It is also much less computationally expensive than creating a large number of mock catalogues. However, the question remains whether perturbation theory can explain nonlinear gravity clustering accurately enough. With the development of previous studies, a study on the cosmological analysis of BOSS galaxies using the covariance matrices computed by perturbation theory in the power spectrum case has already been done ([Wadekar et al.(2020)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.102.123521)). [Sugiyama et al. 2020](https://academic.oup.com/mnras/article/497/2/1684/5867793) has shown that the perturbation theory approach can reproduce the bispectrum covariance matrix measured by the Patchy mock with reasonable accuracy, but there is still room for improvement. Based on the results of [Sugiyama et al. 2020](https://academic.oup.com/mnras/article/497/2/1684/5867793), one goal of the `Hitomi` project is to build an analytical model of the covariance matrix that is reliable enough to be applied to the cosmological analysis of galaxy bispectra (3PCFs).

1. power spectrum covariance

    - gaussian part

        - linear power spectrum (***possible but not yet***)
        - gaussian shot-noise (***possible but not yet***)
        - linear bias (***possible but not yet***)
        - non-linear power spectrum (**not yet**)
        - survey window effect (**not yet**)
                 
    - non-gaussian part

        - trispectrum consisting of perturbations of the density and velocity fields up to the third-order (***possible but not yet***)
        - non-gaussian shot-noise (***possible but not yet***)
        - linear bias (***possible but not yet***)
        - non-linear bias parameters up to the third-order (**not yet**)
        - super-sample effect (**not yet**)
        - loop corrections (**not yet**)

2. cross-covariance between the power spectrum and bisectrum

    - non-gaussian part

        - 5-point spectrum consisting of perturbations of the density and velocity fields up to the fourth-order (***possible but not yet***)
        - non-gaussian shot-noise (***possible but not yet***)
        - linear bias (***possible but not yet***)
        - non-linear bias parameters up to the fourth-order (**not yet**)
        - super-sample effect (**not yet**)
        - loop corrections (**not yet**)

2. bispectrum covariance

    - gaussian part

        - linear power spectrum (***possible but not yet***)
        - gaussian shot-noise (***possible but not yet***)
        - linear bias (***possible but not yet***)
        - non-linear power spectrum (**not yet**)
        - survey window effect (**not yet**)

    - non-gaussian part

        - 6-point spectrum consisting of perturbations of the density and velocity fields up to the fifth-order (***possible but not yet***)
        - non-gaussian shot-noise (***possible but not yet***)
        - linear bias (***possible but not yet***)
        - non-linear bias parameters up to the fifth-order (**not yet**)
        - super-sample effect (**not yet**)
        - loop corrections (**not yet**)

## [5] Analysis

The current `Hitomi` project envisages four types of data analysis: i.e., BAO, RSD, modified gravity, and PNG.

`Hitomi` uses [montepython3](https://github.com/brinckmann/montepython_public/tree/3.4) for parameter estimation.

----------------------

The BAO analysis, widely performed using power spectrum and 2PCF, mainly constrains the Hubble parameter and the angular diameter distance. Furthermore, the RSD analysis can also constrain the growth rate parameter fsigma8 that characterizes the amplitude of the velocity field. The addition of the bispectrum and 3PCF is considered to improve the constraints on these three parameters. Indeed, [Sugiyama et al. 2021](https://academic.oup.com/mnras/article/501/2/2862/6019894) have shown that adding of the 3PCF can improve the constraint on the Hubble parameter by ~ 30%.

1. BAO analysis (***possible but not yet***)

2. RSD analysis (***possible but not yet*** on large scales and **not yet** on small scales)

----------------------

As of February 08, 2022, the primary focus of the `Hitomi` project is to test the modified theory of gravity using the anisotropic component of the 3PCF, i.e., the nonlinear shape of the velocity field. 

In the case of power spectra and 2PCFs, the anisotropic component along the line of sight, e.g., the quadrupole component, produced by RSD, can increase the cosmological information because the density fluctuations of galaxies that correspond to the isotropic component are biased quantities, while the velocity field is considered to be the unbiased physical quantity on large scales. Specifically, the quadrupole component of the power spectrum or 2PCF can measure the growth rate function, fsigma8, which can be used to test the modified theory of gravity.

The same is true for bispectra and 3PCFs; the unique feature of `Hitomi` is that it can treat the anisotropic components of bispectra and 3PCFs both observationally and theoretically. Thus, the nonlinear modified gravity test with the anisotropic 3PCF is an analysis that takes full advantage of `Hitomi` characteristics.

3. modified gravity ([implemented](DEMO_MG.md))

----------------------

The ability to test PNGs is also a feature of the bispectrum and 3PCF. In the case of the power spectrum or 2PCF, PNGs are mainly constrained through a scale-dependent bias. This effect is proportional to 1/k^2 and disappears at small scales. In the case of the bispectrum or 3PCF, it should be possible to test the PNG signal even at small scales where systematic errors can be suppressed.

As shown by [Sugiyama et al. 2020](https://academic.oup.com/mnras/article/497/2/1684/5867793), the statistical error of the bispectrum is dominated by non-Gaussian errors. However, [Shirasaki et al. 2021](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.023506) has shown using simulations that the reconstruction can dramatically reduce these non-Gaussian errors. The reconstruction also reduces the bispectrum signal due to gravitational nonlinearity, but the lowest order PNG signal, which is not coupled by gravity, is unchanged before and after the reconstruction. Therefore, the verification of PNG using the reconstructed bispectrum (or 3PCF) is underway in the `Hitomi` project.
 
4. primordial non-gaussianity (**not yet**)

## [6] Figure

`Hitomi` plans to release all the code used to plot the actual figures used in our papers. This release will make it easier for users to reproduce our results.

It is surprisingly time-consuming to draw beautiful and clear figures. It is a significant part of the research, but it should not be the most crucial part. If users like our figure format, they can save their research time by using it.

## Future Prospects


## Licence

The license of the source code is MIT; see LICENSE.txt. As an additional clause, when using the code in a scientific publication, users are also required to cite the papers as specified in [Bibtex entry](#bibtex-entry). 

## Bibtex Entry

When using `Hitomi` in a publication, please acknowledge the code by citing the following papers. Users should also properly cite the various libraries required by `Hitomi`, such as `Class`, `montepython`, `FFTLog`, etc.

Please cite the following paper when measuring the power spectrum or 2PCF multipoles.

    @article{Sugiyama:2017ggb,
        author = "Sugiyama, Naonori S. and Shiraishi, Maresuke and Okumura, Teppei",
        title = "{Limits on statistical anisotropy from BOSS DR12 galaxies using bipolar spherical harmonics}",
        eprint = "1704.02868",
        archivePrefix = "arXiv",
        primaryClass = "astro-ph.CO",
        doi = "10.1093/mnras/stx2333",
        journal = "Mon. Not. Roy. Astron. Soc.",
        volume = "473",
        number = "2",
        pages = "2737--2752",
        year = "2018"
    }

Please cite the following paper when measuring the bispectrum or 3PCF multipoles.

    @article{Sugiyama:2018yzo,
        author = "Sugiyama, Naonori S. and Saito, Shun and Beutler, Florian and Seo, Hee-Jong",
        title = "{A complete FFT-based decomposition formalism for the redshift-space bispectrum}",
        eprint = "1803.02132",
        archivePrefix = "arXiv",
        primaryClass = "astro-ph.CO",
        doi = "10.1093/mnras/sty3249",
        journal = "Mon. Not. Roy. Astron. Soc.",
        volume = "484",
        number = "1",
        pages = "364--384",
        year = "2019"
    }

Please cite the following paper when calculating the analytical models of the bispectrum covariance matrix.

    @article{Sugiyama:2019ike,
        author = "Sugiyama, Naonori S. and Saito, Shun and Beutler, Florian and Seo, Hee-Jong",
        title = "{Perturbation theory approach to predict the covariance matrices of the galaxy power spectrum and bispectrum in redshift space}",
        eprint = "1908.06234",
        archivePrefix = "arXiv",
        primaryClass = "astro-ph.CO",
        doi = "10.1093/mnras/staa1940",
        journal = "Mon. Not. Roy. Astron. Soc.",
        volume = "497",
        number = "2",
        pages = "1684--1711",
        year = "2020"
    }

Please cite the following paper when calculating the analytical models of the bispectrum or 3PCF multipoles.

    @article{Sugiyama:2020uil,
        author = "Sugiyama, Naonori S. and Saito, Shun and Beutler, Florian and Seo, Hee-Jong",
        title = "{Towards a self-consistent analysis of the anisotropic galaxy two- and three-point correlation functions on large scales: application to mock galaxy catalogues}",
        eprint = "2010.06179",
        archivePrefix = "arXiv",
        primaryClass = "astro-ph.CO",
        doi = "10.1093/mnras/staa3725",
        journal = "Mon. Not. Roy. Astron. Soc.",
        volume = "501",
        number = "2",
        pages = "2862--2896",
        year = "2021"
    }

Please cite the following paper when measuring the reconstructed bispectrum or 3PCF multipoles or calculating their analytical models.

    @article{Shirasaki:2020vkk,
        author = "Shirasaki, Masato and Sugiyama, Naonori S. and Takahashi, Ryuichi and Kitaura, Francisco-Shu",
        title = "{Constraining primordial non-Gaussianity with postreconstructed galaxy bispectrum in redshift space}",
        eprint = "2010.04567",
        archivePrefix = "arXiv",
        primaryClass = "astro-ph.CO",
        doi = "10.1103/PhysRevD.103.023506",
        journal = "Phys. Rev. D",
        volume = "103",
        number = "2",
        pages = "023506",
        year = "2021"
    }

## Acknowledgement

NS thanks to Mike Shengbo Wang for his careful reading of Hitomi's codes and his detailed comments.

