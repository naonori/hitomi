
-------------

This study aims to test the modified theory of gravity by using the two-point and three-point correlation functions of galaxies. In particular, we investigate two classes of scalar-tensor theories: DHOST and Horndeski theories.
 
In these two theories, the linear growth of dark matter density fluctuations in the late-time Universe is scale-independent. In other words, the scale-dependence of the dark matter linear fluctuations has already been determined in the early universe, when CMB was emitted. Therefore, in this study, we fix the shape of the linear power spectrum of dark matter to the LCDM model by Planck. Strictly speaking, the parameters of the LCDM model should be varied to the extent that Planck allows, but for simplicity, we do not do so here.

The angular diameter distance and the Hubble parameter, determined by the Alcock-Paczyn'ski (AP) effect, also have information from the theory of gravity, but we do not consider them here. The reason is that the primary purpose is to investigate the characteristic nonlinear scale dependence of the modified gravity theory through the galaxy 3PCF.

Therefore, this code computes the 2PCF and 3PCF multipoles normalized by sigma8 using the LCDM parameters of Planck and stores the results. By pre-calculating the 2PCF and 3PCF in this way, the MCMC code only needs to load the results when performing parameter estimation, which dramatically reduces the computation time.

To account for the damping of the BAO signal due to nonlinear effects, we adopt the 2PCF template model given by Eisenstein et al. 2007 and the 3PCF template model given by Sugiyama et al. 2021.

As for future development, there is room for calculation of linear power spectra considering DHOST theories, calculation of power spectra including nonlinear effects other than BAO, and data analysis including the AP effect, but these are left for future works. 

---------------
