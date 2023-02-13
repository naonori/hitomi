
Here, we perform parameter estimation using the two and three-point correlation functions in Horndeski theories. The following summarizes the details of this analysis.

0. Please copy and paste `MontePython.py` from `WORK/cosmo/montepython/montepython` to `galaxy` and `mock`.
1. We use only large scales above 80 Mpc/h
2. We vary six parameters: linear bias (b1), linear growth rate (xi_f = ln_{Omega_m} f), non-linear local bias (b2), non-linear tidal bias (bs2), sigma8, velocity tidal term (xi_t = ln_{Omega_m} E_t)
3. We use monopole and quadrupole 2PCFs (xi0 and xi2), two monopole 3PCFs (zeta000 and zeta110), and two quadrupole 3PCFs (zeta202 and zeta112).
4. We include the effect of nonlinear BAO decay in the template model.
5. We ignore the AP effect.
6. We divide the LOWZCMASSTOT sample obtained from BOSS into four samples: NGC at =0.38, SGC at =0.38, NGC at =0.61, SGC at =0.61. Then, we perform parameter estimation in the following nine patterns.
    - NGC at =0.38
    - SGC at =0.38
    - NGC at =0.61
    - SGC at =0.61
    - NGC + SGC at =0.38
    - NGC + SGC at =0.61
    - NGC at =0.38 + 0.61
    - SGC at =0.38 + 0.61
    - NGC + SGC at 0.38 + 0.61
7. `galaxy` performs the analysis using actual BOSS data, while `mock` repeats the same analysis performed in `galaxy` on 100 Patchy mock catalogs.
8. `galaxy_rescaled` computes the case where the covariance matrix for the samples at z=0.38 are rescaled so that the resulting p-values become to be acceptable.

