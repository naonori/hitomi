
Here, we perform parameter estimation using only the two-point correlation function. The following summarizes the details of this analysis.

0. Please copy and paste `MontePython.py` from `WORK/cosmo/montepython/montepython` to `galaxy` and `mock`.
1. We use only large scales above 80 Mpc/h
2. We vary only two parameters: linear bias and linear growth rate.
3. We include the effect of nonlinear BAO decay in the template model.
4. We ignore the AP effect.
5. We divide the LOWZCMASSTOT sample obtained from BOSS into four samples: NGC at =0.38, SGC at =0.38, NGC at =0.61, SGC at =0.61. Then, we perform parameter estimation in the following six patterns.
    - NGC at =0.38
    - SGC at =0.38
    - NGC at =0.61
    - SGC at =0.61
    - NGC + SGC at =0.38
    - NGC + SGC at =0.61
6. `galaxy` performs the analysis using actual BOSS data, while `mock` repeats the same analysis performed in `galaxy` on 100 Patchy mock catalogs.

