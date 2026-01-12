This folder contains the supporting materials for the paper Two-dimensional data examination by exploratory functional data analysis to improve detection of scattered soil contamination by Cu-bearing pesticides (Matys Grygar et al., 2026).
The folder contains:
    - review_code.R: the file to reproduce figures from the paper;
    - Initially_proc_data.csv: .csv file with the raw underlying data;
    - Cu_arithmetic_new.csv*, B_Cu_geometric.csv*: B-spline coefficients for the Cu arithmetic and geometric marginals. These B-spline coefficients, combined with a proper B-spline basis (in review_code.R), result in the L^2_0 representation of the given marginals.
    - B_overall_kutna_hora.csv*, B_overall_plzen_sever.csv*, B_overall_rakovnik.csv*, B_overall_usti_nad_orlici.csv*: the bivariate B-spline representation (not decomposed) of a sample of specific districts. 


*The smoothing of the raw data and its decomposition was performed using the codes in https://github.com/skorst01/Bivariate-Compositional-Splines. The choice of parameters is specified in ...
