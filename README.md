# HomCAR

This repository contains the R code to reproduce the analyses of the paper entitled *"A proposal for homoskedastic modelling with conditional auto-regressive distributions"* (Martínez-Beneito et al., 2025).

The results of the models presented in this paper were obtained using the R-INLA stable version 24.06.27 on R 4.4.1

## Table of contents

-   [R code](#r-code)
-   [Acknowledgements](#acknowledgements)
-   [References](#references)


# R code

This folder contains the following files:

-   [**Carto_files.Rdata**](https://github.com/spatialstatisticsupna/HomCAR/blob/master/Rcode/Carto_files.Rdata)

    An .Rdata file with `sf` objects containing the cartographies of municipalities in the autonomous regions of the Valencian region (`CartoVR`), Castile and Leon (`CartoCL`), Aragon (`CartoAR`) and Castile-La Mancha (`CartoCM`).

-   [**Results_ICARmodels.Rdata**](https://github.com/spatialstatisticsupna/HomCAR/blob/master/Rcode/Results_ICARmodels.Rdata)

    An .Rdata file with the posterior mean estimates of log-risks for 100 causes of death across municipalities in the Valencian region (`Res.VR`), Castile and Leon (`Res.CL`), Aragon (`Res.AR`) and Castile-La Mancha (`Res.CM`), obtained from models fitted using the traditional ICAR prior.

-   [**Results_HomCARmodels.Rdata**](https://github.com/spatialstatisticsupna/HomCAR/blob/master/Rcode/Results_HomCARmodels.Rdata)

    An .Rdata file with the posterior mean estimates of log-risks for 100 causes of death across municipalities in the Valencian region (`Res.VR.HomCAR`), Castile and Leon (`Res.CL.HomCAR`), Aragon (`Res.AR.HomCAR`) and Castile-La Mancha (`Res.CM.HomCAR`), obtained from models fitted using the proposed HomCAR prior.
    
-   [**Section2.R**](https://github.com/spatialstatisticsupna/HomCAR/blob/master/Rcode/Section2.R)

    R script to reproduce the results from Section 2 (*The ICAR distribution and its marginal variances*) of the paper, computing the marginal variances of the ICAR distribution for municipalities in four autonomous regions of Spain and comparing them with the empirical variances of the posterior mean log-risk estimates across 100 causes of death analyzed in each region.
    
-   [**Section4.R**](https://github.com/spatialstatisticsupna/HomCAR/blob/master/Rcode/Section4.R)

    R script to reproduce the results from Section 4 (*A practical assessment of the HomCAR distribution*) of the paper, presenting a simulation study that compares the performance of the proposed HomCAR distribution with the traditional ICAR prior
    
    
# Acknowledgements

This research was supported by the projects PID2020-113125RBI00/MCIN/AEI/10.13039/501100011033, PID2022-136455NBI00/MCIN/AEI/10.13039/501100011033 and PID2024-155382OBI00/MICIU/AEI/10.13039/501100011033 funded by the Ministry of Science, Innovation and Universities of Spain and FEDER (UE) and grant CIAICO PID2022-136455NB-I00 funded by Dirección General de Ciencia e Investigación (Generalitat Valenciana).

![plot](https://github.com/spatialstatisticsupna/HomCAR/blob/main/miciu-aei.png)


# References

Martínez-Beneito, M.A., Adin, A., Goicoa, T., and Ugarte, M.D. (2025). A proposal for homoskedastic modelling with conditional auto-regressive distributions. *Statistics and Medicine*, accepted for publication.
