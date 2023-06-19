## dfphase1 1.2.0 (June 2023)

- shewhart: when argument limits=NA (i.e., when the control limits are
  not precomputed), the control limits are computed using a permutation
  approach also for the rank based statistics; hence, they guarantee the
  desired FAP also when the in-control distribution is not absolutely
  continuous.

- shewhart: added control charts based on the Lepage and Cucconi control
  statistics.

- shewhart: default for argument aggregation is now 'mean' (and not
  'median' as before). This was introduced for 'textbook obedience'. 

- mshewhart: added the possibility to take into account the subgroup
  structure of the data when computing the scatter matrix used to
  (inner) standardise the spatial signs and ranks and the multivariate
  signed ranks.

## dfphase1 1.1.4 (September 2021)

- Minor changes to the Rcpp sources to cope with Rcpp's STRICT_R_HEADERS and
  USE_FC_LEN_T/FCONE declarations

- Added CITATION  

## dfphase1 1.1.3 (February 2020)

- Minor changes to the manual pages (references show the DOI and use the
  RD \doi macro,...).

## dfphase1 1.1.2 (February 2020)

- `rsp` uses the *Rcpp* interface as the rest of the package (and not
  the *.C* interface as before).

- Correction of some small problems resulting in some warnings

- Eliminated the dependency from the deprecated (since C++14) `std::shuffle`

## dfphase1 1.1.1 (January 2017)

- Minor modifications for making the package installable on Solaris

## dfphase1 1.1.0 (January 2017) 

- Inclusion of the "mphase1" method based on forward and adaptive LASSO search
       
## dfphase1 1.0.1 (September 2016)

- Minor modifications to enhance portability 


## dfphase1 1.0.0 (August 2016)

- First public release




