# RefInt
RefInt: A curated reference interactome with true positives and negatives for benchmarking protein-protein interactions. 

## Benchmarking your scored PPIs?
To generate ROC curve, a list of positives and negatives is needed. Initially,
from manually curated protein complexs, PPIs within the complex were considered as positives, whilst PPIs between complexes were negatives. But, are they true negatives? Probably not.

## Strategy
1. Maunally curated protein complexes

* Homo sapiens: The comprehensive resource of mammalian protein complexes [(CORUM)](http://mips.helmholtz-muenchen.de/corum/)

* Mus musculus: The comprehensive resource of mammalian protein complexes [(CORUM)](http://mips.helmholtz-muenchen.de/corum/)

## Available Data
**Homo sapiens**
[RefIntHs](https://github.com/zqzneptune/RefInt/blob/master/Data/RefIntHs.RDS)
[RefComplex](https://github.com/zqzneptune/RefInt/blob/master/Data/CORUM_03.09.2018_Hs_Complex_n_2357.RDS)

**Mus musculus**
[RefIntMm](https://github.com/zqzneptune/RefInt/blob/master/Data/RefIntMm.RDS)
[RefComplex](https://github.com/zqzneptune/RefInt/blob/master/Data/CORUM_03.09.2018_Mm_Complex_n_529.RDS)
