[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/thesias/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/thesias/badges/downloads.svg)](https://anaconda.org/bioconda/thesias)
[![Build Status](https://travis-ci.org/daissi/thesias.svg?branch=master)](https://travis-ci.org/daissi/thesias)
<p align="center"><img src="https://raw.githubusercontent.com/daissi/thesias/master/misc/LogoThesias.png" alt="THESIAS"></p>

**THESIAS** (**Testing Haplotype EffectS In Association Studies**) (1) is a popular software for carrying haplotype association analysis in unrelated individuals.
This program is based on the maximum likelihood model described in (2,4) and is linked to the SEM algorithm (2,3).
THESIAS allows the simultaneous estimation of haplotype frequencies and of their associated effects on the phenotype of interest.
In this version, quantitative, qualitative, categorical and survival outcomes can be studied.
Covariate-adjusted haplotype effects as well as haplotype x covariate interactions can be investigated.

## Installation
### From sources
```
make
make install
```

### Bioconda
```
conda install -c bioconda thesias
```

### Debian / Ubuntu
```
apt install thesias
```

### MacOS / MS Windows
Please refer to the binaries [archive](https://github.com/daissi/thesias/releases/download/3.1/ThesiasPackage.zip).

## Changes
- 2020/01/16: **3.1.1**
  - Clarify license with GPL-3+.
  - Update contact email address.
  - New makefile.

- 2007: **3.1**
  - [Download archive](https://github.com/daissi/thesias/releases/download/3.1/ThesiasPackage.zip)

## License
THESIAS is licensed under [GPL-3+](LICENSE).

## References
1. D.A. Trégouët and V. Garelle, *Bioinformatics* (2007), [A new JAVA interface implementation of THESIAS: testing haplotype effects in association studies](https://doi.org/10.1093/bioinformatics/btm058).
2. D.A. Trégouët and L. Tiret, *European Journal of Human Genetics* (2004), [Cox proportional hazards survival regression in haplotype-based association analysis using the Stochastic-EM algorithm](https://doi.org/10.1038/sj.ejhg.5201238).
3. D.A. Trégouët et al., *Annals of Human Genetics* (2004), [A new algorithm for haplotype‐based association analysis: the Stochastic‐EM algorithm](https://doi.org/10.1046/j.1529-8817.2003.00085.x).
4. D.A. Trégouët et al., *Human Molecular Genetics* (2002), [Specific haplotypes of the P-selectin gene are associated with myocardial infarction](https://doi.org/10.1093/hmg/11.17.2015).
