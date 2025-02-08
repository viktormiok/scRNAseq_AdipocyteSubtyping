![](https://img.shields.io/badge/language-Python-orange.svg) ![version](https://img.shields.io/badge/GiHub_version-1.1.0-519dd9) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/viktormiok/scRNAseq_AdipocyteSubtyping) ![GitHub issues](https://img.shields.io/github/issues/viktormiok/scRNAseq_AdipocyteSubtyping)

![dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-orange)  	![commit](https://img.shields.io/github/last-commit/viktormiok/scRNAseq_AdipocyteSubtyping) ![GitHub](https://img.shields.io/github/license/viktormiok/scRNAseq_AdipocyteSubtyping)

[![Edit with Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/viktormiok/scRNAseq_AdipocyteSubtyping) 



- [Overview](#overview)
- [Implementation](#implementation)
- [Data and notebooks](#data-and-notebooks)
- [License](#license)
- [References](#references)

## Overview

Pathological changes in adipose tissue are a major underlying cause of systemic insulin resistance and metabolic syndrome. However, the mechanisms underlying this functional impairment are unknown due to an inability to identify and target dysfunctional adipocytes. Based on our previous work and ongoing projects in the laboratory, we conclude that white and brown adipose tissues are functionally and developmentally diverse and that differences between cell populations within these different depots can be detected at the cell surface. 
Thus, our primary objective is identifying individual cell populations and signaling pathways responsible for initiating local and subsequent systemic insulin resistance. 


## Implementation

Python pipeline for analysis and visualization of scRNA-seq data from adolescent (2w) and adult (8w) mice of brown, perigonadal, and subcutaneous adipose depots.

<img src="https://user-images.githubusercontent.com/22052679/148747517-60266170-396b-43e2-82cd-a226077820dc.png" align="top" height="700" width="600">

## Data and notebooks
All the data and notebooks required for performing temporal integrative genomics analysis and publishing in the reference articles have been deposited in the National Center for Biotechnology Information Gene Expression Omnibus (GEO). They are accessible through the GEO Series accession numbers:
| Publication     | GEO number | Notebook |
| ------------- | ------------- | ------------- |
| [Suwandhi et. al.](https://doi.org/10.1038/s41467-021-21826-9)  | [__`GSE164350`__](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4117045)  | [asc1_white_vs_beige_adipocyte_fate.ipynb](https://github.com/viktormiok/scRNAseq_AdipocyteSubtyping/blob/main/notebooks/asc1_white_vs_beige_adipocyte_fate.ipynb) |
| [Karlina et. al.](https://www.life-science-alliance.org/content/4/1/e202000924)  | [__`GSE161447`__](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138079)  |[brown_adipocytes_subtypes.ipynb](https://github.com/viktormiok/scRNAseq_AdipocyteSubtyping/blob/main/notebooks/brown_adipocytes_subtypes.ipynb) |
| [Yan et. al.](https://diabetesjournals.org/diabetes/article/72/Supplement_1/1666-P/150177)  | [__`GSE...`__]  |[adipocyte_scRNAseq_imuneCells.ipynb](https://github.com/viktormiok/scRNAseq_Adipocyte_ImmuneCells/blob/main/notebooks/adipocyte_scRNAseq_imuneCells.ipynb) |
| Altun et. al.| [__`GSE...`__]  |[igf1_preadipocytes_full_analysis_v0_.ipynb](https://github.com/viktormiok/scRNAseq_AdipocyteSubtyping/blob/main/notebooks/asc1_white_vs_beige_adipocyte_fate.ipynb) |

To access one of the data sets for instance GSE164350 you need to run the code below. Unpacking the data requires tar and gunzip, which should already be available on most systems.

```
cd ../  #To get to the main GitHub repo folder
mkdir -p data/tigaR_data_analysis/
cd data/tigaR_data_analysis/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE164350/suppl/GSE164350_RAW.tar
mkdir GSE164350_RAW
tar -C GSE164350_RAW -xvf GSE164350_RAW.tar
gunzip GSE164350_RAW/*_Regional_*
```

## License

__`scRNAseq_AdipocyteSubtyping`__ repository is distributed under the GPL-3.0 License. Please read the license before using information from the repository distributed in the `LICENSE` file.

# References

Publications:

- Suwandhi, L., Altun, I., Karlina, R., Miok, V., Wiedemann, T., Fischer, D., Walzthoeni, T., Lindner, C., Böttcher, Anika., Heinzmann, S.S., Israel, A., Khalil, A. E. M. M., Braun, A., Pramme-Steinwachs, I., Burtscher, I., Schmitt-Kopplin, P., Heinig, M., Elsner, M., Lickert, Heiko., Theis, F. J., Ussar, S. (2021), "[Asc-1 regulates white versus beige adipocyte fate in a subcutaneous stromal cell population](https://doi.org/10.1038/s41467-021-21826-9)", *Nature communications*, 12(1), 1-12.

- Karlina, R., Lutter, D., Miok, V., Fischer, D., Altun, I., Schöttl, T., Schorpp, K., Israel, A., Cero, Cheryl., Johnson, J.W., Kapser-Fischer, I., Böttcher, A., Keipert, S., Feuchtinger, A., Graf, E., Strom, T., Walch, A., Lickert, H., Walzthoeni, T., Heinig, M., Theis, F.J., García-Cáceres, C., Cypess, A.M., Ussar, S. (2021). "[Identification and characterization of distinct brown adipocyte subtypes in C57BL/6J mice](https://www.life-science-alliance.org/content/4/1/e202000924)". *Life science alliance*, 4(1).

- Ussar, S., Karlina, R., Lutter, D., Miok, V., Fischer, D., Altun, I., Schöttl, T., Schorpp, K., Israel, A., Cero, Cheryl., Johnson, J.W., Kapser-Fischer, I., Böttcher, A., Keipert, S., Feuchtinger, A., Graf, E., Strom, T., Walch, A., Lickert, H., Walzthoeni, T., Heinig, M., Theis, F.J., García-Cáceres, C., Cypess, A.M., Ussar, S. (2021). "[Identification and characterization of distinct murine brown adipocyte lineages](https://www.biorxiv.org/content/10.1101/2020.08.24.264416v1.abstract)". *bioRxiv*, 4(1).

- Yan, X., Miok, V., Karlina, R., Böttcher, A., Lutter, D., Lickert, H., Garcia Caceres, C., Ussar, S. (2023) "[ScRNAseq–Based Analysis of the Stromal Immune Cell Composition in Murine Preweaning and Adult White Adipose Tissues](https://diabetesjournals.org/diabetes/article/72/Supplement_1/1666-P/150177)", *Diabetes* 72 (Supplement_1), 1666-1666.

