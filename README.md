# scRNA-seq of B cells in primary Sjögren's syndrome

***
### Citation

This is a public repository containing scripts used in the publication:

Arvidsson G, Czarnewski P, Johansson A, Raine A, Imgenberg-Kreuz J, Nordlund J, Nordmark G, Syvänen AC. [Multi-modal single cell sequencing of B cells in primary Sjögren's Syndrome](https://pubmed.ncbi.nlm.nih.gov/37610265/). Arthritis Rheumatol. 2023 Aug 23. doi: 10.1002/art.42683. Epub ahead of print. PMID: 37610265.

***
### Conda environment

Two files for setting up the analysis conda environment:
- `conda/environment_pSS.yml`  
- `conda/environment_pSS_MacOS.yml`


### Running the analysis
1. Clone the repository\
```
git clone https://github.com/Molmed/pSS_scRNAseq.git
```

2. Create and activate the conda environment\
```
cd pSS_scRNAseq

conda activate base
conda install -c conda-forge mamba

mamba env create -n env_pSS -f env_pSS_MacOS.yml
conda activate env_pSS
```

***
### Datasets

Datasets used in the study:

| Technology | Dataset | Reference | Accession No |
|------------|---------|-----------|--------------|
| 10X scRNAseq 5' v1.1 | Gene expression of B cells from pSS patients | present study | [GSE214974](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214974) |
| 10X scVDJseq | VDJ sequencing of B cell receptor from pSS patients | present study | [10.17044/scilifelab.21269619]() |
| 10X scRNAseq 5' | Gene expression from sorted B cell subtypes | [Stewart et al. Frontiers in Immunology 2020](https://pubmed.ncbi.nlm.nih.gov/33815362/) | [E-MTAB-9544](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9544/) |
| PCR product sequencing | IGH CDR3 sequences from pSS patients | [Corsiero et al. Plos One 2014](https://pubmed.ncbi.nlm.nih.gov/25535746/) | #### |
| PCR product sequencing | IGH CDR3 sequences from pSS patients | [Visser et al. Frontiers in Immunology 2020](https://pubmed.ncbi.nlm.nih.gov/32760405/) | GenBank Accession Numbers: MT151387–MT151600 |
| Protein sequencing of purified RF antibodies | IGH CDR3 sequences from RF antibodies in pSS patients | [Wang et al. Arthritis Rheumatol. 2018](https://pubmed.ncbi.nlm.nih.gov/29697211/) | #### |
