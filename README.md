# Masterthesis

This masterthesis was written as part of the Bioinformatik (M. Sc.) program at the Friedrich-Schiller-Universität Jena. The project was supervised by Prof. Dr. Manja Marz and Kevin Lamkiewicz.

## Introduction



## Repository

The repository contains the following folders and files.

| Content | Description |
| -- | -- |
| Graphics/ | Folder containing all graphics created with `Inkscape` |
| PCA/ | Folder containing all PK and PD method comparison results |
| Results/ | Folder containing the result of every segments final clustering |
| Thesis/ | Folder containing the components of written the thesis |
| UMAP/ | Folder containing all UK and UD method comparison results |
| A.fasta | The FASTA file used in the project |
| Clustering.py | The *Influenza A Virus* clustering tool |
| Environment.yml | The configuration file for recreation of the used environment |
| PCA.ipynb | The pipeline used for the PK and PD method comparison results |
| README.md | The instructions for usage of the clustering tool |
| Thesis.pdf | The thesis in PDF format |
| UMAP.ipynb | The pipeline used for the UK and UD method comparison results |

## Installation

`git clone https://github.com/ahenoch/Masterthesis.git`

`conda env create -f Environment.yml`

`conda activate Masterthesis`

For the default execution discussed in the thesis, the default execution of the tool with the addition of the `-p 50` parameter has to be used. Every FASTA containing sequenced genomes of the *Influenza A Virus* can be used when containing the appropriate header. The `-t 12` can be exchanged for a given number of threads to use.

`python3 Clustering.py -i A.fasta -p 50 -t 12`

![Vectorization Pipeline](/Graphics/Vectorization.pdf)

![Clustering Pipeline](/Graphics/Clustering.pdf)

![Visualization Pipeline](/Graphics/Tree.pdf)

## Manual

| Parameter | | Description |
| -- | -- | -- |
| -i | --infile | path to input file (e. g. A.fasta) |
| -o | --outfolder | path to input file (default: Results) |
| -s | --segments | segments to run the pipeline on (default: 1 2 3 4 5 6 7 8) |
| -c | --custom_header | if using FASTA with a custom header, every part of it has to be declared (default: accession strain segment protein genus subtype date host curation genome) |
| -m | --metric | metric to use in the pipeline (default: cosine) |
| -mc | --min\_clust | min\_cluster_size parameter for HDBSCAN (default: 2) |
| -ms | --min\_sample | min\_samples parameter for HDBSCAN (default: 1) |
| -n | --neigbors | n\_neighbors parameter for UMAP (default: 100) |
| -u | --umap | n\_components parameter for UMAP (optional) |
| -p | --pca | n\_components parameter for PCA (optional) |
| -k | --max\_kneedle | search area for Kneedle Algorithm (default: 500) |
| -r | --render | format of output graphics (default pdf) |
| -t | --threads | number of threads to use for the pairwise cluster validation |
| -h | --help | open the help page |