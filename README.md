# Masterthesis

This masterthesis was written as part of the Bioinformatik (M. Sc.) program at the Friedrich-Schiller-Universität Jena. The project was supervised by Prof. Dr. Manja Marz and Kevin Lamkiewicz.

## Introduction



## Repository

The repository contains the following folders and files. The main results are contained in the Results folder.

| Content | Description |
| -- | -- |
| Graphics/ | Folder containing all graphics created with `Inkscape` |
| PCA/ | Folder containing all PK and PD method comparison results |
| Results/ | Folder containing the result of every segments final clustering |
| Thesis/ | Folder containing the components of written the thesis |
| UMAP/ | Folder containing all UK and UD method comparison results |
| A.fasta | The FASTA file used in the project |
| Clustering.py | The *Influenza A Virus* clustering tool as described above |
| Environment.yml | The configuration file for recreation of the used environment |
| PCA.ipynb | The pipeline used for the PK and PD method comparison results |
| Thesis.pdf | The thesis in PDF format |
| UMAP.ipynb | The pipeline used for the UK and UD method comparison results |

## Installation

The tool used in the project can be installed as described in the following. 

`git clone https://github.com/ahenoch/Masterthesis.git`

`conda env create -f Environment.yml`

`conda activate Masterthesis`

For the default execution discussed in the thesis, the default execution of the tool with the addition of the `-p 50` parameter has to be used. Every FASTA containing sequenced genomes of the *Influenza A Virus* can be used when containing the appropriate header. The `-t 12` can be exchanged for a given number of threads to use.

`python3 Clustering.py -i A.fasta -p 50 -t 12`

The tool combines three pipeline generated in the project:

- Vectorization Pipeline, that can be found ![here](/Graphics/Vectorization.pdf)
- Clustering Pipeline, that can be found ![here](/Graphics/Clustering.pdf)
- Visualization Pipeline, that can be found ![here](/Graphics/Tree.pdf)

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

## Results 

The main results of the project include method comparisons to find the method most appropriate to cluster the *Influenza A Virus* as described in the ![thesis](/Masterthesis.pdf). A new classification was build using the method in the final step of the project, that can be found ![here](/Results/cluster.csv). All clusters are visualized as trees and validation plots.

| Segment | Clustertree | Validation |
| -- | -- | -- |
| 1 | ![here](/Graphics/Clustertree_Segment_1.pdf) | ![here](/Graphics/Cluster_Difference_Segment_1.pdf) |
| 2 | ![here](/Graphics/Clustertree_Segment_2.pdf) | ![here](/Graphics/Cluster_Difference_Segment_2.pdf) |
| 3 | ![here](/Graphics/Clustertree_Segment_3.pdf) | ![here](/Graphics/Cluster_Difference_Segment_3.pdf) |
| 4 | ![here](/Graphics/Clustertree_Segment_4.pdf) | ![here](/Graphics/Cluster_Difference_Segment_4.pdf) |
| 5 | ![here](/Graphics/Clustertree_Segment_5.pdf) | ![here](/Graphics/Cluster_Difference_Segment_5.pdf) |
| 6 | ![here](/Graphics/Clustertree_Segment_6.pdf) | ![here](/Graphics/Cluster_Difference_Segment_6.pdf) |
| 7 | ![here](/Graphics/Clustertree_Segment_7.pdf) | ![here](/Graphics/Cluster_Difference_Segment_7.pdf) |
| 8 | ![here](/Graphics/Clustertree_Segment_8.pdf) | ![here](/Graphics/Cluster_Difference_Segment_8.pdf) |