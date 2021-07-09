# Masterthesis

This thesis was written as part of the Bioinformatik (M. Sc.) program at the Friedrich-Schiller-Universität Jena. The whole project was supervised by Prof. Dr. Manja Marz and Kevin Lamkiewicz.

## Abstract

Reoccurring local outbreaks of new, highly pathogenic strains of the *Influenza A Virus* (IAV), picture a unnoticed but still persisting major danger to the whole human population, that reached global extend with millions of deaths several times in the past. Due to the lack of a real cure, resort to vaccines producing varying levels of immunization with yearly expiration is inevitable. For better preparation on possible future pandemics, enlarging the knowledge of the IAV is crucial. High evolution-rates by more drastic mutation mechanisms of the IAV, with an aged classification containing little insight, complicate accurate novel research though. This thesis, thus, serves the elaboration of a pipeline usable on all existing and subsequently sequenced genomes of IAV, to not only renew the classification but also keep it up-to-date. Instead of being alignment based, this method utilizes the better scalablity of *k*-mer frequency vectors in a novel hybrid clustering approach, connecting hierarchical with density-based clustering. Most appropriate parameters were thoroughly tested and selected by different validation techniques. For best preservation of the high amount of information included in the vectors used in the clustering, different tools were tested for the best representation in a clusterable dimension. Thereby, a workflow combining a careful selected vector clustering method with appropriate parameters, posterior to a dimension reduction, preserving a high amount of information, was proposed. 

## Repository

The repository contains the following folders and files. The main results are contained in Results.

| Content | Description |
| -- | -- |
| Graphics/ | Folder containing all graphics created with `Inkscape` |
| PCA/ | Folder containing all PK and PD method comparison results |
| Results/ | Folder containing the result of every segments final clustering |
| Thesis/ | Folder containing the components of written the thesis |
| UMAP/ | Folder containing all UK and UD method comparison results |
| Clustering.py | The *Influenza A Virus* clustering tool as described above |
| Environment.yml | The configuration file for recreation of the used environment |
| PCA.ipynb | The pipeline used for the PK and PD method comparison results |
| Thesis.pdf | The thesis in PDF format |
| UMAP.ipynb | The pipeline used for the UK and UD method comparison results |

## Installation

The tool used in the project can be installed and used as described in the following. 

`git clone https://github.com/ahenoch/Masterthesis.git`

`conda env create -f Environment.yml`

`conda activate Masterthesis`

For the default execution discussed in the thesis, the default execution of the tool with the addition of the `-p 50` parameter has to be used. Every FASTA containing sequenced genomes of the *Influenza A Virus* can be used when containing the appropriate header. The `-t 12` can be exchanged for a given number of threads to use.

`python3 Clustering.py -i A.fasta -p 50 -t 12`

The FASTA file exceeds the size allowed to be placed in the repository. The following table can be used as input in the [nucleotide search interface](https://www.fludb.org/brc/influenza_sequence_search_segment_display.spg?method=ShowCleanSearch&decorator=influenza) of the Influenza Research Database (IRD) to recreate the file in a more recent version. The FASTA can be also downloaded directly from the [here](https://cloud.uni-jena.de/s/fYkQ2NAwjND8oEM) in the version of GenBank Genome Sequence/Annotation Update <= 11/2020, that is used in the project.

| Field | Parameter |
| -- | -- |
| Data Type | Genome Segments |
| Virus Type | A |
| Complete Genome | Complete Genome Only |
| Select Segments | All |
| Complete | All |

The header of the FASTA has to be modified on there according to 

1. accession 
2. strain 
3. segment 
4. protein 
5. genus 
6. subtype 
7. date 
8. host 
9. curation

to be used with the expected result. The tool combines three pipeline generated in the project:

- ![Preprocessing Pipeline](/Graphics/Vectorization.pdf)
- ![Clustering Pipeline](/Graphics/Clustering.pdf)
- ![Postprocessing Pipeline](/Graphics/Tree.pdf)

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
| -r | --render | format of output graphics (default: pdf) |
| -t | --threads | number of threads to use for the pairwise cluster validation (default: 12)|
| -h | --help | open the help page |

## Results 

The main results of the project include method comparisons to find the method most appropriate to cluster the *Influenza A Virus* as described in the ![thesis](/Thesis.pdf). A new classification was build using the method in the final step of the project, that includes cluster assignment for all accessions and can be found ![here](/Results/cluster.csv). All clusters are visualized as trees and validation plots.

| Segment | Clustertree | Validation |
| -- | -- | -- |
| 1 | ![here](/Results/Clustertree_Segment_1.pdf) | ![here](/Results/Cluster_Difference_Segment_1.pdf) |
| 2 | ![here](/Results/Clustertree_Segment_2.pdf) | ![here](/Results/Cluster_Difference_Segment_2.pdf) |
| 3 | ![here](/Results/Clustertree_Segment_3.pdf) | ![here](/Results/Cluster_Difference_Segment_3.pdf) |
| 4 | ![here](/Results/Clustertree_Segment_4.pdf) | ![here](/Results/Cluster_Difference_Segment_4.pdf) |
| 5 | ![here](/Results/Clustertree_Segment_5.pdf) | ![here](/Results/Cluster_Difference_Segment_5.pdf) |
| 6 | ![here](/Results/Clustertree_Segment_6.pdf) | ![here](/Results/Cluster_Difference_Segment_6.pdf) |
| 7 | ![here](/Results/Clustertree_Segment_7.pdf) | ![here](/Results/Cluster_Difference_Segment_7.pdf) |
| 8 | ![here](/Results/Clustertree_Segment_8.pdf) | ![here](/Results/Cluster_Difference_Segment_8.pdf) |

## Next Steps

All future steps are described in the thesis. The most important ones are listed here:

- k-mer frequency vectors considering evolutionary aspects with the BLOSUM matrix
- Kneedle Algorithm without area and poly calculation to obliterate even the slightest bias
- Deeper examination of the best dimension to reduce the vectors with `PCA` to
- Commenting the Clustering.py script
 
---

Special thanks to Kevin Lamkiewicz for supervision!