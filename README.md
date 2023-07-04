# sonsohdezz_TFG
Repository containing files and scripts created for Final Degree Project in the Biotechnology Grade at the Technical University of Madrid (UPM) in july of 2023.

<details>

<summary>  Display to check more information about the project. </summary>

#### Title
> Evaluation of human gut microbiome taxonomic progiles derived from the eggNOG-mapper functional annotation tool. 
#### Author
> Sonsoles Hernández Piñel
#### Tutors
> - Carlos Pérez Cantalapiedra
> - Joaquín Giner Lamia
#### Institution
> [Technical University of Madrid (UPM)](https://www.upm.es/internacional)
#### Collaborating Institution
> [Centre for Biotechnology and Plant Genomics (CBGP)](https://www.cbgp.upm.es/index.php/en/about-us)

  
</details>



<details>

<summary> Display to check more information about the files and their use. </summary>

#### 

| File  | Description | Use | 
| ------------- | ------------- | ------------- |
| `contigs.py`  | Genome fragmentation script, which uses as parameters the desired contig size (in bp) and the desired displacement window size (in bp), and returns a fasta file with the created contigs. To run it, the path to the file containing the genome to be fragmented is to be introduced, as well as the path to the directory where you want to store the result file. | `python contigs.py -i <inputfile> -o <outputdir> -c <contigs_size> -w <scrollingwindow_size>` |
| `eggnog_phyloseq.py`  | Script that creates the three required input files for Phyloseq from the taxonomic annotation file obtained with eggNOG-mapper. To run it, you must introduce the path to the folder containing the taxonomies obtained for the different samples, and the path to the folder where you want to store the result files.| `python eggnog_phyloseq.py -i <inputdir> -o <outputdir>` |
| `lineage.py`  | Script that generates a table with the taxonomies of the different contigs of a metagenomic sample from the functional annotation file obtained with eggNOG-mapper. To run it, you must introduce the path to the folder containing the .tsv annotation file(s) obtained after running eggNOG-mapper, and the path to the directory where you want to store the result file(s). | `python lineage.py -i <inputdir> -o <outputdir>` |
| `motus_phyloseq.py`  | Transformation script that creates the three required input files for Phyloseq from the taxonomic annotation file obtained with mOTUs v3. To run it, you must introduce the path to the folder containing the taxonomies obtained for the different samples, and the path to the folder where you want to store the result files. | `python motus_phyloseq.py -i <inputdir> -o <outputdir>` |
| `phyloseq.R` | Script that generates the graphs of interest for the analysis of taxonomic results obtained with functional analysis tools, such as eggNOG-mapper or mOTUs v3. | 


</details>


 
