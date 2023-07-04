# sonsohdezz_TFG
Repository containing files and scripts created for Final Degree Project in the Biotechnology Grade at the Technical University of Madrid (UPM) in july of 2023.

<details>

<summary> Final Degree Project Information </summary>

\
* **Title:** "EVALUATION OF HUMAN GUT MICROBIOME TAXONOMIC PROFILES DERIVED FROM THE EGGNOG-MAPPER FUNCTIONAL ANNOTATION TOOL."
* **Author:** Sonsoles Hernández Piñel
* **Tutors:** Carlos Pérez Cantalapiedra & Joaquín Giner Lamia
* **Institution:** [Technical University of Madrid (UPM)](https://www.upm.es/internacional)
* **Collaborating Institution:** [Centre for Biotechnology and Plant Genomics (CBGP)]([https://pages.github.com/](https://www.cbgp.upm.es/index.php/en/about-us))

  
</details>


| File  | Description |
| ------------- | ------------- |
| `contigs.py`  | Genome fragmentation script, which uses as parameters the desired contig size (in bp) and the desired displacement window size (in bp), and returns a fasta file with the contigs created.  |
| `eggnog_phyloseq.py`  | Transformation script that creates the three required input files in phyloseq from the taxonomic annotation file obtained with eggNOG-mapper.  |
| `lineage.py`  | Script that generates a table with the taxonomies of a metagenomic sample from the functional annotation file obtained with eggNOG-mapper.|
| `motus_phyloseq.py`  | Transformation script that creates the three required input files in phyloseq from the taxonomic annotation file obtained with mOTUs v3.|
| `phyloseq.R` | Script that generates the graphs of interest for the analysis of taxonomic results obtained with functional analysis tools, such as eggNOG-mapper or mOTUs v3. |
 
