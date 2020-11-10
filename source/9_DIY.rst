##############
DO-IT-YOURSELF
##############

In this final hands-on session you will analyse shotgun (reduced) sequencing data generated from an ancient human tooth.The genomic library built from the DNA extract was sequenced on an Illumina platform in paired-end mode. Your task is:   

1. Process the raw reads (remove adapters, merge the reads, section 2). 
2. Align the reads to the human mitochondrial DNA (mtDNA) reference sequence, assess the damage of DNA molecules, call the variants (sections 4-5-6).  
3. Run the metagenomic screning of the DNA extract with Kraken using the Minikraken database (section 3).

After reads pre-processing it is up to you whether first aliging the reads or screening the metagenomic content. 

- **Option 1**: You can use your ``vcf`` file to assign an haplogroup to the human samples that you analysed. Some useful tools for haplogroup assignation:  
  
    - Check the variant positions in Phylotree (http://www.phylotree.org/)  
    - Load the ``vcf`` file in Haplogrep (https://haplogrep.uibk.ac.at)
   
- **Option 2**: Run again the metagenomic screening with a Custom Database of Kraken (provided by us), and compare the results with those obtained with Minikraken.