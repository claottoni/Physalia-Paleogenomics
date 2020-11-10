#########
R session
#########

.. highlight:: r

In this hands-on we will use R to run analyses and create charts from the abundance tables that we generated with Kraken2.

First of all, in R, we must set up the folder which contains the abundance table:   

::

  setwd("~/Documents/WORK/2-PROJECTS/HiddenFoods/HiddenFoods_analysis/krk2_genome_length_normalized_abundances_July2019/species")
  
Then we can import the abundance table files: 
::

  # --> you can edit the files, eg remove last column 'Detail!
  # dataframes in R must do not allow identical row names, force row names as numbering by using row.names=NULL.
  #EOF issue with single quotes (') in species names, use read.delim
  HiddenFoods_Jan2019.krk2.norm.species = read.delim("HiddenFoods_Jan2019.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  HiddenFoods_May2019.krk2.norm.species = read.delim("HiddenFoods_May2019.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Warinner2014.krk2.norm.species = read.delim("Warinner2014.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Velsko_modern_calculus.krk2.norm.species = read.delim("Velsko_modern_calculus.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Velsko_ancient_calculus_total.krk2.norm.species = read.delim("Velsko_ancient_calculus_total.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Plaque.krk2.norm.species = read.delim("Plaque.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Mann2018.krk2.norm.species = read.delim("Mann2018.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Baboons_Ozga.krk2.norm.species = read.delim("Baboons_Ozga.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Weyrich2017.krk2.norm.species = read.delim("Weyrich2017.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Baboons_Egypt_all.krk2.norm.species = read.delim("Baboons_Egypt_all.krk2.rcf.abundance.species.flt10.teeth.env.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Baboons_Egypt_nofilter.krk2.norm.species = read.delim("Baboons_Egypt_nofilter.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Baboons_Egypt_teeth.krk2.norm.species = read.delim("Baboons_Egypt_teeth.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Skin.krk2.norm.species = read.delim("Skin.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  ObregonTito_gut.krk2.norm.species = read.delim("ObregonTito_gut.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Soil.krk2.norm.species = read.delim("Soil.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Chimps.krk2.norm.species = read.delim("Chimps_Ozga.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Brealey_animals.krk2.norm.species = read.delim("Brealey_animals.krk2.rcf.abundance.species.norm2.final", header=T, fill=T, row.names=NULL, sep="\t")
  Brealey_animals.krk2.norm.species = Brealey_animals.krk2.norm.species[,-c(2:5,7,8)]		#remove Blanks from filtrated Brealey, Gb1reg/flt, Gb2reg/flt
  Eisenhofer_Japan_flt.krk2.norm.species = read.table("Eisenhofer_Japan.krk2.rcf.abundance.species.norm2.final.flt", header=T, fill=T, row.names=NULL, sep="\t")
  Eisenhofer_Japan_flt.krk2.norm.species = Eisenhofer_Japan_flt.krk2.norm.species[,-c(19,20)]	#remove EBC-2 in  Eisenhofer which has no classified species, and Details column.




*******
Barplot
*******


*****
UPGMA
*****


****
nMDS
****
