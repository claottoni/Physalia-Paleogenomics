#########
R session
#########


In this hands-on we will use R to run analyses and create charts from the abundance tables that we generated with Bracken.
We will a sequence dataset of ancient dental calculus from a recent study (Velsko et al, 2019), which consists of paired-end sequence data. 
After filtering and collapsing the reads, classify the reads of each sample with Kraken2 and estimate the abudances with Bracken. 
After that, run the following R script to merge the species aundances of all the samples in one table. The script needs as argument the path to the folder containing the bracken results. 
Move to the bracken results folder and type the following command: 
::
  
  Rscript brackenToAbundanceTable.R . 

The script will generate two abundance tables, ``taxa_abundance_bracken_IDs.txt``, which contains the species names as NCBI IDs, and ``taxa_abundance_bracken_names.txt``, which contains the actual species names. 

Since the Minikraken database was built from complete Bacterial, Archaeal and Viral genomes, we must make a normalization for the genome leghts of each species. To do that will use a Python script that takes three arguments: 

1) The abundance table
2) The table of genome legths
3) The name of the output file

To run the normalization type the command: 
::

  python gL-normalizer-lite.py taxa_abundance_bracken_names.txt prokaryotes_viruses_organelles.table taxa_abundance_bracken_names_normalized.txt

Then we can download the normalized table in your local machine (e.g. with ´´scp´´) and open R. 

.. highlight:: r

First of all, in R, we must set up the folder which contains the abundance table:   
::

  setwd("/path/to/your/folder")
    
Then we can import the abundance table files, setting the 1st column as row names. Then we remove the 1st, 2nd and 3rd column, which will not be used for the analysis. 
::

  table.species = read.delim("taxa_abundance_bracken_names_normalized.txt", header=T, fill=T, row.names=1, sep="\t")
  table.species.final = table.species[,-c(1:3)]

Then we filter the species in the table for their abundance by removing those that are represented below a threshold of 0.02%. 
First we define a function (that we call low.count.removal)
::

  low.count.removal = function(
                          data, 		# OTU count data frame of size n (sample) x p (OTU)
                          percent=0.02	# cutoff chosen
                          ){
      keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
      data.filter = data[,keep.otu]
      return(list(data.filter = data.filter, keep.otu = keep.otu))
  }

We run the function on our table, setting up the threshold at 0.02%: 
::

  result.filter = low.count.removal(t(table.species.final), percent=0.02)

We generate a table with the filtered data: 
::

  table.species.final.flt = result.filter$data.filter
  
In the next step, we normalize the reads for the sequencing depth. This means that we will account for the total reads generated for each sample, and normalize the species abundance for that number. 
To do that we will define a **Total Sum Squared** function 
::

  TSS.divide = function(x){
   x/sum(x)
  }

The function is applied to the table, and each row must represent a sample. For this reason we transpose the table.
::

  table.species.final.flt.tss = t(apply(table.species.final.flt, 1, TSS.divide))

We have just generated a table of species abundances of ancient dental calculus samples, normalized for genome lenghts and sequencing depths. 
We can now include in our analysis a dataset of normalized species abundances generated with Kraken2 representing other microbiomes. You can find this table in the server, download it and import it in R.  
::
 
  microbiomes.literature = read.table("microbiomes_literature.txt", header=T, fill=T, row.names=1, sep="\t")

The following commands are used to merge the two tables, the dental calculus dataset that you generated and the other microbiomes dataset.
::

  table.total = merge(t(table.species.final.flt.tss), t(microbiomes.literature), by=0, all=TRUE)
  table.total[is.na(table.total)] <- 0
  row.names(table.total) = table.total[,1]
  table.total.final = t(table.total[,-1])


*****
UPGMA
*****

Once generated the final including both datasets (dental calculus and other microbiomes), we run an UPGMA cluster analysis. We must first install the ``vegan`` and `ape`` package in R.
::

  install.packages("vegan")
  install.packages("ape")
  library(vegan)
  library(ape)

Then we use vegan to calculate the **Bray-Curtis** distances, and run the cluster analysis.
::

  bray_dist = vegdist(table.total.final, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)
  bray_dist.clust = hclust(bray_dist, method="average", members = NULL)

Finally, we plot the dendrogram: 
::

  plot(as.phylo(bray_dist.clust), type = "unrooted", cex = 0.5, lab4ut="axial", no.margin=T, show.tip.label=T, label.offset=0.02, edge.color = "gray", edge.width = 1, edge.lty = 1)

To visualize better our samples, we can define colors. We will assign group labels on each sample: 
::

  labels = c("Velsko-ancient","Velsko-ancient","Velsko-ancient","Velsko-ancient","Velsko-ancient","Velsko-ancient","Velsko-ancient","Velsko-ancient","Velsko-ancient","Velsko-ancient",
						"Velsko-modern","Velsko-modern","Velsko-modern","Velsko-modern","Velsko-modern","Velsko-modern","Velsko-modern","Velsko-modern",
						"Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth",
						"Soil","Soil","Soil","Soil","Soil","Soil","Soil",
						"Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth",
						"Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque",
						"Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque",
						"Plaque","Plaque","Plaque","Plaque","Plaque",
						"Skin","Skin","Skin","Skin",
						"Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut",
						"Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut",
						"Gut","Gut","Gut","Gut",
						"Skin","Skin","Skin","Skin","Skin",
						"Plaque","Plaque")

To have a better look at the correspondence of data we can create a dataframe: 
::

  table.total.final.df = as.data.frame(table.total.final)
  table.total.final.df$group = labels
  
We assign colors to each label: 
::

  coul=c("#E41A1C",		#Ancient calculus		
		"#419681",		#Ancient tooth					
		"#4DAF4A",		#Gut				
		"lightgray",	#Plaque	
		"#984EA3",		#Skin		
		"#FF7F00",		#Soil		
		"goldenrod",	#Velsko-ancient		
		"#994C00")		#Velsko-modern

And finally, we plot again the dendrogram, this time by customizing the tips assigning color-coded labels:
::

  plot(as.phylo(bray_dist.clust), type = "unrooted", cex = 0.5, lab4ut="axial", no.margin=T, show.tip.label=T, label.offset=0.02, edge.color = "gray", edge.width = 1, edge.lty = 1)
  tiplabels(pch=19, col = coul[factor(labels)], bg = coul[factor(labels)], cex=1, lwd=1)          

And we can add a legend:
::

  legend("topleft", legend = sort(unique(labels)), bty = "n", col = coul, pch = 19, pt.cex=1, cex=0.6, pt.lwd=1)


