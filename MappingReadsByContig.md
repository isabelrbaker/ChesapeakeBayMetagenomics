This code will help us determine gene abundance by how reads map to coding sequences on the contigs.

The pros: mapping will be more sensitive.
The cons: as we are mapping to sample-specific contig sequences, we can't compare between samples as easily.

For this, we will need the reads mapped to each contig. Luckily, JGI did this for us already with bbmap. The SAM files can be found here:

```
/home/chesapeake2019metagenomes/ReadsMappedtoContigs/*/QC_and_Genome_Assembly/*/pairedMapped.sam.gz
```

We will be mapping to these files using a tool called featureCounts.


```
ml anaconda
conda activate
conda create -n featurecounts
conda activate featurecounts
conda install bioconda::subread
conda deactivate #this is to exit the conda environment
```

Now, to figure out where our genes are (where we want to see mapped reads), we need to provide our gff file to featurecounts in SAF format. 
NOTE that I had to replace all contig IDs with their local names.

```
files<-sort(list.files("../../scafID_gff", full.names = TRUE,recursive=TRUE)) #read in the gff

for(f in files){
  df<-read.gff(f)
  geneID<-str_extract(df$attributes, "ID=.*")
  geneID<-gsub("ID=","",geneID)
  geneID<-gsub("\\;.*","",geneID)
  
  chr<-df$seqid
  
  start<-df$start
  end<-df$end
  strand<-df$strand
  
  gtf<-data.frame(GeneID=geneID,Chr=chr,Start=start,End=end,Strand=strand)
  
  #now, create new file name
  fn<-paste("../../scafID_SAF/",sub(".gff","",sub(".*/", "", f)),".SAF",sep="")
  
  write.table(gtf,fn,sep="\t",quote=F,row.names = F)
  
  
}

```

Now, all of our SAF files can be found here:

```
/home/chesapeake2019metagenomes/ReadsMappedtoContigs/SAF
```

Now, we need to sort the SAM files into BAM files


```

#!/bin/bash
#
#SBATCH --job-name=samtoolssort        # Job name
#SBATCH --output=output.txt           # Standard output file
#SBATCH --error=error.txt             # Standard error file
#SBATCH --mem=100G
#SBATCH --ntasks=48           # Tasks


cd ~/chesapeake2019metagenomes/



#run samtools command to convert sam files to sorted bam files
#load samtools
ml python/3.11.6 samtools/1.15.1

#run samtools command to convert sam files to sorted bam files

for folder in ReadsMappedtoContigs/Ga*
        do
	name=$(basename "${folder}")

        samtools sort -o "${folder}"/"${name}".sorted.bam -O bam "${folder}"/QC_and_Genome_Assembly/*/pairedMapped.sam.gz

done


```
Ok, now we can finally count reads!

```
#!/bin/bash
#
#SBATCH --job-name=featureCounts      # Job name
#SBATCH --output=output.txt           # Standard output file
#SBATCH --error=error.txt             # Standard error file
#SBATCH --mem=100G
#SBATCH --ntasks=48           # Tasks



cd ~/chesapeake2019metagenomes/

ml anaconda

conda activate featurecounts

for folder in ReadsMappedtoContigs/Ga*
	do

	name=$(basename "${folder}")

        featureCounts "${name}".sorted.bam -a ReadsMappedtoContigs/SAF/"${name}"_annotation_scafID.SAF -B -F SAF -o ReadsMappedtoContigs/Counts/"${name}"_featureCounts.txt -p -P

done

```

Now that we have counted all the reads, we can compile them into one big data frame in R.

```

setwd("..")
path<-"featureCounts_contigs"
files<-sort(list.files(path, pattern=".txt", full.names = TRUE,recursive=FALSE)) #read in the


counts<-data.frame(geneID=character(0),scaf=character(0),start=numeric(0),end=numeric(0),
                   strand=character(0),length=character(0),reads=numeric(0))

for(f in files){
  df<-read.table(f,header=T,stringsAsFactors = F,sep="\t")
  colnames(df)<-c("geneID","scaf","start","end","strand","length","reads")
  counts<-rbind(counts, df)
}

counts$start<-NULL
counts$end<-NULL
counts$strand<-NULL

scaf.data<-scaf.data %>%
  left_join(counts)




#ok, now to convert to fpkm
#fpkm=(#reads mapped to gene * 10^9)/(#total reads mapped to sequences in sample*length of gene)

counts<-scaf.data[!is.na(scaf.data$reads),]

nrow(counts)

countSums<- counts %>%
  group_by(site) %>% 
  summarise(totalCounts = sum(reads))

counts<- counts %>% 
  left_join(countSums)

counts$reads<-as.numeric(counts$reads)
counts$totalCounts<-as.numeric(counts$totalCounts)
counts$length<-as.numeric(counts$length)

counts$fpkm<-(counts$reads*(10^9))/(counts$totalCounts*counts$length)


write.table(counts, 
            "ReadCounts_Conting_FeatureCounts.txt",
            sep="\t",quote=F,row.names = T)


```
