library(rjson)
#AGP 500 metagenome submission
#Run this on the commandline to generate a list of all mzXMLs in a massive ID
#lftp ftp://massive.ucsd.edu -e "cd MSV000080179/peak/mzXMLs/; ls; exit;" |cut -f 21 -d " ">all_AGP_metabolome_samples.txt

#Read the file and chop the identifier up so it matches the sequencing data
mb <- read.table("~/Downloads/all_AGP_metabolome_samples.txt", stringsAsFactors = F)
mb$ID <- sapply(mb$V1, function(x) strsplit(x, "_")[[1]][1])
mg$ID <- sub("^0+", "", sub("10317\\.","", mg$sample_name))

#Read the Qiita ENA ID list
mg <- read.table("~/Downloads/all_metagenome_ena_ids.txt", sep = "\t", fill = T, stringsAsFactors = F, header = 1)

#Merge the tables
mgmg <- merge(mb, mg, by= "ID")
#Make a link out of the mzxml file
mgmg$Metabo_URL <- paste0("ftp://massive.ucsd.edu/MSV000080179/peak/mzXMLs/", mgmg$V1)
#Exclude one duplicate sample
mgmg <- subset(mgmg, ID != 28742)

#Open the json file from the IOMEGA website
mgj <- fromJSON(file = "~/Downloads/test2.json")

#Loop over the list of genomes to populate the genomes
for (i in 1:nrow(mgmg)){
  mgj$genomes[[i]] <- list(genome_ID = list(genome_type = "metagenome",ENA_NCBI_accession = mgmg$experiment_accession[i] ),genome_label = mgmg$ID[i], BioSample_accession = mgmg$experiment_accession[i],
                           publications = "29795809" )

}

#Loop over the table again to populate the list of genome-metabolome links. All labels (sample, extraction etc.) must match what was put in the form
for (i in 1:nrow(mgmg)){
  mgj$genome_metabolome_links[[i]] <- list(genome_label = mgmg$ID[i],
                                           metabolomics_file = mgmg$Metabo_URL[i],
                                           sample_preparation_label = "Metagenome",
                                           extraction_method_label = "Ethanol extraction",
                                           instrumentation_method_label = "TOF")
}
#Write json output to be uploaded to IOMEGA
write(toJSON(mgj), "~/Downloads/agp_500_metagenome_submission.json")
