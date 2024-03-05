library("dplyr")

##Import CP-RSEM results from each sample pair for genes
CPgenes1 <- read.table("RSEM.genes.results1.tsv", header=T);
CPgenes2 <- read.table("RSEM.genes.results2.tsv", header=T);
CPgenes3 <- read.table("RSEM.genes.results3.tsv", header=T);

#isoforms
CPisoforms1 <- read.table("RSEM.isoforms.results1.tsv", header=T);
CPisoforms2 <- read.table("RSEM.isoforms.results2.tsv", header=T);
CPisoforms3 <- read.table("RSEM.isoforms.results3.tsv", header=T);


##generating narrow down tables if required
CPgenes1Narrow <- CPgenes1[c("gene_id","transcript_id.s.","expected_count", "TPM")];
CPgenes2Narrow <- CPgenes2[c("gene_id","transcript_id.s.","expected_count", "TPM")];
CPgenes3Narrow <- CPgenes3[c("gene_id","transcript_id.s.","expected_count", "TPM")];

#isoforms
CPisoforms1Narrow <- CPisoforms1[c("transcript_id","gene_id", "expected_count", "TPM")];
CPisoforms2Narrow <- CPisoforms2[c("transcript_id","gene_id", "expected_count", "TPM")];
CPisoforms3Narrow <- CPisoforms3[c("transcript_id","gene_id", "expected_count", "TPM")];

##fused genes dataframe for ease of use (assumed files have the same transcript order)
CPfused <- CPgenes1Narrow[c("gene_id","transcript_id.s.")];
CPfused["TPM_CP1"] <- CPgenes1Narrow$TPM;
CPfused["TPM_CP2"] <- CPgenes2Narrow$TPM;
CPfused["TPM_CP3"] <- CPgenes3Narrow$TPM;
CPfused["Expected_count1"] <-CPgenes1Narrow$expected_count;
CPfused["Expected_count2"] <-CPgenes2Narrow$expected_count;
CPfused["Expected_count3"] <-CPgenes3Narrow$expected_count;

##fused isoforms dataframe for ease of use (assumed files have the same transcript order)
CPfusedIso <- CPisoforms1Narrow[c("transcript_id","gene_id")];
CPfusedIso["TPM_CP1"] <- CPisoforms1Narrow$TPM;
CPfusedIso["TPM_CP2"] <- CPisoforms2Narrow$TPM;
CPfusedIso["TPM_CP3"] <- CPisoforms3Narrow$TPM;
CPfusedIso["Expected_count1"] <-CPisoforms1Narrow$expected_count;
CPfusedIso["Expected_count2"] <-CPisoforms2Narrow$expected_count;
CPfusedIso["Expected_count3"] <-CPisoforms3Narrow$expected_count;

#selected subset specifically of only the genes identified in the selected dataset
idListDf12 <- read.table("AllIds.allids");
idList1 <- idListDf12[[1]];
row.names(CPfused) <- CPfused$gene_id;
indexesGen <- subset(CPfused, rownames(CPfused) %in% idList1);

##writing the datatable to a CSV for analysis and manual additions
write.csv(indexesGen, "/home/h/Documents/SilicoProject/SilicoR/GeneTable.csv", row.names = FALSE)






