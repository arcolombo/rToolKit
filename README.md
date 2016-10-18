# rToolKit
We are interested in analyzing the repeat patterns in Buenrostro clonal data of pHSC, LSC and Blasts.  


# Are there differences of repeat Expression between stages of pHSC,LSC, and Blast?
  The stage level analysis shows that the differential expression between LSC-pHSC is down regulated with respect to repeats (TMM and FDR controlled).  The composition of repeat biotypes in the DE list for LSC-pHSC consists of ALU 22-29%, ERV1 29-30%, and L1 19-16%.   
 For Blast-LSC we see ERV1 26% ,L1 24% ,and ALU 17%.

 When examining the overall expression of individual patients (TMM and TPM normalized) we see slightly less expression in each individual patient of repeats in LSC overall.  Many LTR groups of the topMADs have high repeat expression in LSC. however in most individual patients LSC are expressed lower than pHSC and Blast, namely in Alu repeats specifically.  
  From the Differential Expression list we observe down regulated expression in repeat elements of LSC compared to pHSC and Blast. 
```
  library(repeatToolKit)
  load("aml.RData")
  amlX<-kexpByStage(aml)
  patientPlot_cpm_library_norm(amlX,patientID="SU353",normType="TMM")
  patientPlot_cpm_library_norm(amlX,patientID="SU070",normType="TMM")
  patientPlot_cpm_library_norm(amlX,patientID="SU209",normType="TMM")
  patientPlot_cpm_library_norm(amlX,patientID="SU209",normType="TMM")
  patientPlot_cpm_library_norm(amlX,patientID="SU444",normType="TMM")
  patientPlot_cpm_library_norm(amlX,patientID="SU496",normType="TMM")
  patientPlot_cpm_library_norm(amlX,patientID="SU575",normType="TMM")
  patientPlot_cpm_library_norm(amlX,patientID="SU583",normType="TMM")
  patientPlot_cpm_library_norm(amlX,patientID="SU654",normType="TMM")
  stageWiseAnalysis(amlX,byLevel="transcript")

```
 We observe many different expressions from LTRs, the TopDE list consists of ERV1, with minor contributions of ERV3,ERVK,ERVL.


# What are the Pathway Activities for pHSC,LSC,and Blasts?
We observe for pHSC-LSC HALLMARK INFLAMMATORY RESPONSE logFC 0.2034681, p.value 1.091457e-02, FDR 2.742192e-02
We observe that for Blast-LSC HALLMARK INFLAMMATORY RESPONSE has logFC 0.50284256 with a pvalue of 2.168942e-02 , FDR 4.016559e-02
 The activity of the inflammation is low for the entire cell types that were not cleared by the immune system (still alive).  we imagine that for dead cells that were cleared, the immune system epsilon near the event of clearance, the inflammatory response would be globally much higher. 
```
  setwd(""/home/arcolombo/Documents/Arkas-Paper-Data/AML-bonemarrow-LSCs/patient-plot-data/wgcna_data/repeat_short_list_biotypes")
 dbn<-"qusageDbLite.cpm.sqlite"
 kexpEnrich(qusageDbLite(dbn),contrast="phsc")
 kexpEnrich(qusageDbLite(dbn),contrast="blast")

```
 
# Samples Clustered By Gene Expression Correlate To Specific Repeat Elements
Here we cluster samples based on gene expression profile hierarchial cluster, andsee that there are positive association of repeat elements of ALU, ERV1, ERV3,ERVL,ERVK, and LTR Retrotransposons to pHSC and Blast as a group.

![Correlation and RE](/inst/images/TxBiotype_Correlation_Samples-1.png)

