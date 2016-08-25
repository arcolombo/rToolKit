# rToolKit
This is an rToolKit but covertly means repeatToolKit.

this package is for analyzing patient specific kallisto experiments downstream from arkas.
aml is a full kallistoExperiment which has 6 samples that are patient specific, patientPlot can show a specific patient ID and explore repeat biotypes and clonal stages.
```
library(repeatToolKit)
data("aml")

patientPlot(aml,patientID="SU353")

```
