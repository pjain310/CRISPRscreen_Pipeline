# Methods Comparison Pipeline for pooled CRISPR screens

This methods comparison pipeline was created to determine the best tools/set of tools to be used to analyze pooled CRISR screen data. This pipeline runs different screen anlaysis tools (namely, Mageck, MageckMLE, CB2, PBNPA) on counts data, and also makes initial counts distribution and correlation plots for the data.


### Prerequisites

This pipeline has the following software requirements:

1. Python version 3.x
2. R version 3.5 or higher
3. R packages:
   3.1 PBNPA
   3.2 CB2
   3.3 tibble
   3.4 magrittr
   3.5 dplyr
   3.6 ggplot2
   3.7 ineq
   3.8 VennDiagram
4. Mageck version 0.5+


### Running the script

The pipeline can be run by writing the following command

```
./pipeline.sh -i ./sample_input/Evers-shRNA_CRISPR.txt -s ./sample_input/sample_map_Evers.txt -o ./sample_output/results_Evers -t pbnpa

```

## Author

* **Prerna Jain** - (https://github.com/pjain310)
Please feel free to reach out at 08.prerna@gmail.com in case you have any questions!
