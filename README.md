# NIA CARD ANG CAVIAR QTL Fine Mapping Pipeline
This repository includes a pipeline for streamlined fine mapping of potentially causal QTL variants using [CAVIAR](http://genetics.cs.ucla.edu/caviar/index.html).

The pipeline is run in four steps:
1. Preparation of up to 100 most significantly associated variants per phenotype, making sure to include at least one structural variant when necessary, for table of phenotype/variant pairs passing FDR-corrected p-value (q-value) filter (typically 0.05).
2. Calculating linkage disequilibrium (LD) matrices for each variant set selected above.
3. Running CAVIAR on each variant set per phenotype.
4. Appending CAVIAR post causal probabilities to initial significant phenotype/variant pair table.
