# Convert G4 hunter output from the website into BED files in the genomic coordinates

This repository include the script and demodata to convert G4hunter output results into bed file in the genomic coordinate.

There are a few parameters that require input file locations and tune input/output file format.

```
g4_outfile:                 Path to the output files from G4 hunter
region_file:                Path to the region file, can be your peak calling results or the gene region file. First column is chromosoeom, second column is region starting site, third column is region ending site. The rest columns will not be used so can be nothing or anything.
separate_by:                The separator of region file, can be "\t" or " " or ","
first_position_ignore:      TRUE or FALSE to indicate whether the starting position in the region starting site is included in your sequence (fasta) file. TRUE for not included, FALSE for included.
out_file:                   Path to store the transformed output in BED format.
output_include_origin_peak: TRUE or FALSE to indicate whether you want to save the original region information. FALSE to only keep the motif information.
```

Demon runs: In terminal or in a shell script.

```
g4_outfile="....../demoData/demo_g4hunter_output.csv"
region_file <- "....../demoData/AllgenesChrStartEndNameStrand.csv"
separate_by <- ","
first_position_ignore <- FALSE
out_file <- "....../path/to/the/output/folder/demo_g4huter_motif.bed"
output_include_origin_peak <- FALSE

Rscript convert_G4hunter_out_to_bedfiles.R "${g4_outfile}" "${region_file}" "${separate_by}" "${first_position_ignore}" "${out_file}" "${output_include_origin_peak}"

```
