## data information

This notebook loads the expression array and methylation array dataset and match for the patient identifier.

The resulting data were written as tsv in `data/Figueroa/matchedata`.

## preprocessing

This notebook preprocessed the two matched dataset as follows
- expression array was filtered for probes that has a gene annotated and participates in biological process related to transcription and methylation
- methylation array was successively normalized (both across rows and columns)

The resulting data are written as tsv in `data/Figueroa/processeddata`


## PhenoGraph clustering

This notebook performs PhenoGraph clustering on the successively normalized data of methylation array and write the community assign to `data/Figueroa/clusters`


## Eigen methylation

This notebook computes the eigengene methylation profile of each community identified by PhenoGraph clustering. Example Figuere of the module trajectories and the eigengene methylation profile can be found in `figures/`