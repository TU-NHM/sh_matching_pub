SH matching analysis

Description

Developed as part of EOSC-Nordic (https://www.eosc-nordic.eu/) project (task 5.2.1: Cross-border data processing workflows), UNITE SH matching analysis is a digital service for the global species discovery from eDNA (environmental DNA). SH matching service is based on the UNITE (https://unite.ut.ee) datasets hosted in PlutoF (https://plutof.ut.ee). Its output includes information about what species are present in eDNA samples, are they potentially undescribed new species, where are they found in other studies, are they alien or threatened species, etc. The output will provide DOI (Digital Object Identifier) based stable identifiers for the communicating species found in eDNA. DOIs are connected to the taxonomic backbone of PlutoF and GBIF (https://www.gbif.org). In this way every DOI is accompanied by a taxon name which is still widely used for the communication of species. In the case of undescribed species, DOIs will soon be issued by the PlutoF system (only if SH matching service integrated with the PlutoF platform is used for the analysis). SH matching service covers all Eukaryota by using rDNA ITS marker sequences accompanied by sample metadata.


Citing

When using this resource, please cite as:

Abarenkov K, Kõljalg U, Nilsson RH (2022) UNITE Species Hypotheses Matching Analysis. Biodiversity Information Science and Standards 6: e93856. https://doi.org/10.3897/biss.6.93856


SH matching output

The output archive includes the following files:

1. matches/ directory - a folder containing the main output files describing matches between the user's dataset and existing SHs, as well as newly formed compound clusters and SHs, and an HTML summary.
    a. matches_out_[distance_threshold].csv - a data matrix describing the user's sequences assigned to existing SHs or forming new SHs.
    b. matches_out_[distance_threshold].html - an HTML summary of the data contained in the corresponding CSV files.
    c. matches_out_all.txt - a data matrix describing the user's sequences assigned to existing SHs or forming new SHs, with all distance thresholds merged.
    d. maches_out_taxonomy.csv - a list of query sequences with consensus taxonomy, taking into account the taxonomy of SHs across all distance thresholds.
2. err_[run_id].log - information and error messages from different analysis steps, mainly intended for debugging and support;
3. excluded_[run_id].txt - a list of sequences excluded from the analysis, along with the reasons for their exclusion;
4. source_[run_id]_names - a translation table between the user's sequence identifiers and the internal system identifiers used throughout the analysis;
5. source_[run_id]_fastanames - a list of duplicate sequences found in the user's input data;
6. krona_03.html - a Krona HTML file that can be interactively visualized in a web browser. The pie chart shows all taxonomic units for which taxonomic information is available from matching SHs and compound clusters.
