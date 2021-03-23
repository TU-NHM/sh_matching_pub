SH matching analysis

Description

Developed as part of EOSC-Nordic project (task 5.2.1: Cross-border data processing workflows), UNITE SH matching analysis is a digital service for the global species discovery from environmental DNA. SH matching service is based on the UNITE datasets (https://unite.ut.ee) hosted in PlutoF (https://plutof.ut.ee). Its output includes information about what species are present in eDNA samples, are they potentially undescribed new species, where are they found in other studies, are they alien or threatened species, etc. The output will provide DOI (Digital Object Identifier) based stable identifiers for the communicating species found in eDNA. DOIs are connected to the taxonomic backbone of PlutoF and GBIF (https://www.gbif.org). In this way every DOI is accompanied by a taxon name which is still widely used for the communication of species. In the case of undescribed species, DOIs will soon be issued by the PlutoF system. SH matching service covers all Eukaryota by using rDNA ITS marker sequences and accompanied sample metadata.

SH matching output

Output archive includes the following files:

1. matches/ directory - folder including the main output files describing user's dataset matches to existing SHs, new compound clusters and SHs and HTML output
	a. matches_out_*.csv - data matrix describing user's sequences belonging to existing SHs or forming new SHs. These only include sequences which were placed into existing compound clusters;
	b. matches_1_out_*.csv - data matrix describing user's sequences splitting into new SHs. These include sequences placed outside of the current SH system and existing compound clusters. These sequences should be critically reviewed because they often include large proportion of chimeric and/or low quality sequences;
	c. matches_out_*.html - HTML output summarizing the data in matches CSV files.
2. err_[run_id].log - info and error messages in different analysis steps, mostly for debugging and providing support;
3. excluded_[run_id].txt - list of sequences excluded from the analysis with the reason why they were included;
4. source_[run_id]_names - translation table between user's sequence identifiers and system internal identifiers used throughout the analysis;
5. source_[run_id]_fastanames - list of duplicate sequences in user's input data;
6. krona_97.html - Krona HTML file krona_97.html can be interactively visualized with web browser. Piechart shows all those taxonomic units for which there is taxonomic information from matching SHs and compound clusters.
