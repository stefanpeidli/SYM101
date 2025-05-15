# How to download copy number data from DepMap

- Go to DepMap portal: https://depmap.org/portal/
- Click on "Downloads" in the top menu.
- Select "Custom Downloads".
- Choose the datasets you need
- Select cell lines of interest (HT-2 + the ones in this study, i.e. H1299, H1975, A549)
- Select options:
    - Exclude columns and rows of NA's from download files: yes
    - Add cell line metadata to download: yes
    - Merge into a single file: yes
- Click "Download" to get the data in CSV format.

To get expression data for all transcripts, go to "All Data" and download "OmicsExpressionAllGenesTPMLogp1Profile.csv".
A version subsetted to the relevant cell lines is available in the `resources/DepMap` directory.