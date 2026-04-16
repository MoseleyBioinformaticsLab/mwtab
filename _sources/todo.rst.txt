TODO
====

Add --limit or --ignore option to validate to filter out certain types of errors/warnings. Need to first create some classifications to tag them with.

Think about extending METABOLITES and EXTENDED blocks with an "Attributes" line like "Factors" in DATA block as a way to add more information about the columns themselves.
Hunter also wanted to consider adding things like the _factors properties into the JSON as well. For example, the _factors could be added into ['MS_METABOLITE_DATA'] under a 'Factors' key.

Think about adding an "UNASSIGNED" data block for the datasets we found that have a results_file instead of having the data in the mwTab file.
Pretty sure most of these if not all are all unnassigned data where there are basically bins and no metabolite assinments.


Had some new undocumented requirements pop up when submitting the Helsley data. The following is an email from Eoin about it:
Hi Travis,
You should have received an automated email regarding the completion of your submission.
There were a couple of problems with the mwtab files. Our submission system requires a single common study design block for all analyses that have the same set of headings for all samples (even if some of them are blank) Headings like Sex and Resected tissue type have to include all samples. 
Also, the MW needs a mandatory column which indicates the source of each sample-I added one and designated Liver as the source of each patient sample. The “Show all samples link” in the main page (https://dev.metabolomicsworkbench.org:22222/data/subject_fetch.php?STUDY_ID=ST004733&STUDY_TYPE=MS&RESULT_TYPE=1&Access=IloQ2417) shows the proper layout.
The other issue was that there were a large number of unassigned annotations (m/z_rt features) that may not be added to an mwtab file. These are submitted as tab-delimited text files and saved with the raw data. I pulled these out of the mwtab files and saved them, one for each analysis:
ST004733_AN008002_Results.txt (124.1K)
ST004733_AN008003_Results.txt (249.1K)
ST004733_AN008004_Results.txt (135K)
ST004733_AN008005_Results.txt (254.2K)
ST004733_AN008006_Results.txt (109.5K)
The can be seen  under the “Download raw/supplementary data” on the main page.

Need to address these. New check in validation for all SSF to have the same factors. 
Must have a new required factor "Sample source" for every sample. It might not have to have that name, but that is 
the name most common in the latest deposited datasets.
Add check in validation that Metabolite names aren't just numbers or rt_mz or something like that. 
Those now have to be submitted as results files.

