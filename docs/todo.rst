TODO
====

Add options to validate CLI and validate method in mwtab to save out the new JSON file.

Add --limit or --ignore option to validate to filter out certain types of errors/warnings. Need to first create some classifications to tag them with.

Think about extending METABOLITES and EXTENDED blocks with an "Attributes" line like "Factors" in DATA block as a way to add more information about the columns themselves.
Hunter also wanted to consider adding things like the _factors properties into the JSON as well. For example, the _factors could be added into ['MS_METABOLITE_DATA'] under a 'Factors' key.

Think about adding an "UNASSIGNED" data block for the datasets we found that have a results_file instead of having the data in the mwTab file.
Pretty sure most of these if not all are all unnassigned data where there are basically bins and no metabolite assinments.
