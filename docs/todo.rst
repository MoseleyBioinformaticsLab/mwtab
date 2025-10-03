

Add some information into the documentation about how some of the validations are pretty strict to account for the major of cases, 
but if warranted could be ignored. For example, COLUMN_PRESSURE in the CHROMATOGRAPHY section is looking for a single number or 
range of numbers followed by a unit, but there might be some situations where the method is complex and thus the column pressure 
is not static, so something like "60 bar at starting conditions. 180 bar at %A" would be the value.

Add something to the documentation explaining how NMR_BINNED datasets will have both 'Metabolite' and 'Bin range(ppm)' keys in their data.

Make sure there is something in the documentation about how the results file is in nice parts for the mwtab object, but is not like that in the JSON.

Once the documentation is built, make sure the different style doc strings are being picked up correctly. Christian used a different style from me.

Note in the documentation that is not recommended to download JSON directly. Better to download mwTab and convert.

Note in documentation about compatability mode and that it is expected to be in compatability mode for validation. Might get unexpected results if not.

Add documentation about duplicate_keys parameter to MWTabFile class.

Add documentation about MWTabfile.study_id, analysis_id, and header.

Think about extending METABOLITES and EXTENDED blocks with an "Attributes" line like "Factors" in DATA block as a way to add more information about the columns themselves.
Hunter also wanted to consider adding things like the _factors properties into the JSON as well. For example, the _factors could be added into ['MS_METABOLITE_DATA'] under a 'Factors' key.

