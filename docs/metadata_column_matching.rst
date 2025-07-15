Metadata Column Matching
========================
The mwTab format gives a lot of freedom to uploaders in what kind of metadata associated 
with the assigned metabolites can be given. While this has its merits, a big drawback 
is the amount of different names given for the same column of information. For example, 
using both 'm/z' and 'moverz'. Another drawback is the variety of value formats for 
some columns. For example, a column containing HMDB IDs might leave the 'HMDB' prefix 
off of the beginning of the ID, or the 'HMDB' prefix might be lower case. To address 
these issues a significant amount of work was done to data mine the most used and 
most useful columns and using this information code was created to be able to identify the columns under 
their various names as well as evaluate the format of their values. 

This code is utilized in the validation to inform users about the state of their 
metadata columns in the METABOLITES section. In order to do this we chose standardized 
column names and value formats based on the information obtained during data mining, 
so there were some executive decisions we made in deciding what the standard names and 
values are or should be. We did our best to base them on the data already in the 
Metabolomics Workbench.

We did our best to test this code on the data in the Metabolomics Workbench, but it 
is not perfect. Due to the nature of regular expressions and trying to match a variety 
of names and values, inevitably there are going to be false positives and negatives. 
You may see a validation warning about a column name or value in error. If you do 
find incorrect warnings or issues that are not being warned about in the validation, 
please visit our GitHub Issues_ and create a new issue for it.




Reusing Match Code
~~~~~~~~~~~~~~~~~~
Although the mwtab package is largely meant to be used through its CLI, the code 
created for column name and value matching could have many more use cases outside 
of this package. This section will explain the different parts of the code and how 
they might be used outside of this package.

The 3 significant parts are:
1. ColumnFinder class
2. Regular expressions and associated functions
3. Dictionary of ColumnFinders for specific mwTab columns


ColumnFinder Class
------------------
Column matching has 2 components, matching column names and matching column values, 
therefore, a ColumnFinder has a NameMatcher and ValueMatcher to handle those functions. 
All NameMatcher attributes are lists of strings or lists of lists of strings and all 
are used in its only method, dict_match. All ValueMatcher attributes are strings and 
all are used in its only method, series_match. ColumnFinder takes both a NameMatcher 
and ValueMatcher, along with a standard_name attribute. The "standard_name" attribute 
is not used by any method and is simply there to tie the instance of the class to a 
name. The :doc:`api` page for metadata_column_matching.py contains more detailed information 
and examples on each class.


Regular Expressions
-------------------







.. _Issues: https://github.com/MoseleyBioinformaticsLab/mwtab/issues