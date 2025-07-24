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
name. The :ref:`api:Metadata Column Matching` page for metadata_column_matching.py 
contains more detailed information and examples on each class.


Regular Expressions
-------------------
In order to create some of the complex regular expressions needed to validate column 
names and values, a modular approach to constructing them was taken. This means the 
smaller building block regular expressions could have use outside of this particular 
one. There is also a function used heavily in creating these regular expressions, 
make_list_regex. The :ref:`api:Metadata Column Matching` page for metadata_column_matching.py 
contains more detailed information about using this function. There are too many 
regular expressions to explain each one, so the source code for creating them is 
included below. The names are largely self explanatory, and being able to see the 
regular expression for each name can help clear up confusion.

.. literalinclude:: ../src/mwtab/metadata_column_matching.py
   :start-after: # Comment for Sphinx to pull in regular expressions.
   :end-before: # Comment for Sphinx to find the end of regular expressions.


column_finders Dictionary
-------------------------
The column_finders dictionary is the culmination of the ColumnFinder class and 
regular expressions for the purpose of finding and evaluating/validating columns 
in Metabolomics Workbench datasets. The keys are standardized column names, and the 
values are a ColumnFinder class to find and validate that column. The entire library 
of datasets in the Metbolomics Workbench was used to determine the most relevant and 
most abundant columns that should go into the column_finders dictionary. The column_finders 
dictionary is likely to be useful for other similar metabolomics data, but may vary 
outside of that usecase. Some columns, such as database ID columns, like PubChem or 
KEGG, are likely to be more widely useable, but any user would have to test and 
make that determination for themselves. The list of standardized column names available 
in the dictionary is given below.

.. literalinclude:: standard_column_names.json




.. _Issues: https://github.com/MoseleyBioinformaticsLab/mwtab/issues