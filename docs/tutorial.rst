The mwtab Tutorial
==================

The :mod:`mwtab` package provides classes and other facilities for parsing,
accessing, and manipulating data stored in ``mwTab`` and ``JSON`` representation
of ``mwTab`` formats.

Also, the :mod:`mwtab` package provides simple command-line interface to convert
between ``mwTab`` and its ``JSON`` representation as well as validate consistency
of the files.


Brief mwTab Format Overview
~~~~~~~~~~~~~~~~~~~~~~~~~~~


   .. note::

      For full official specification see the following link (``mwTab file specification``):
      http://www.metabolomicsworkbench.org/data/tutorials.php


The ``mwTab`` formatted files consist of multiple blocks. Each new block starts with ``#``.

   * Some of the blocks contain only "key-value"-like pairs.

   .. code-block:: none

      #METABOLOMICS WORKBENCH STUDY_ID:ST000001 ANALYSIS_ID:AN000001
      VERSION             	1
      CREATED_ON          	2016-09-17
      #PROJECT
      PR:PROJECT_TITLE                 	FatB Gene Project
      PR:PROJECT_TYPE                  	Genotype treatment
      PR:PROJECT_SUMMARY               	Experiment to test the consequence of a mutation at the FatB gene (At1g08510)
      PR:PROJECT_SUMMARY               	the wound-response of Arabidopsis

   .. note::

      ``*_SUMMARY`` "key-value"-like pairs are typically span through multiple lines.


   * ``#SUBJECT_SAMPLE_FACTORS`` block is specially formatted, i.e. it contains header
     specification and tab-separated values.

   .. code-block:: none

      #SUBJECT_SAMPLE_FACTORS:         	SUBJECT(optional)[tab]SAMPLE[tab]FACTORS(NAME:VALUE pairs separated by |)[tab]Additional sample data
      SUBJECT_SAMPLE_FACTORS           	-	LabF_115873	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded
      SUBJECT_SAMPLE_FACTORS           	-	LabF_115878	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded
      SUBJECT_SAMPLE_FACTORS           	-	LabF_115883	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded
      SUBJECT_SAMPLE_FACTORS           	-	LabF_115888	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded
      SUBJECT_SAMPLE_FACTORS           	-	LabF_115893	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded
      SUBJECT_SAMPLE_FACTORS           	-	LabF_115898	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded


   * ``#MS_METABOLITE_DATA`` (results) block contains ``Samples`` identifiers, ``Factors`` identifiers
     as well as tab-separated data between ``*_START`` and ``*_END``.

   .. code-block:: none

      #MS_METABOLITE_DATA
      MS_METABOLITE_DATA:UNITS	Peak height
      MS_METABOLITE_DATA_START
      Samples	LabF_115904	LabF_115909	LabF_115914	LabF_115919	LabF_115924	LabF_115929	LabF_115842	LabF_115847	LabF_115852	LabF_115857	LabF_115862	LabF_115867	LabF_115873	LabF_115878	LabF_115883	LabF_115888	LabF_115893	LabF_115898	LabF_115811	LabF_115816	LabF_115821	LabF_115826	LabF_115831	LabF_115836
      Factors	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded	Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded
      1_2_4-benzenetriol	1874.0000	3566.0000	1945.0000	1456.0000	2004.0000	1995.0000	4040.0000	2432.0000	2189.0000	1931.0000	1307.0000	2880.0000	2218.0000	1754.0000	1369.0000	1201.0000	3324.0000	1355.0000	2257.0000	1718.0000	1740.0000	3472.0000	2054.0000	1367.0000
      1-monostearin	987.0000	450.0000	1910.0000	549.0000	1032.0000	902.0000	393.0000	705.0000	100.0000	481.0000	265.0000	120.0000	1185.0000	867.0000	676.0000	569.0000	579.0000	387.0000	1035.0000	789.0000	875.0000	224.0000	641.0000	693.0000
      ...
      MS_METABOLITE_DATA_END

   * ``#METABOLITES`` metadata block contains a header specifying fields and
     tab-separated data between ``*_START`` and ``*_END``.

   .. code-block:: none

      #METABOLITES
      METABOLITES_START
      metabolite_name	moverz_quant	ri	ri_type	pubchem_id	inchi_key	kegg_id	other_id	other_id_type
      1,2,4-benzenetriol	239	522741	Fiehn	10787		C02814	205673	BinBase
      1-monostearin	399	959625	Fiehn	107036		D01947	202835	BinBase
      2-hydroxyvaleric acid	131	310750	Fiehn	98009			218773	BinBase
      3-phosphoglycerate	299	611619	Fiehn	724		C00597	217821	BinBase
      ...
      METABOLITES_END

   * ``#NMR_BINNED_DATA`` metadata block contains a header specifying fields and
     tab-separated data between ``*_START`` and ``*_END``.

   .. code-block:: none

      #NMR_BINNED_DATA
      NMR_BINNED_DATA_START
      Bin range(ppm)	CDC029	CDC030	CDC032	CPL101	CPL102	CPL103	CPL201	CPL202	CPL203	CDS039	CDS052	CDS054
      0.50...0.56	0.00058149	1.6592	0.039301	0	0	0	0.034018	0.0028746	0.0021478	0.013387	0	0
      0.56...0.58	0	0.74267	0	0.007206	0	0	0	0	0	0	0	0.0069721
      0.58...0.60	0.051165	0.8258	0.089149	0.060972	0.026307	0.045697	0.069541	0	0	0.14516	0.057489	0.042255
      ...
      NMR_BINNED_DATA_END

   * Order of metadata and data blocks (MS)

   .. code-block:: none

      #METABOLOMICS WORKBENCH
      VERSION             	1
      CREATED_ON          	2016-09-17
      ...
      #PROJECT
      ...
      #STUDY
      ...
      #SUBJECT
      ...
      #SUBJECT_SAMPLE_FACTORS:         	SUBJECT(optional)[tab]SAMPLE[tab]FACTORS(NAME:VALUE pairs separated by |)[tab]Additional sample data
      ...
      #COLLECTION
      ...
      #TREATMENT
      ...
      #SAMPLEPREP
      ...
      #CHROMATOGRAPHY
      ...
      #ANALYSIS
      ...
      #MS
      ...
      #MS_METABOLITE_DATA
      MS_METABOLITE_DATA:UNITS	peak area
      MS_METABOLITE_DATA_START
      ...
      MS_METABOLITE_DATA_END
      #METABOLITES
      METABOLITES_START
      ...
      METABOLITES_END
      #END


Using mwtab as a Library
~~~~~~~~~~~~~~~~~~~~~~~~


Importing mwtab Package
-----------------------

If the :mod:`mwtab` package is installed on the system, it can be imported:

>>> import mwtab


Constructing MWTabFile Generator
--------------------------------

The :mod:`~mwtab.fileio` module provides the :func:`~mwtab.fileio.read_files`
generator function that yields :class:`~mwtab.mwtab.MWTabFile` instances. Constructing a
:class:`~mwtab.mwtab.MWTabFile` generator is easy - specify the path to a local ``mwTab`` file,
directory of files, archive of files:

>>> import mwtab
>>>
>>> mwtfile_gen = mwtab.read_files("ST000001_AN000001.txt")  # single mwTab file
>>>
>>> mwtfiles_gen = mwtab.read_files("ST000001_AN000001.txt", "ST000002_AN000002.txt") # several mwTab files
>>>
>>> mwtdir_gen = mwtab.read_files("mwtabfiles_dir")   # directory of mwTab files
>>>
>>> mwtzip_gen = mwtab.read_files("mwtabfiles.zip")  # archive of mwTab files
>>>
>>> mwturl_gen = mwtab.read_files("1", "2")            # ANALYSIS_ID of mwTab file
>>>


Processing MWTabFile Generator
------------------------------

The :class:`~mwtab.mwtab.MWTabFile` generator can be processed in several ways:

   * Feed it to a for-loop and process one file at a time:

   >>> for mwtfile in mwtdir_gen:
   >>>     print("STUDY_ID:", mwtfile.study_id)       # print STUDY_ID
   >>>     print("ANALYSIS_ID", mwtfile.analysis_id)  # print ANALYSIS_ID
   >>>     print("SOURCE", mwtfile.source)            # print source
   >>>     for block_name in mwtfile:                 # print names of blocks
   >>>          print("\t", block_name)

   .. note:: Once the generator is consumed, it becomes empty and needs to be created again.

   * Since the :class:`~mwtab.mwtab.MWTabFile` generator behaves like an iterator,
     we can call the :py:func:`next` built-in function:

   >>> mwtfile1 = next(mwtdir_gen)
   >>> mwtfile2 = next(mwtdir_gen)
   >>> ...

   .. note:: Once the generator is consumed, :py:class:`StopIteration` will be raised.

   * Convert the :class:`~mwtab.mwtab.MWTabFile` generator into a :py:class:`list` of
     :class:`~mwtab.mwtab.MWTabFile` objects:

   >>> mwtfiles_list = list(mwtdir_gen)
   >>>

   .. warning:: Do not convert the :class:`~mwtab.mwtab.MWTabFile` generator into a
                :py:class:`list` if the generator can yield a large number of files, e.g.
                several thousand, otherwise it can consume all available memory.


Accessing Data From a Single MWTabFile
--------------------------------------

Since a :class:`~mwtab.mwtab.MWTabFile` is a Python :py:class:`collections.OrderedDict`,
data can be accessed and manipulated as with any regular Python :py:class:`dict` object
using bracket accessors.


   * Accessing top-level "keys" in :class:`~mwtab.mwtab.MWTabFile`:

   >>> list(mwtfile.keys())
   [
       'METABOLOMICS WORKBENCH',
       'PROJECT',
       'STUDY',
       'SUBJECT',
       'SUBJECT_SAMPLE_FACTORS',
       'COLLECTION',
       'TREATMENT',
       'SAMPLEPREP',
       'CHROMATOGRAPHY',
       'ANALYSIS', 'MS',
       'MS_METABOLITE_DATA',
       'METABOLITES'
   ]


   * Accessing individual blocks in :class:`~mwtab.mwtab.MWTabFile`:

   >>> mwtfile["PROJECT"]
   OrderedDict([
       ('PROJECT_TITLE', 'Intestinal Samples II pre/post transplantation'),
       ('PROJECT_TYPE', 'Human intestinal samples'),
       ('PROJECT_SUMMARY', 'Intestinal Samples II pre/post transplantation'),
       ('INSTITUTE', 'University of California, Davis'),
       ('DEPARTMENT', 'Davis Genome Center'),
       ('LABORATORY', 'Fiehn'),
       ('LAST_NAME', 'Fiehn'),
       ('FIRST_NAME', 'Oliver'),
       ('ADDRESS', '451 E. Health Sci. Drive, Davis, California 95616, USA'),
       ('EMAIL', 'ofiehn@ucdavis.edu'),
       ('PHONE', '-')
   ])


   * Accessing individual "key-value" pairs within blocks:

   >>> mwtfile["PROJECT"]["INSTITUTE"]
   'University of California, Davis'


   * Accessing data in ``#SUBJECT_SAMPLE_FACTORS`` block:

   >>> mwtfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]
   [
       OrderedDict([
           ('subject_type', '-'),
           ('local_sample_id', 'LabF_115873'),
           ('factors', 'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded'),
           ('additional_sample_data', '')]),
       OrderedDict([('subject_type', '-'),
           ('local_sample_id', 'LabF_115878'),
           ('factors', 'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded'),
           ('additional_sample_data', '')]),
       ...
   ]


   >>> mwtfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"][0]["local_sample_id"]
   'LabF_115873'


   * Accessing data in ``#MS_METABOLITE_DATA`` block:

   >>> mwtfile["MS_METABOLITE_DATA"]["METABOLITE_DATA:UNITS"]
   'peak area'


   >>> mwtfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Samples"]
   [
       'LabF_115904',
       'LabF_115909',
       'LabF_115914',
       'LabF_115919',
       'LabF_115924',
       'LabF_115929',
       'LabF_115842',
       'LabF_115847',
       'LabF_115852',
       'LabF_115857',
       'LabF_115862',
       'LabF_115867',
       'LabF_115873',
       'LabF_115878',
       'LabF_115883',
       'LabF_115888',
       'LabF_115893',
       'LabF_115898',
       'LabF_115811',
       'LabF_115816',
       'LabF_115821',
       'LabF_115826',
       'LabF_115831',
       'LabF_115836'
   ]


   >>> mwtfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Factors"]
   [
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:fatb-ko KD; At1g08510 | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded',
       'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Wounded'
   ]


   >>> mwtfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["DATA"]
   [
       OrderedDict([
           ('metabolite_name', '1_2_4-benzenetriol'),
           ('LabF_115904', '1874.0000'),
           ('LabF_115909', '3566.0000'),
           ('LabF_115914', '1945.0000'),
           ('LabF_115919', '1456.0000'),
           ('LabF_115924', '2004.0000'),
           ('LabF_115929', '1995.0000'),
           ('LabF_115842', '4040.0000'),
           ('LabF_115847', '2432.0000'),
           ('LabF_115852', '2189.0000'),
           ('LabF_115857', '1931.0000'),
           ('LabF_115862', '1307.0000'),
           ('LabF_115867', '2880.0000'),
           ('LabF_115873', '2218.0000'),
           ('LabF_115878', '1754.0000'),
           ('LabF_115883', '1369.0000'),
           ('LabF_115888', '1201.0000'),
           ('LabF_115893', '3324.0000'),
           ('LabF_115898', '1355.0000'),
           ('LabF_115811', '2257.0000'),
           ('LabF_115816', '1718.0000'),
           ('LabF_115821', '1740.0000'),
           ('LabF_115826', '3472.0000'),
           ('LabF_115831', '2054.0000'),
           ('LabF_115836', '1367.0000')]),
       OrderedDict([
           ('metabolite_name', '1-monostearin'),
           ('LabF_115904', '987.0000'),
           ('LabF_115909', '450.0000'),
           ('LabF_115914', '1910.0000'),
           ('LabF_115919', '549.0000'),
           ('LabF_115924', '1032.0000'),
           ('LabF_115929', '902.0000'),
           ('LabF_115842', '393.0000'),
           ('LabF_115847', '705.0000'),
           ('LabF_115852', '100.0000'),
           ('LabF_115857', '481.0000'),
           ('LabF_115862', '265.0000'),
           ('LabF_115867', '120.0000'),
           ('LabF_115873', '1185.0000'),
           ('LabF_115878', '867.0000'),
           ('LabF_115883', '676.0000'),
           ('LabF_115888', '569.0000'),
           ('LabF_115893', '579.0000'),
           ('LabF_115898', '387.0000'),
           ('LabF_115811', '1035.0000'),
           ('LabF_115816', '789.0000'),
           ('LabF_115821', '875.0000'),
           ('LabF_115826', '224.0000'),
           ('LabF_115831', '641.0000'),
           ('LabF_115836', '693.0000')]),
       ...
   ]


Manipulating Data From a Single MWTabFile
-----------------------------------------

In order to change values within :class:`~mwtab.mwtab.MWTabFile`, descend into
the appropriate level using square bracket accessors and set a new value.

   * Change regular "key-value" pairs:

   >>> mwtfile["PROJECT"]["PHONE"]
   '-'
   >>> mwtfile["PROJECT"]["PHONE"] = "1-530-754-8258"
   >>> mwtfile["PROJECT"]["PHONE"]
   '1-530-754-8258'

   * Change ``#SUBJECT_SAMPLE_FACTORS`` values:

   >>> mwtfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"][0]
   OrderedDict([
       ('subject_type', '-'),
       ('local_sample_id', 'LabF_115873'),
       ('factors', 'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded'),
       ('additional_sample_data', '')
   ])
   >>> mwtfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"][0]["additional_sample_data"] = "Additional details"
   >>> mwtfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"][0]
   OrderedDict([
       ('subject_type', '-'),
       ('local_sample_id', 'LabF_115873'),
       ('factors', 'Arabidopsis Genotype:Wassilewskija (Ws) | Plant Wounding Treatment:Control - Non-Wounded'),
       ('additional_sample_data', 'Additional details')
   ])

Printing a MWTabFile and its Components
---------------------------------------

   * Print entire file in ``mwTab`` format.

   >>> mwtfile.print_file(file_format="mwtab")
   #METABOLOMICS WORKBENCH STUDY_ID:ST000001 ANALYSIS_ID:AN000001
   VERSION             	1
   CREATED_ON          	2016-09-17
   #PROJECT
   PR:PROJECT_TITLE                    	FatB Gene Project
   PR:PROJECT_TYPE                     	Genotype treatment
   PR:PROJECT_SUMMARY                  	Experiment to test the consequence of a mutation at the FatB gene (At1g08510)
   PR:PROJECT_SUMMARY                  	the wound-response of Arabidopsis
   PR:INSTITUTE                        	University of California, Davis
   PR:DEPARTMENT                       	Davis Genome Center
    ...

  * Print entire file in ``JSON`` format.

   >>> mwtfile.print_file(file_format="json")
   {
       "METABOLOMICS WORKBENCH": {
           "HEADER": "#METABOLOMICS WORKBENCH STUDY_ID:ST000001 ANALYSIS_ID:AN000001",
           "STUDY_ID": "ST000001",
           "ANALYSIS_ID": "AN000001",
           "VERSION": "1",
           "CREATED_ON": "2016-09-17"
       },
       "PROJECT": {
           "PROJECT_TITLE": "FatB Gene Project",
           "PROJECT_TYPE": "Genotype treatment",
           "PROJECT_SUMMARY": "Experiment to test the consequence of a mutation at the FatB gene (At1g08510)\nthe wound-response of Arabidopsis",
           "INSTITUTE": "University of California, Davis",
           "DEPARTMENT": "Davis Genome Center", ...
       },
       ...
   }

   * Print single block in ``mwTab`` format.

   >>> mwtfile.print_block("STUDY", file_format="mwtab")
   ST:STUDY_TITLE                      	Fatb Induction Experiment (FatBIE)
   ST:STUDY_TYPE                       	Genotype treatment
   ST:STUDY_SUMMARY                    	This experiment tests the consequence of a mutation at the FatB gene
   ST:STUDY_SUMMARY                    	in the wound-response of Arabidopsis. The FatB mutant allele (fatb KD J.
   ST:STUDY_SUMMARY                    	(Plant Cell 2003, Vol 15, 1020-1033)) was obtained from Dr. Katayonn Dehesh,
   ST:STUDY_SUMMARY                    	of California, Davis, Davis, CA. This allele is in the Ws background.The
   ST:STUDY_SUMMARY                    	growth conditions are as follows: 1. Seeds (between 14 and 16) are sown on
   ST:STUDY_SUMMARY                    	in 100 x 100 x 15mm square Falcon Petri Dishes (Fisher Scientific, catalogue
   ST:STUDY_SUMMARY                    	Seeds were arranged on the plates in a single horizontal line at the 1-cm mark
   ST:STUDY_SUMMARY                    	the top of the plate.2. Each plate contains between 20 and 25-ml of sterile MS
   ST:STUDY_SUMMARY                    	containing 0.1% (w/v) sucrose.3. Prior to sowing, seeds were sterilized by
   ST:STUDY_SUMMARY                    	for 1 minute at room temperature with a 300-l solution of 50% (v/v) ethanol,
   ST:STUDY_SUMMARY                    	solution was removed and replaced with a 300-l solution consisting of 1% (v/v)
   ST:STUDY_SUMMARY                    	20 (Fischer BioReagents, catalogue #BP33750), and 50% (v/v) bleach solution
   ST:STUDY_SUMMARY                    	and incubated at room temperature for 10-minutes. The seeds were then washed
   ST:STUDY_SUMMARY                    	three changes of 0.3-ml of sterile water.
   ST:INSTITUTE                        	University of California, Davis
   ST:DEPARTMENT                       	Davis Genome Center
   ST:LABORATORY                       	Fiehn
   ST:LAST_NAME                        	Kind
   ST:FIRST_NAME                       	Tobias
   ST:ADDRESS                          	451 E. Health Sci. Drive, Davis, CA 95616, USA
   ST:EMAIL                            	tkind@ucdavis.edu
   ST:PHONE                            	-
   ST:SUBMIT_DATE                      	2013-01-15
   ST:NUM_GROUPS                       	4
   ST:TOTAL_SUBJECTS                   	24

   * Print single block in ``JSON`` format.

   >>> mwtfile.print_block("STUDY", file_format="json")
   {
       "STUDY_TITLE": "Fatb Induction Experiment (FatBIE)",
       "STUDY_TYPE": "Genotype treatment",
       "STUDY_SUMMARY": "This experiment tests the consequence of a mutation at the FatB gene\nin the wound-response of Arabidopsis. The FatB mutant allele (fatb KD J.\n(Plant Cell 2003, Vol 15, 1020-1033)) was obtained from Dr. Katayonn Dehesh,\nof California, Davis, Davis, CA. This allele is in the Ws background.The\ngrowth conditions are as follows: 1. Seeds (between 14 and 16) are sown on\nin 100 x 100 x 15mm square Falcon Petri Dishes (Fisher Scientific, catalogue\nSeeds were arranged on the plates in a single horizontal line at the 1-cm mark\nthe top of the plate.2. Each plate contains between 20 and 25-ml of sterile MS\ncontaining 0.1% (w/v) sucrose.3. Prior to sowing, seeds were sterilized by\nfor 1 minute at room temperature with a 300-l solution of 50% (v/v) ethanol,\nsolution was removed and replaced with a 300-l solution consisting of 1% (v/v)\n20 (Fischer BioReagents, catalogue #BP33750), and 50% (v/v) bleach solution\nand incubated at room temperature for 10-minutes. The seeds were then washed\nthree changes of 0.3-ml of sterile water.",
       "INSTITUTE": "University of California, Davis",
       "DEPARTMENT": "Davis Genome Center",
       "LABORATORY": "Fiehn",
       "LAST_NAME": "Kind",
       "FIRST_NAME": "Tobias",
       "ADDRESS": "451 E. Health Sci. Drive, Davis, CA 95616, USA",
       "EMAIL": "tkind@ucdavis.edu",
       "PHONE": "-",
       "SUBMIT_DATE": "2013-01-15",
       "NUM_GROUPS": "4",
       "TOTAL_SUBJECTS": "24"
   }


Writing data from a MWTabFile object into a file
------------------------------------------------
Data from a :class:`~mwtab.mwtab.MWTabFile` can be written into file
in original ``mwTab`` format or in equivalent JSON format using
:meth:`~mwtab.mwtab.MWTabFile.write()`:

   * Writing into a ``mwTab`` formatted file:

   >>> with open("ST000001_AN000001_modified.txt", "w") as outfile:
   ...     mwtfile.write(outfile, file_format="mwtab")
   >>>

   * Writing into a ``JSON`` file:

   >>> with open("ST000001_AN000001_modified.json", "w") as outfile:
   ...     mwtfile.write(outfile, file_format="json")
   >>>


Converting mwTab Files
----------------------

``mwTab`` files can be converted between the ``mwTab`` file format and their ``JSON``
representation using the :mod:`mwtab.converter` module.

One-to-one file conversions
***************************

   * Converting from the ``mwTab`` file format into its equivalent ``JSON`` file format:

   .. code-block:: python

      from mwtab.converter import Converter

      # Using valid ANALYSIS_ID to access file from URL: from_path="1"
      converter = Converter(from_path="1", to_path="AN000001.json",
                            from_format="mwtab", to_format="json"))
      converter.convert()


   * Converting from JSON file format back to ``mwTab`` file format:

   .. code-block:: python

      from mwtab.converter import Converter

      converter = Converter(from_path="AN000001.json", to_path="AN000001.txt",
                            from_format="json", to_format="mwtab"))
      converter.convert()


Many-to-many files conversions
******************************

   * Converting from the directory of ``mwTab`` formatted files into their equivalent
     ``JSON`` formatted files:

   .. code-block:: python

      from mwtab.converter import Converter

      converter = Converter(from_path="mwtdir_mwtab",
                            to_path="mwtdir_json",
                            from_format="mwtab",
                            to_format="json")
      converter.convert()

   * Converting from the directory of ``JSON`` formatted files into their equivalent
     ``mwTab`` formatted files:

   .. code-block:: python

      from mwtab.converter import Converter

      converter = Converter(from_path="mwtdir_json",
                            to_path="mwtdir_mwtab",
                            from_format="json",
                            to_format="mwtab"))
      converter.convert()


.. note:: Many-to-many files and one-to-one file conversions are available.
          See :mod:`mwtab.converter` for full list of available conversions.


Command-Line Interface
~~~~~~~~~~~~~~~~~~~~~~

The mwtab Command-Line Interface provides the following functionality:
   * Convert from the ``mwTab`` file format into its equivalent ``JSON`` file format and vice versa.
   * Validate the ``mwTab`` formatted file.

.. code-block:: none

   The mwtab command-line interface
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Usage:
       mwtab -h | --help
       mwtab --version
       mwtab convert (<from-path> <to-path>) [--from-format=<format>] [--to-format=<format>] [--validate] [--verbose]
       mwtab validate <from-path> [--verbose]

   Options:
       -h, --help                      Show this screen.
       --version                       Show version.
       --verbose                       Print what files are processing.
       --validate                      Validate the mwTab file.
       --from-format=<format>          Input file format, available formats: mwtab, json [default: mwtab].
       --to-format=<format>            Output file format, available formats: mwtab, json [default: json].
       --mw-rest=<url>                 URL to MW REST interface [default: "http://www.metabolomicsworkbench.org/rest/study/analysis_id/{}/mwtab/txt"].

Converting ``mwTab`` files in bulk
----------------------------------

One-to-one file conversions
***************************

   * Convert from a local file in ``mwTab`` format to a local file in ``JSON`` format:

   .. code:: bash

      $ python3 -m mwtab convert ST000001_AN000001.txt ST000001_AN000001.json \
                --from-format=mwtab --to-format=json

   * Convert from a local file in ``JSON`` format to a local file in ``mwTab`` format:

   .. code:: bash

      $ python3 -m mwtab convert ST000001_AN000001.json ST000001_AN000001.txt \
                --from-format=json --to-format=mwtab

   * Convert from a compressed local file in ``mwTab`` format to a compressed local file in ``JSON`` format:

   .. code:: bash

      $ python3 -m mwtab convert ST000001_AN000001.txt.gz ST000001_AN000001.json.gz \
                --from-format=mwtab --to-format=json

   * Convert from a compressed local file in ``JSON`` format to a compressed local file in ``mwTab`` format:

   .. code:: bash

      $ python3 -m mwtab convert ST000001_AN000001.json.gz ST000001_AN000001.txt.gz \
                --from-format=json --to-format=mwtab

   * Convert from a uncompressed URL file in ``mwTab`` format to a compressed local file in ``JSON`` format:

   .. code:: bash

      $ python3 -m mwtab convert 1 ST000001_AN000001.json.bz2 \
                --from-format=mwtab --to-format=json

   .. note:: See :mod:`mwtab.converter` for full list of available conversions.


Many-to-many files conversions
******************************

   * Convert from a directory of files in ``mwTab`` format to a directory of files in ``JSON`` format:

   .. code-block:: none

      $ python3 -m mwtab convert mwtabfiles_mwtab mwtabfiles_json \
                --from-format=mwtab --to-format=json

   * Convert from a directory of files in ``JSON`` format to a directory of files in ``mwTab`` format:

   .. code-block:: none

      $ python3 -m mwtab convert mwtabfiles_json mwtabfiles_mwtab \
                --from-format=json --to-format=mwtab

   * Convert from a directory of files in ``mwTab`` format to a zip archive of files in ``JSON`` format:

   .. code-block:: none

      $ python3 -m mwtab convert mwtabfiles_mwtab mwtabfiles_json.zip \
                --from-format=mwtab --to-format=json

   * Convert from a compressed tar archive of files in ``JSON`` format to a directory of files in ``mwTab`` format:

   .. code-block:: none

      $ python3 -m mwtab convert mwtabfiles_json.tar.gz mwtabfiles_mwtab \
                --from-format=json --to-format=mwtab

   * Convert from a zip archive of files in ``mwTab`` format to a compressed tar archive of files in ``JSON`` format:

   .. code-block:: none

      $ python3 -m mwtab convert mwtabfiles_mwtab.zip mwtabfiles_json.tar.bz2 \
                --from-format=mwtab --to-format=json

   .. note:: See :mod:`mwtab.converter` for full list of available conversions.


Validating ``mwTab`` files
--------------------------

The :mod:`mwtab` package provides the :func:`~mwtab.validator.validate_file` function
that can validate files based on a ``JSON`` schema definition. The :mod:`mwtab.mwschema`
contains schema definitions for every block (section) of ``mwTab`` formatted file, i.e.
it lists the types of attributes (e.g. :py:class:`str` as well as specifies which keys are
optional and which are required).

   * To validate file(s) simply call the ``validate`` command and provide path to file(s):

   .. code-block:: none

      $ python3 -m mwtab validate path
