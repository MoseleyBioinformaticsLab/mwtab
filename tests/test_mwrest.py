import pytest
from mwtab.mwrest import GenericMWURL, analysis_ids, study_ids, generate_mwtab_urls, MWRESTFile


def test_study_analysis():
    an_ids = analysis_ids()
    assert an_ids
    st_ids = study_ids()
    assert st_ids


@pytest.mark.parametrize("kwds", [
    ({"context": "study",
      "input_item": "analysis_id",
      "input_value": "AN000002",
      "output_item": "mwtab",
      'output_format': "txt"}),
    ({"context": "study",
      "input_item": "study_id",
      "input_value": "ST000001",
      "output_item": "mwtab",
      'output_format': "txt"}),
    ({"base_url": "https://www.test.org/rest/",
      "context": "study",
      "input_item": "study_id",
      "input_value": "ST000001",
      "output_item": "mwtab",
      'output_format': "txt"}),
])
def test_mwrest(kwds):
    test_mwurl = GenericMWURL(kwds)
    assert test_mwurl.url == test_mwurl.base_url + "/".join([
        kwds["context"],
        kwds["input_item"],
        kwds["input_value"],
        kwds["output_item"],
        kwds.get("output_format") or ""
    ])


@pytest.mark.parametrize("kwds", [
    ({"context": "study",
      "input_item": "analysis_id",
      "input_value": "ST000001",
      "output_item": "mwtab",
      'output_format': "txt"}),
    ({"context": "moverz",
      "input_item": "LIPIDS",
      "m/z_value": 49,
      "ion_type_value": "M+H",
      "m/z_tolerance_value": 0.1,
      'output_format': "txt"}),
    ({"context": "exactmass",
      "LIPID_abbreviation": "Test",
      "ion_type_value": "M+H"}),
])
def test_fail_mwrest(kwds):
    try:
        GenericMWURL(kwds)
        assert False
    except Exception as e:
        assert type(e) == ValueError


def test_generate_mwtab_urls():
    assert next(generate_mwtab_urls(['ST000001'])) == 'https://www.metabolomicsworkbench.org/rest/study/study_id/ST000001/mwtab/txt'

def test_generate_mwtab_urls_exception(mocker):
    mocker.patch('mwtab.mwrest.fileio._return_correct_yield', side_effect = [Exception])
    with pytest.raises(Exception):
        next(generate_mwtab_urls(['ST000001']))

def test_GenericMWURL_missing_context():
    with pytest.raises(KeyError, match = r'context'):
        GenericMWURL({"input_item": "analysis_id",
                      "input_value": "AN000001",
                      "output_item": "mwtab",
                      'output_format': "txt"})

def test_GenericMWURL_bad_context():
    with pytest.raises(ValueError, match = r'Error: Invalid/missing context'):
        GenericMWURL({"context": "asdf",
                      "input_item": "analysis_id",
                      "input_value": "AN000001",
                      "output_item": "mwtab",
                      'output_format': "txt"})

def test_GenericMWURL_missing_input_item():
    with pytest.raises(KeyError, match = r"Missing input item\(s\): \{'input_item'\}"):
        GenericMWURL({"context": "study",
                      "input_value": "AN000001",
                      "output_item": "mwtab",
                      'output_format': "txt"})

def test_GenericMWURL_bad_input_item():
    with pytest.raises(ValueError, match = r"Invalid input item"):
        GenericMWURL({"context": "study",
                      "input_item": "asdf",
                      "input_value": "AN000001",
                      "output_item": "mwtab",
                      'output_format': "txt"})

@pytest.mark.parametrize("kwds, message", [
    ({"context": "study",
      "input_item": "analysis_id",
      "input_value": "AN000001",
      "output_item": ['asdf'],
      'output_format': "txt"}, 
      r'Invalid output items. Study only takes a single output item.'),
    
    ({"context": "compound",
      "input_item": "formula",
      "input_value": "C3H6O2",
      "output_item": ['asdf'],
      'output_format': "txt"},
      r"Invalid output item\(s\): \{'asdf'\}"),
    
    ({"context": "study",
      "input_item": "analysis_id",
      "input_value": "AN000001",
      "output_item": "asdf",
      'output_format': "txt"},
      r'Invalid output item'),
])
def test_GenericMWURL_bad_output_item(kwds, message):
    with pytest.raises(ValueError, match = message):
        GenericMWURL(kwds)


@pytest.mark.parametrize("kwds, error, message", [
    ({"context": "moverz",
      "input_item": "LIPIDS",
      "ion_type_value": 'M+H',
      'm/z_tolerance_value': ".01"}, 
      KeyError,
      r"Missing input item\(s\): \{'m/z_value'\}"),
    
    ({"context": "moverz",
      "input_item": "asdf",
      "m/z_value": "58",
      "ion_type_value": 'M+H',
      'm/z_tolerance_value': ".01"},
      ValueError,
      r"Invalid input item"),
    
    ({"context": "moverz",
      "input_item": "LIPIDS",
      "m/z_value": "58",
      "ion_type_value": 'asdf',
      'm/z_tolerance_value': "txt"},
      ValueError,
      r'Invalid ion type value'),
    
    ({"context": "moverz",
      "input_item": "LIPIDS",
      "m/z_value": "58",
      "ion_type_value": 'M+H',
      'm/z_tolerance_value': "2"},
      ValueError,
      r'm/z tolerance value outside of range: 0.0001-1'),
])
def test_GenericMWURL_moverz_validate(kwds, error, message):
    with pytest.raises(error, match = message):
        GenericMWURL(kwds)


@pytest.mark.parametrize("kwds, error, message", [
    ({"context": "exactmass",
      "LIPID_abbreviation": "LIPIDS"}, 
      KeyError,
      r"Missing input item\(s\): \{'ion_type_value'\}"),
    
    ({"context": "exactmass",
      "LIPID_abbreviation": "LIPIDS",
      "ion_type_value": 'asdf'},
      ValueError,
      r"Invalid ion type value"),
])
def test_GenericMWURL_exactmass_validate(kwds, error, message):
    with pytest.raises(error, match = message):
        GenericMWURL(kwds)


@pytest.mark.parametrize("kwds, error, message", [
    ({"context": "study",
      "input_item": "study_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid Metabolomics Workbench \(MW\) study ID \(ST<6-digit integer>\)"),
    
    ({"context": "study",
      "input_item": "study_title",
      "input_value": 1,
      "output_item": 'summary',
      'output_format': "txt"},
      ValueError,
      r"Invalid study title \(<string>\)"),
    
    ({"context": "study",
      "input_item": "metabolite_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid Metabolomics Workbench metabolite ID for a study \(ME<6-digit integer>\)"),
    
    ({"context": "compound",
      "input_item": "regno",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid Metabolomics Workbench Metabolite database ID \(<integer>\)"),
    
    ({"context": "compound",
      "input_item": "inchi_key",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid InChIKey \(27-character string\)"),
    
    ({"context": "compound",
      "input_item": "lm_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid LIPID MAPS ID \(LM<2-character LIPID MAPS category><8-10 character string>\)"),
    
    ({"context": "compound",
      "input_item": "pubchem_cid",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid PubChem Compound ID \(<integer>\)"),
    
    ({"context": "compound",
      "input_item": "hmdb_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid Human Metabolome Database ID \(HMDB<integer>\)"),
    
    ({"context": "compound",
      "input_item": "kegg_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid KEGG compound ID \(CO<integer>\)"),
    
    ({"context": "compound",
      "input_item": "chebi_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid ChEBI compound id \(<integer>\)"),
    
    ({"context": "compound",
      "input_item": "metacyc_id",
      "input_value": 1,
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid METACYC compound ID \(<string>\)"),
    
    ({"context": "compound",
      "input_item": "abbrev",
      "input_value": 1,
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid : Lipid bulk abbreviation \(<string>\)"),
    
    ({"context": "protein",
      "input_item": "mgp_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid Human Metabolome Gene/Protein \(MGP\) database gene ID \(MGP<6-digit integer>\)"),
    
    ({"context": "protein",
      "input_item": "gene_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid Entrez gene ID \(<integer>\)"),
    
    ({"context": "protein",
      "input_item": "taxid",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid NCBI taxonomy ID \(<integer>\)"),
    
    ({"context": "protein",
      "input_item": "mrna_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid mRNA ID \(NM_<integer>\)"),
    
    ({"context": "protein",
      "input_item": "refseq_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid mRNA ID \(NP_<integer>\)"),
    
    ({"context": "protein",
      "input_item": "protein_gi",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid NCBI protein GI \(<integer>\)"),
    
    ({"context": "protein",
      "input_item": "uniprot_id",
      "input_value": "AN000001",
      "output_item": 'summary',
      'output_format': "txt"}, 
      ValueError,
      r"Invalid UniProt ID \(see uniprot.org/help/accession_numbers\)"),
])
def test_GenericMWURL_validate_input(kwds, error, message):
    with pytest.raises(error, match = message):
        GenericMWURL(kwds)


def test_MWRESTFile_write():
    class IOError_maker():
        def write(self, text):
            raise IOError
    
    mwfile = MWRESTFile('asdf')
    with pytest.raises(IOError, match = r'"filehandle" parameter must be writable\.'):
        mwfile.write(IOError_maker())


def test_MWRESTFile_is_json():
    mwfile = MWRESTFile('asdf')
    
    assert mwfile._is_json(b'{"key1":1}') == {'key1':1}
    
    with pytest.raises(TypeError, match = r"Expecting <class 'str'> or <class 'bytes'>, but <class 'list'> was passed"):
        mwfile._is_json([])
    
    assert mwfile._is_json('[]]') == False













