"""Tests for mwtab.mwtab module, specifically the download_results_file method."""

import pytest
from mwtab.mwtab import MWTabFile
from unittest.mock import MagicMock


class TestDownloadResultsFile:
    """Test suite for MWTabFile.download_results_file() method."""
    
    def _create_mwfile_with_ids(self, st_id="ST000001", an_id="AN000001"):
        """Helper to create a MWTabFile with study and analysis IDs."""
        mwfile = MWTabFile("test_source")
        mwfile["METABOLOMICS WORKBENCH"] = {
            "STUDY_ID": st_id,
            "ANALYSIS_ID": an_id,
            "VERSION": "1",
            "CREATED_ON": "2024-01-01"
        }
        return mwfile
    
    def test_download_results_file_ms_success(self, tmp_path, mocker):
        """Test successful download of MS results file."""
        mwfile = self._create_mwfile_with_ids("ST000001", "AN000001")
        mwfile["MS"] = {
            "MS_RESULTS_FILE": {
                "filename": "ST000001_AN000001_Results.txt",
                "UNITS": "Peak area",
                "Has m/z": "Yes",
                "Has RT": "Yes",
                "RT units": "Minutes"
            }
        }
        
        mock_response = MagicMock()
        mock_response.content = b"mock results data content"
        mock_get = mocker.patch('requests.Session.get', return_value=mock_response)
        
        result = mwfile.download_results_file(output_dir=str(tmp_path))
        
        expected_path = str(tmp_path / "ST000001_AN000001_Results.txt")
        assert result == expected_path
        mock_get.assert_called_once_with(
            "https://www.metabolomicsworkbench.org/studydownload/ST000001_AN000001_Results.txt"
        )
        
        with open(expected_path, 'rb') as f:
            assert f.read() == b"mock results data content"
    
    def test_download_results_file_nmr_success(self, tmp_path, mocker):
        """Test successful download of NMR results file."""
        mwfile = self._create_mwfile_with_ids("ST000002", "AN000002")
        mwfile["NM"] = {
            "NMR_RESULTS_FILE": {
                "filename": "ST000002_AN000002_Results.txt"
            }
        }
        
        mock_response = MagicMock()
        mock_response.content = b"nmr results data"
        mock_get = mocker.patch('requests.Session.get', return_value=mock_response)
        
        result = mwfile.download_results_file(output_dir=str(tmp_path))
        
        expected_path = str(tmp_path / "ST000002_AN000002_Results.txt")
        assert result == expected_path
        mock_get.assert_called_once_with(
            "https://www.metabolomicsworkbench.org/studydownload/ST000002_AN000002_Results.txt"
        )
    
    def test_download_results_file_no_results_info(self, tmp_path):
        """Test that None is returned when no results file info exists."""
        mwfile = self._create_mwfile_with_ids("ST000003", "AN000003")
        
        result = mwfile.download_results_file(output_dir=str(tmp_path))
        assert result is None
    
    def test_download_results_file_missing_study_id(self, tmp_path):
        """Test that ValueError is raised when study_id is missing."""
        mwfile = MWTabFile("test")
        mwfile["METABOLOMICS WORKBENCH"] = {"ANALYSIS_ID": "AN000001"}
        mwfile["MS"] = {"MS_RESULTS_FILE": {"filename": "test.txt"}}
        
        with pytest.raises(ValueError, match="study_id"):
            mwfile.download_results_file(output_dir=str(tmp_path))
    
    def test_download_results_file_missing_analysis_id(self, tmp_path):
        """Test that ValueError is raised when analysis_id is missing."""
        mwfile = MWTabFile("test")
        mwfile["METABOLOMICS WORKBENCH"] = {"STUDY_ID": "ST000001"}
        mwfile["MS"] = {"MS_RESULTS_FILE": {"filename": "test.txt"}}
        
        with pytest.raises(ValueError, match="analysis_id"):
            mwfile.download_results_file(output_dir=str(tmp_path))
    
    def test_download_results_file_uses_provided_session(self, tmp_path, mocker):
        """Test that provided session is used instead of creating new one."""
        mwfile = self._create_mwfile_with_ids("ST000004", "AN000004")
        mwfile["MS"] = {"MS_RESULTS_FILE": {"filename": "test.txt"}}
        
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.content = b"data"
        mock_session.get.return_value = mock_response
        
        result = mwfile.download_results_file(output_dir=str(tmp_path), session=mock_session)
        
        mock_session.get.assert_called_once_with(
            "https://www.metabolomicsworkbench.org/studydownload/ST000004_AN000004_Results.txt"
        )
        assert result is not None
    
    def test_download_results_file_http_error(self, tmp_path, mocker):
        """Test that HTTP errors are properly raised."""
        mwfile = self._create_mwfile_with_ids("ST000005", "AN000005")
        mwfile["MS"] = {"MS_RESULTS_FILE": {"filename": "test.txt"}}
        
        mock_response = MagicMock()
        mock_response.raise_for_status.side_effect = Exception("404 Not Found")
        mock_get = mocker.patch('requests.Session.get', return_value=mock_response)
        
        with pytest.raises(Exception, match="404"):
            mwfile.download_results_file(output_dir=str(tmp_path))
    
    def test_download_results_file_default_filename(self, tmp_path, mocker):
        """Test that default filename is used when results file dict has no filename key."""
        mwfile = self._create_mwfile_with_ids("ST000006", "AN000006")
        mwfile["MS"] = {"MS_RESULTS_FILE": {"UNITS": "Peak area"}}
        
        mock_response = MagicMock()
        mock_response.content = b"data"
        mock_get = mocker.patch('requests.Session.get', return_value=mock_response)
        
        result = mwfile.download_results_file(output_dir=str(tmp_path))
        
        expected_path = str(tmp_path / "ST000006_AN000006_Results.txt")
        assert result == expected_path