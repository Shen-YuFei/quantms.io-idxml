"""
Test module for IdXML functionality.
"""

import unittest
from pathlib import Path
import tempfile
import shutil
import pandas as pd

from quantmsio.core.idxml import IdXML


class TestIdXML(unittest.TestCase):
    """Test cases for IdXML parser functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_data_dir = Path("tests/examples/idxml")
        self.test_idxml_file = (
            self.test_data_dir
            / "HF2_8379_RST_1_phos_281118_consensus_fdr_filter_pep_luciphor.idXML"
        )

        # Create temporary directory for test outputs
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Clean up test fixtures."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_idxml_initialization(self):
        """Test IdXML parser initialization."""
        if self.test_idxml_file.exists():
            idxml = IdXML(self.test_idxml_file)
            self.assertIsInstance(idxml, IdXML)
            self.assertEqual(idxml.idxml_path, self.test_idxml_file)
        else:
            self.skipTest("Test IdXML file not found")

    def test_protein_parsing(self):
        """Test protein identification parsing."""
        if self.test_idxml_file.exists():
            idxml = IdXML(self.test_idxml_file)
            protein_count = idxml.get_protein_count()
            self.assertGreaterEqual(protein_count, 0)
        else:
            self.skipTest("Test IdXML file not found")

    def test_peptide_parsing(self):
        """Test peptide identification parsing."""
        if self.test_idxml_file.exists():
            idxml = IdXML(self.test_idxml_file)
            psm_count = idxml.get_psm_count()
            self.assertGreaterEqual(psm_count, 0)
        else:
            self.skipTest("Test IdXML file not found")

    def test_dataframe_conversion(self):
        """Test conversion to DataFrame."""
        if self.test_idxml_file.exists():
            idxml = IdXML(self.test_idxml_file)
            df = idxml.to_dataframe()
            self.assertIsInstance(df, pd.DataFrame)

            if not df.empty:
                # Check required columns
                required_columns = [
                    "sequence",
                    "peptidoform",
                    "modifications",
                    "precursor_charge",
                    "posterior_error_probability",
                    "is_decoy",
                    "calculated_mz",
                    "observed_mz",
                    "additional_scores",
                    "mp_accessions",
                    "rt",
                    "reference_file_name",
                    "spectrum_ref",
                    "q_value",
                    "consensus_support",
                ]

                for col in required_columns:
                    self.assertIn(col, df.columns)
        else:
            self.skipTest("Test IdXML file not found")

    def test_parquet_conversion(self):
        """Test conversion to parquet format."""
        if self.test_idxml_file.exists():
            idxml = IdXML(self.test_idxml_file)
            output_path = self.temp_dir / "test_output.parquet"

            # Convert to parquet
            idxml.to_parquet(output_path)

            # Check if file was created
            self.assertTrue(output_path.exists())

            # Check file size
            self.assertGreater(output_path.stat().st_size, 0)
        else:
            self.skipTest("Test IdXML file not found")

    def test_modification_parsing(self):
        """Test modification parsing from peptide sequences."""
        idxml = IdXML(self.test_idxml_file) if self.test_idxml_file.exists() else None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test modification parsing with a sample sequence
        test_sequence = "GSGEKPVSAPGDDT(Phospho)ES(Phospho)LHSQGEEEFDMPQPPHGHVLHR"
        modifications = idxml._parse_modifications(test_sequence)

        self.assertIsInstance(modifications, list)
        self.assertEqual(len(modifications), 2)  # Two phospho modifications

        for mod in modifications:
            self.assertIn("modification_name", mod)
            self.assertIn("position", mod)
            self.assertIn("localization_probability", mod)
            self.assertEqual(mod["modification_name"], "Phospho")

    def test_scan_number_extraction(self):
        """Test scan number extraction from spectrum reference."""
        idxml = IdXML(self.test_idxml_file) if self.test_idxml_file.exists() else None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test scan number extraction
        test_spectrum_ref = "controllerType=0 controllerNumber=1 scan=37144"
        scan_number = idxml._extract_scan_number(test_spectrum_ref)
        self.assertEqual(scan_number, "37144")

        # Test with no scan number
        test_spectrum_ref_no_scan = "controllerType=0 controllerNumber=1"
        scan_number_no_scan = idxml._extract_scan_number(test_spectrum_ref_no_scan)
        self.assertEqual(scan_number_no_scan, test_spectrum_ref_no_scan)

    def test_theoretical_mz_calculation(self):
        """Test theoretical m/z calculation."""
        idxml = IdXML(self.test_idxml_file) if self.test_idxml_file.exists() else None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test with a simple peptide
        test_sequence = "PEPTIDE"
        test_charge = 2

        calculated_mz = idxml._calculate_theoretical_mz(test_sequence, test_charge)
        self.assertIsInstance(calculated_mz, float)
        self.assertGreater(calculated_mz, 0)


if __name__ == "__main__":
    unittest.main()
