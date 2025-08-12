"""
IdXML parser for quantmsio package.
This module provides functionality to parse OpenMS IdXML files and convert them to quantms.io PSM format.
"""

import xml.etree.ElementTree as ET
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from typing import Union, Optional, Dict, List, Tuple
import logging
import re

from quantmsio.core.common import PSM_SCHEMA
from quantmsio.core.format import PSM_FIELDS
from quantmsio.utils.pride_utils import generate_scan_number


class IdXML:
    """
    Parser for OpenMS IdXML files.

    This class provides functionality to parse IdXML files and convert them to quantms.io PSM format.
    IdXML is an OpenMS format for storing peptide and protein identifications.
    """

    def __init__(self, idxml_path: Union[Path, str]):
        """
        Initialize the IdXML parser.

        :param idxml_path: Path to the IdXML file
        """
        self.idxml_path = Path(idxml_path)
        self._protein_map = {}
        self._peptide_identifications = []
        self._parse_protein_identifications()
        self._parse_peptide_identifications()

    def _parse_protein_identifications(self) -> None:
        """
        Parse protein identifications from the IdXML file.
        Extracts protein accessions and their target/decoy status.
        """
        try:
            tree = ET.parse(self.idxml_path)
            root = tree.getroot()

            # Find all ProteinIdentification elements
            for protein_id in root.findall(".//ProteinIdentification"):
                for protein_hit in protein_id.findall(".//ProteinHit"):
                    accession = protein_hit.get("accession", "")
                    if accession:
                        # Check if it's a decoy protein
                        is_decoy = 0
                        target_decoy_param = protein_hit.find(
                            './/UserParam[@name="target_decoy"]'
                        )
                        if target_decoy_param is not None:
                            is_decoy = (
                                1 if target_decoy_param.get("value") == "decoy" else 0
                            )
                        elif accession.startswith("DECOY_"):
                            is_decoy = 1

                        self._protein_map[accession] = {
                            "is_decoy": is_decoy,
                            "accession": accession,
                        }

        except ET.ParseError as e:
            logging.error(f"Error parsing IdXML file: {e}")
            raise
        except Exception as e:
            logging.error(
                f"Unexpected error while parsing protein identifications: {e}"
            )
            raise

    def _parse_peptide_identifications(self) -> None:
        """
        Parse peptide identifications from the IdXML file.
        Extracts peptide hits with their associated information.
        """
        try:
            tree = ET.parse(self.idxml_path)
            root = tree.getroot()

            # Find all PeptideIdentification elements
            for peptide_id in root.findall(".//PeptideIdentification"):
                mz = float(peptide_id.get("MZ", 0))
                rt = float(peptide_id.get("RT", 0))
                spectrum_ref = peptide_id.get("spectrum_reference", "")

                # Extract scan number from spectrum reference
                scan_number = self._extract_scan_number(spectrum_ref)

                for peptide_hit in peptide_id.findall(".//PeptideHit"):
                    peptide_data = self._parse_peptide_hit(
                        peptide_hit, mz, rt, scan_number
                    )
                    if peptide_data:
                        self._peptide_identifications.append(peptide_data)

        except ET.ParseError as e:
            logging.error(f"Error parsing IdXML file: {e}")
            raise
        except Exception as e:
            logging.error(
                f"Unexpected error while parsing peptide identifications: {e}"
            )
            raise

    def _parse_peptide_identifications(self) -> None:
        """
        Parse peptide identifications from the IdXML file.
        Extracts peptide hits with their associated information.
        """
        try:
            tree = ET.parse(self.idxml_path)
            root = tree.getroot()

            # Find all PeptideIdentification elements
            for peptide_id in root.findall(".//PeptideIdentification"):
                mz = float(peptide_id.get("MZ", 0))
                rt = float(peptide_id.get("RT", 0))
                spectrum_ref = peptide_id.get("spectrum_reference", "")

                # Extract scan number from spectrum reference
                scan_number = self._extract_scan_number(spectrum_ref)

                for peptide_hit in peptide_id.findall(".//PeptideHit"):
                    peptide_data = self._parse_peptide_hit(
                        peptide_hit, mz, rt, scan_number
                    )
                    if peptide_data:
                        self._peptide_identifications.append(peptide_data)

        except ET.ParseError as e:
            logging.error(f"Error parsing IdXML file: {e}")
            raise
        except Exception as e:
            logging.error(
                f"Unexpected error while parsing peptide identifications: {e}"
            )
            raise

    def _parse_peptide_hit(
        self, peptide_hit: ET.Element, mz: float, rt: float, scan_number: str
    ) -> Optional[Dict]:
        """
        Parse individual peptide hit information.

        :param peptide_hit: XML element containing peptide hit data
        :param mz: Precursor m/z value
        :param rt: Retention time value
        :param scan_number: Scan number from spectrum reference
        :return: Dictionary containing parsed peptide hit data
        """
        try:
            sequence = peptide_hit.get("sequence", "")
            if not sequence:
                return None

            charge = int(peptide_hit.get("charge", 1))
            score = float(peptide_hit.get("score", 0))

            # Extract protein references
            protein_refs = peptide_hit.get("protein_refs", "")
            protein_accessions = []
            if protein_refs:
                # Parse protein references (comma-separated)
                ref_list = [ref.strip() for ref in protein_refs.split(",")]
                for ref in ref_list:
                    if ref in self._protein_map:
                        protein_accessions.append(self._protein_map[ref]["accession"])

            # Extract modifications
            modifications = self._parse_modifications(sequence)

            # Generate peptidoform notation in MSstats format
            peptidoform = self._generate_peptidoform_msstats(sequence, modifications)

            # Extract additional scores and parameters
            additional_scores = []
            q_value = None
            posterior_error_probability = None
            consensus_support = None

            for user_param in peptide_hit.findall(".//UserParam"):
                param_name = user_param.get("name", "")
                param_value = user_param.get("value", "")

                if param_name == "q-value":
                    try:
                        q_value = float(param_value)
                    except ValueError:
                        pass
                elif param_name == "Posterior Error Probability_score":
                    try:
                        posterior_error_probability = float(param_value)
                    except ValueError:
                        pass
                elif param_name == "consensus_support":
                    try:
                        consensus_support = float(param_value)
                    except ValueError:
                        pass
                else:
                    # Add other scores to additional_scores
                    try:
                        score_value = float(param_value)
                        additional_scores.append(
                            {"score_name": param_name, "score_value": score_value}
                        )
                    except ValueError:
                        pass

            # Determine if this is a decoy peptide
            is_decoy = 0
            target_decoy_param = peptide_hit.find('.//UserParam[@name="target_decoy"]')
            if target_decoy_param is not None:
                is_decoy = 1 if target_decoy_param.get("value") == "decoy" else 0

            # Calculate theoretical m/z (simplified calculation)
            calculated_mz = self._calculate_theoretical_mz(sequence, charge)

            return {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "modifications": modifications,
                "precursor_charge": charge,
                "posterior_error_probability": posterior_error_probability or 0.0,
                "is_decoy": is_decoy,
                "calculated_mz": calculated_mz,
                "observed_mz": mz,
                "additional_scores": additional_scores,
                "mp_accessions": protein_accessions,
                "rt": rt,
                "reference_file_name": self._extract_filename_from_spectrum_ref(
                    scan_number
                ),
                "spectrum_ref": scan_number,
                "q_value": q_value or 0.0,
                "consensus_support": consensus_support or 0.0,
            }

        except Exception as e:
            logging.warning(f"Error parsing peptide hit: {e}")
            return None

    def _parse_modifications(self, sequence: str) -> List[Dict]:
        """
        Parse modifications from peptide sequence.

        :param sequence: Peptide sequence with modification annotations
        :return: List of modification dictionaries
        """
        modifications = []

        # Simple regex to find modifications like (Phospho), (Oxidation), etc.
        mod_pattern = r"\(([^)]+)\)"
        matches = list(re.finditer(mod_pattern, sequence))

        for match in matches:
            mod_name = match.group(1)
            position = match.start()

            # Find the amino acid position this modification applies to
            aa_pos = 0
            for i, char in enumerate(sequence[:position]):
                if char.isalpha():
                    aa_pos += 1

            modifications.append(
                {
                    "modification_name": mod_name,
                    "position": aa_pos,
                    "localization_probability": 1.0,  # Default value
                }
            )

        return modifications

    def _generate_peptidoform_msstats(
        self, sequence: str, modifications: List[Dict]
    ) -> str:
        """
        Generate peptidoform in MSstats notation from sequence and modifications.

        :param sequence: Peptide sequence
        :param modifications: List of modification dictionaries
        :return: Peptidoform in MSstats notation (e.g., PEPTIDE(Oxidation))
        """
        if not modifications:
            return sequence

        # Create a copy of the sequence to modify
        result = list(sequence)

        # Sort modifications by position (descending) to avoid index shifting
        sorted_mods = sorted(modifications, key=lambda x: x["position"], reverse=True)

        for mod in sorted_mods:
            position = mod["position"]
            mod_name = mod["modification_name"]

            # Insert modification after the amino acid at the specified position
            if 0 <= position < len(result):
                result.insert(position + 1, f"({mod_name})")

        return "".join(result)

    def _extract_scan_number(self, spectrum_ref: str) -> str:
        """
        Extract scan number from spectrum reference string.

        :param spectrum_ref: Spectrum reference string
        :return: Extracted scan number
        """
        # Pattern: controllerType=0 controllerNumber=1 scan=37144
        scan_match = re.search(r"scan=(\d+)", spectrum_ref)
        if scan_match:
            return scan_match.group(1)
        return spectrum_ref

    def _extract_filename_from_spectrum_ref(self, spectrum_ref: str) -> str:
        """
        Extract filename from spectrum reference.

        :param spectrum_ref: Spectrum reference string
        :return: Extracted filename
        """
        # For now, return a default filename since IdXML doesn't always contain file info
        return "unknown_file"

    def _calculate_theoretical_mz(self, sequence: str, charge: int) -> float:
        """
        Calculate theoretical m/z for a peptide sequence.
        This is a simplified calculation.

        :param sequence: Peptide sequence
        :param charge: Charge state
        :return: Theoretical m/z value
        """
        # Simplified amino acid masses (monoisotopic)
        aa_masses = {
            "A": 71.03711,
            "R": 156.10111,
            "N": 114.04293,
            "D": 115.02694,
            "C": 103.00919,
            "E": 129.04259,
            "Q": 128.05858,
            "G": 57.02146,
            "H": 137.05891,
            "I": 113.08406,
            "L": 113.08406,
            "K": 128.09496,
            "M": 131.04049,
            "F": 147.06841,
            "P": 97.05276,
            "S": 87.03203,
            "T": 101.04768,
            "W": 186.07931,
            "Y": 163.06333,
            "V": 99.06841,
        }

        # Calculate peptide mass
        peptide_mass = 0
        for aa in sequence:
            if aa in aa_masses:
                peptide_mass += aa_masses[aa]

        # Add N-terminal H and C-terminal OH
        peptide_mass += 1.007825 + 17.00274

        # Calculate m/z
        if charge > 0:
            return (peptide_mass + (charge - 1) * 1.007825) / charge
        else:
            return peptide_mass

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert parsed peptide identifications to a pandas DataFrame.

        :return: DataFrame containing PSM data
        """
        if not self._peptide_identifications:
            return pd.DataFrame()

        df = pd.DataFrame(self._peptide_identifications)

        # Ensure all required columns are present
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
            if col not in df.columns:
                df[col] = None

        return df

    def to_parquet(self, output_path: Union[Path, str]) -> None:
        """
        Convert IdXML data to parquet format and save to file.

        :param output_path: Output file path for parquet file
        """
        df = self.to_dataframe()
        if df.empty:
            logging.warning("No peptide identifications found to convert")
            return

        # Convert to pyarrow table
        table = pa.Table.from_pandas(df)

        # Write to parquet file
        pq.write_table(table, output_path)
        logging.info(f"Successfully converted IdXML to parquet: {output_path}")

    def get_psm_count(self) -> int:
        """
        Get the total number of PSMs found in the IdXML file.

        :return: Number of PSMs
        """
        return len(self._peptide_identifications)

    def get_protein_count(self) -> int:
        """
        Get the total number of proteins found in the IdXML file.

        :return: Number of proteins
        """
        return len(self._protein_map)
