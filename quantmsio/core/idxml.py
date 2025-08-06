import logging
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from typing import Union
from pathlib import Path
from pyopenms import IdXMLFile
from quantmsio.utils.file_utils import close_file
from quantmsio.core.common import PSM_SCHEMA

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)


class IdXML:
    def __init__(self, idxml_file: Union[Path, str]):
        """
        Initialize the IdXML handler.
        :param idxml_file: Path to the IdXML file
        """
        self.idxml_file = idxml_file
        self.peptide_ids = []
        self.protein_ids = []
        self._load_idxml_file()

    def _load_idxml_file(self):
        """
        Load the IdXML file using pyOpenMS.
        """
        IdXMLFile().load(str(self.idxml_file), self.protein_ids, self.peptide_ids)
        logging.info(
            f"Loaded IdXML file with {len(self.peptide_ids)} peptide identifications"
        )

    def iter_psm_table(self, chunksize: int = 1000000):
        """
        Iterate through PSM entries in chunks.
        :param chunksize: Number of PSMs to process in each chunk
        """
        data = []
        for peptide_id in self.peptide_ids:
            # Get spectrum reference
            spectrum_id = (
                peptide_id.getMetaValue("spectrum_reference")
                if peptide_id.hasMetaValue("spectrum_reference")
                else None
            )

            # Process each hit in the peptide identification
            for hit in peptide_id.getHits():
                # Extract information from the hit
                sequence = hit.getSequence().toString()
                aa_sequence = hit.getSequence()

                # Get protein accessions
                protein_accessions = []
                if hit.metaValueExists("protein_references"):
                    protein_accessions = [
                        ref.getProteinID() for ref in hit.getProteinReferences()
                    ]
                elif len(hit.getProteinAccessions()) > 0:
                    protein_accessions = list(hit.getProteinAccessions())

                psm_data = {
                    "sequence": sequence,
                    "peptidoform": sequence,
                    "modifications": [],  # Will be filled later
                    "mp_accessions": protein_accessions,
                    "precursor_charge": hit.getCharge(),
                    "calculated_mz": (
                        aa_sequence.getMZ(hit.getCharge()) if aa_sequence else None
                    ),
                    "observed_mz": peptide_id.getMZ(),
                    "rt": peptide_id.getRT(),
                    "global_qvalue": (
                        hit.getMetaValue("q_value")
                        if hit.metaValueExists("q_value")
                        else None
                    ),
                    "posterior_error_probability": hit.getScore(),
                    "is_decoy": (
                        1
                        if (
                            hit.metaValueExists("target_decoy")
                            and hit.getMetaValue("target_decoy") == "decoy"
                        )
                        else 0
                    ),
                    "reference_file_name": (
                        spectrum_id.split(".")[0] if spectrum_id else None
                    ),
                    "scan": (
                        spectrum_id.split("=")[-1]
                        if spectrum_id and "=" in spectrum_id
                        else spectrum_id
                    ),
                    "additional_scores": [],
                    "cv_params": [],
                    "predicted_rt": None,
                    "ion_mobility": None,
                    "number_peaks": None,
                    "mz_array": None,
                    "intensity_array": None,
                }

                # Add modification information
                # N-terminal modifications
                if aa_sequence.hasNTerminalModification():
                    psm_data["modifications"].append(
                        {
                            "modification_name": aa_sequence.getNTerminalModificationName(),
                            "fields": [
                                {"position": 0, "localization_probability": 1.0}
                            ],
                        }
                    )

                # C-terminal modifications
                if aa_sequence.hasCTerminalModification():
                    psm_data["modifications"].append(
                        {
                            "modification_name": aa_sequence.getCTerminalModificationName(),
                            "fields": [
                                {
                                    "position": len(sequence),
                                    "localization_probability": 1.0,
                                }
                            ],
                        }
                    )

                # Residue modifications
                for i in range(aa_sequence.size()):
                    mod = aa_sequence.getResidueModification(i)
                    if mod:
                        psm_data["modifications"].append(
                            {
                                "modification_name": mod.getFullId(),
                                "fields": [
                                    {"position": i + 1, "localization_probability": 1.0}
                                ],
                            }
                        )

                # Add score information
                if hit.metaValueExists("MS:1001330"):
                    psm_data["additional_scores"].append(
                        {
                            "score_name": "xcorr_score",
                            "score_value": float(hit.getMetaValue("MS:1001330")),
                        }
                    )

                if hit.metaValueExists("MS:1001171"):
                    psm_data["additional_scores"].append(
                        {
                            "score_name": "deltacn",
                            "score_value": float(hit.getMetaValue("MS:1001171")),
                        }
                    )

                data.append(psm_data)

                # Yield data when chunksize is reached
                if len(data) >= chunksize:
                    df = pd.DataFrame(data)
                    yield df
                    data = []

        # Yield remaining data
        if data:
            df = pd.DataFrame(data)
            yield df

    def generate_report(self, chunksize: int = 1000000):
        """
        Generate PSM report in quantms.io format.
        :param chunksize: Number of PSMs to process in each chunk
        """
        for df in self.iter_psm_table(chunksize=chunksize):
            df = self._convert_to_parquet_format(df)
            table = pa.Table.from_pandas(df, schema=PSM_SCHEMA)
            yield table

    def _convert_to_parquet_format(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert dataframe to parquet format.
        :param df: Input dataframe
        :return: Converted dataframe
        """
        # Convert data types to match PSM schema
        if "precursor_charge" in df.columns:
            df["precursor_charge"] = pd.to_numeric(
                df["precursor_charge"], errors="coerce"
            ).astype("Int32")
        if "calculated_mz" in df.columns:
            df["calculated_mz"] = pd.to_numeric(
                df["calculated_mz"], errors="coerce"
            ).astype("float32")
        if "observed_mz" in df.columns:
            df["observed_mz"] = pd.to_numeric(
                df["observed_mz"], errors="coerce"
            ).astype("float32")
        if "posterior_error_probability" in df.columns:
            df["posterior_error_probability"] = pd.to_numeric(
                df["posterior_error_probability"], errors="coerce"
            ).astype("float32")
        if "is_decoy" in df.columns:
            df["is_decoy"] = pd.to_numeric(df["is_decoy"], errors="coerce").astype(
                "Int32"
            )
        if "rt" in df.columns:
            df["rt"] = pd.to_numeric(df["rt"], errors="coerce").astype("float32")
        if "scan" in df.columns:
            df["scan"] = df["scan"].astype(str)
        if "global_qvalue" in df.columns:
            df["global_qvalue"] = pd.to_numeric(
                df["global_qvalue"], errors="coerce"
            ).astype("float32")

        return df

    def write_psm_to_file(
        self,
        output_path: str,
        chunksize: int = 1000000,
    ) -> None:
        """
        Write PSM data to a parquet file.
        :param output_path: Path to the output parquet file
        :param chunksize: Number of PSMs to process in each chunk
        """
        pqwriter = None
        for psm_table in self.generate_report(chunksize=chunksize):
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, psm_table.schema)
            pqwriter.write_table(psm_table)
        close_file(pqwriter=pqwriter)
