"""
Example script demonstrating IdXML to PSM conversion functionality.
"""

from pathlib import Path
from quantmsio.core.idxml import IdXML


def main():
    """Main function demonstrating IdXML usage."""

    idxml_file = "tests/examples/idxml/HF2_8379_RST_1_phos_281118_consensus_fdr_filter_pep_luciphor.idXML"

    if not Path(idxml_file).exists():
        print(f"IdXML file not found: {idxml_file}")
        print("Please update the file path in this script.")
        return

    try:
        print("=== IdXML to PSM Conversion Example ===")
        print(f"Processing file: {idxml_file}")

        idxml = IdXML(idxml_file)

        psm_count = idxml.get_psm_count()
        protein_count = idxml.get_protein_count()

        print(f"\nFile Statistics:")
        print(f"  Total PSMs: {psm_count:,}")
        print(f"  Total proteins: {protein_count:,}")

        if psm_count == 0:
            print("\nNo PSMs found in the file.")
            return

        print("\nConverting to DataFrame...")
        df = idxml.to_dataframe()

        print(f"DataFrame created with shape: {df.shape}")
        print(f"Columns: {list(df.columns)}")

        if not df.empty:
            print(f"\nSample PSM data (first 3 rows):")
            for i, row in df.head(10).iterrows():
                print(f"  PSM {i+1}:")
                print(f"    Sequence: {row['sequence']}")
                print(f"    Peptidoform: {row['peptidoform']}")
                print(f"    Modifications: {row['modifications']}")
                print(f"    Precursor charge: {row['precursor_charge']}")
                print(f"    Posterior error probability: {row['posterior_error_probability']:.6f}")
                print(f"    Is decoy: {row['is_decoy']}")
                print(f"    Calculated m/z: {row['calculated_mz']:.4f}")
                print(f"    Observed m/z: {row['observed_mz']:.4f}")
                print(f"    Additional scores: {row['additional_scores']}")
                print(f"    MP accessions: {row['mp_accessions']}")
                print(f"    RT: {row['rt']:.2f}")
                print(f"    Reference file name: {row['reference_file_name']}")
                print(f"    Q-value: {row['q_value']}")
                print(f"    Cv params: {row['cv_params']}")
                print(f"    Scan: {row['scan']}")
                print()

        output_file = "output_psm.parquet"
        print(f"Converting to parquet format: {output_file}")
        idxml.to_parquet(output_file)

        if Path(output_file).exists():
            file_size = Path(output_file).stat().st_size / 1024 / 1024
            print(f"Successfully created parquet file: {output_file}")
            print(f"File size: {file_size:.2f} MB")
        else:
            print("Error: Parquet file was not created.")

    except Exception as e:
        print(f"Error processing IdXML file: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
