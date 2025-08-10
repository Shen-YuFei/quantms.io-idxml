# IdXML to PSM Conversion

This document describes how to use the IdXML to PSM conversion functionality in quantms.io.

## Overview

The IdXML parser allows you to convert OpenMS IdXML files to the standardized quantms.io PSM format. IdXML is an OpenMS format for storing peptide and protein identifications, and this tool converts it to a parquet format that can be used for downstream analysis.

## Features

- Parse OpenMS IdXML files
- Extract peptide and protein identifications
- Convert to quantms.io PSM format
- Support for modifications (e.g., phosphorylation, oxidation)
- Export to parquet format
- Command-line interface for easy integration

## Installation

The IdXML functionality is included in the quantms.io package. No additional dependencies are required.

## Usage

### Command Line Interface

#### Convert IdXML to PSM

```bash
quantmsio convert-idxml --idxml_file data.identifications.idXML --output_folder ./output
```

**Options:**

- `--idxml_file`: Path to the IdXML file (required)
- `--output_folder`: Output directory for the parquet file (required)
- `--output_prefix_file`: Prefix for the output filename (optional, default: "psm")
- `--verbose`: Enable verbose logging (optional)

**Example:**

```bash
quantmsio convert-idxml \
  --idxml_file HF2_8379_RST_1_phos_281118_consensus_fdr_filter_pep_luciphor.idXML \
  --output_folder ./output \
  --output_prefix_file phospho_psm \
  --verbose
```

#### Get IdXML File Information

```bash
quantmsio info-idxml --idxml_file data.identifications.idXML
```

This command displays detailed information about the IdXML file, including:

- Total number of PSMs
- Total number of proteins
- Column information
- Sample data

### Programmatic Usage

#### Basic Conversion

```python
from quantmsio.core.idxml import IdXML
from pathlib import Path

# Initialize parser
idxml = IdXML("data.identifications.idXML")

# Get statistics
psm_count = idxml.get_psm_count()
protein_count = idxml.get_protein_count()
print(f"Found {psm_count} PSMs and {protein_count} proteins")

# Convert to DataFrame
df = idxml.to_dataframe()
print(f"DataFrame shape: {df.shape}")

# Convert to parquet
idxml.to_parquet("output.parquet")
```

#### Accessing Parsed Data

```python
# Get protein information
protein_map = idxml._protein_map
for accession, info in protein_map.items():
    print(f"Protein: {accession}, Is Decoy: {info['is_decoy']}")

# Get peptide identifications
peptide_data = idxml._peptide_identifications
for peptide in peptide_data[:5]:  # First 5 peptides
    print(f"Sequence: {peptide['sequence']}")
    print(f"Charge: {peptide['precursor_charge']}")
    print(f"Modifications: {peptide['modifications']}")
    print("---")
```

## Input Format

The IdXML parser expects OpenMS IdXML files with the following structure:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<IdXML version="1.5">
    <ProteinIdentification>
        <ProteinHit id="PH_1" accession="sp|P12345|PROTEIN" score="0.0">
            <UserParam type="string" name="target_decoy" value="target"/>
        </ProteinHit>
    </ProteinIdentification>

    <PeptideIdentification score_type="score" MZ="709.47" RT="6534.87"
                          spectrum_reference="controllerType=0 controllerNumber=1 scan=37144">
        <PeptideHit score="-1.0" sequence="GSGEKPVSAPGDDT(Phospho)ES(Phospho)LHSQGEEEFDMPQPPHGHVLHR"
                    charge="6" protein_refs="PH_1">
            <UserParam type="string" name="target_decoy" value="target"/>
            <UserParam type="float" name="q-value" value="2.49e-04"/>
        </PeptideHit>
    </PeptideIdentification>
</IdXML>
```

## Output Format

The converted PSM data includes the following fields:

- `sequence`: Peptide sequence
- `peptidoform`: Modified peptide sequence in MSstats notation
- `modifications`: List of modifications with positions
- `precursor_charge`: Charge state
- `posterior_error_probability`: PEP score
- `is_decoy`: Decoy indicator (0=target, 1=decoy)
- `calculated_mz`: Theoretical m/z
- `observed_mz`: Experimental m/z
- `additional_scores`: Additional scoring information
- `mp_accessions`: Protein accessions
- `rt`: Retention time
- `reference_file_name`: Source file name
- `spectrum_ref`: Spectrum reference/scan number
- `q_value`: Q-value
- `consensus_support`: Consensus support score

## Supported Modifications

The parser recognizes common modifications in parentheses:

- `(Phospho)`: Phosphorylation
- `(Oxidation)`: Oxidation
- `(Acetyl)`: Acetylation
- `(Carbamidomethyl)`: Carbamidomethylation

## Error Handling

The parser includes comprehensive error handling:

- XML parsing errors are caught and reported
- Missing or malformed data is handled gracefully
- Warnings are logged for problematic entries
- The process continues even if individual peptides fail to parse

## Performance

- Memory-efficient parsing using ElementTree
- Batch processing for large files
- Optimized data structures for fast access
- Minimal memory footprint during conversion

## Limitations

- Currently supports basic modification parsing
- Theoretical m/z calculation is simplified
- File name extraction from spectrum references may be limited
- Some advanced IdXML features may not be fully supported

## Troubleshooting

### Common Issues

1. **File not found**: Ensure the IdXML file path is correct
2. **Permission denied**: Check file and directory permissions
3. **Memory errors**: For very large files, consider processing in chunks
4. **XML parsing errors**: Verify the IdXML file is valid and not corrupted

### Debug Mode

Enable verbose logging to see detailed parsing information:

```bash
quantmsio convert-idxml --idxml_file data.idXML --output_folder ./output --verbose
```

## Examples

### Basic Conversion

```bash
# Convert a single IdXML file
quantmsio convert-idxml --idxml_file experiment.idXML --output_folder ./results
```

### Batch Processing

```bash
# Process multiple files
for file in *.idXML; do
    quantmsio convert-idxml --idxml_file "$file" --output_folder ./results --output_prefix_file "${file%.idXML}"
done
```

### Integration with Workflows

```python
import subprocess
from pathlib import Path

def convert_idxml_files(input_dir, output_dir):
    """Convert all IdXML files in a directory."""
    input_path = Path(input_dir)
    output_path = Path(output_dir)

    for idxml_file in input_path.glob("*.idXML"):
        output_file = output_path / f"{idxml_file.stem}.parquet"

        cmd = [
            "quantmsio", "convert-idxml",
            "--idxml_file", str(idxml_file),
            "--output_folder", str(output_path),
            "--output_prefix_file", idxml_file.stem
        ]

        subprocess.run(cmd, check=True)
        print(f"Converted {idxml_file} to {output_file}")

# Usage
convert_idxml_files("./input", "./output")
```

## Contributing

To contribute to the IdXML functionality:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Support

For issues and questions:

- Check the existing documentation
- Review the test files for examples
- Open an issue on the GitHub repository
- Contact the development team

## License

This functionality is part of the quantms.io package and follows the same license terms.
