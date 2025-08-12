import click
from pathlib import Path
import logging
from quantmsio.core.idxml import IdXML


@click.group()
def idxml():
    """Commands for IdXML file processing."""
    pass


@idxml.command("convert")
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.argument("output_file", type=click.Path(path_type=Path))
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
def convert_idxml_file(input_file: Path, output_file: Path, verbose: bool):
    """Converts an IdXML file to a PSM parquet file."""
    if verbose:
        logging.basicConfig(level=logging.INFO)

    try:
        # Initialize IdXML parser
        idxml_parser = IdXML(input_file)

        # Get basic information
        psm_count = idxml_parser.get_psm_count()
        protein_count = idxml_parser.get_protein_count()

        click.echo(
            f"Found {psm_count} peptide identifications and {protein_count} proteins"
        )

        # Convert to parquet
        idxml_parser.to_parquet(output_file)

        click.echo(f"Successfully converted {input_file} to {output_file}")
        click.echo(f"Output file size: {output_file.stat().st_size} bytes")

    except Exception as e:
        click.echo(f"Error converting IdXML file: {e}", err=True)
        raise click.Abort()


@idxml.command("info")
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
def info_idxml_file(input_file: Path, verbose: bool):
    """Displays information about an IdXML file."""
    if verbose:
        logging.basicConfig(level=logging.INFO)

    try:
        # Initialize IdXML parser
        idxml_parser = IdXML(input_file)

        # Get basic information
        psm_count = idxml_parser.get_psm_count()
        protein_count = idxml_parser.get_protein_count()

        click.echo(f"IdXML File: {input_file}")
        click.echo(f"Peptide Identifications: {psm_count}")
        click.echo(f"Proteins: {protein_count}")

        # Show sample data
        if psm_count > 0:
            df = idxml_parser.to_dataframe()
            click.echo(f"\nSample data (first 3 rows):")
            click.echo(df.head(3).to_string())

    except Exception as e:
        click.echo(f"Error reading IdXML file: {e}", err=True)
        raise click.Abort()
