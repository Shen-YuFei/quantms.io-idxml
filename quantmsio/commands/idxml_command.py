import click
from pathlib import Path
from typing import Union, Optional
from quantmsio.core.project import create_uuid_filename
from quantmsio.core.idxml import IdXML


@click.command(
    "convert-idxml-psm",
    short_help="Convert psm from IdXML to parquet file in quantms.io",
)
@click.option(
    "--idxml_file",
    help="the IdXML file, this will be used to extract the peptide information",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option(
    "--chunksize",
    help="Read batch size",
    default=1000000,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def convert_idxml_psm(
    idxml_file: Union[Path, str],
    output_folder: str,
    chunksize: int,
    output_prefix_file: Optional[str],
) -> None:
    """
    Convert IdXML file to quantms.io PSM parquet format.
    :param idxml_file: the IdXML file, this will be used to extract the peptide information
    :param output_folder: Folder where the parquet file will be generated
    :param chunksize: Read batch size
    :param output_prefix_file: Prefix of the parquet file needed to generate the file name
    """

    if idxml_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = "psm"

    idxml_manager = IdXML(idxml_file)
    output_path = (
        f"{output_folder}/{create_uuid_filename(output_prefix_file, '.psm.parquet')}"
    )
    idxml_manager.write_psm_to_file(output_path=output_path, chunksize=chunksize)
