# from .birch import birch
import click as original_click
import rich_click as click
from rich.traceback import install

install(
    suppress=[original_click, click],
)

@click.group('aitw-darcy')
def cmd_root():
    """ Command line interface for AITW Darcy simulations """
    pass

@cmd_root.command()
@click.argument('json_file', required=True, type=click.Path(exists=True))
@click.option('--output_dir', type=click.Path(), help='Output directory')
def infiltration(
        json_file: str,
        output_dir: str = '.',
    ):
    """ Run the infiltration simulation command line interface """
    from aitw_darcy.data import Params
    from aitw_darcy.infiltration import run_simulation

    params = Params.from_json(json_file)
    run_simulation(params, output_dir=output_dir)
