import argparse
from . import __version__
from .tools import predict, screen, oracle, validate

TOOLS = [
    ("predict", predict, "Co-fold a single system using Boltz."),
    ("screen", screen, "Co-fold a library using Boltz for virtual screening."),
    ("oracle", oracle, "Use Boltz as an oracle function for single SMILES predictions."),
    ("validate", validate, "Validate Boltz system configuration."),
]

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='boltz-tools',
        description='Collection of helpful tools for performing various co-folding tasks using Boltz.'
    )
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    
    subparsers = parser.add_subparsers(dest='command', required=True)

    for name, module, help_text in TOOLS:
        sp = subparsers.add_parser(name, help=help_text)
        module.add_arguments(sp)
        sp.set_defaults(func=module.main)

    args = parser.parse_args(argv)
    args.func(args)
    