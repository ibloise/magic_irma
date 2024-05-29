import argparse
from magic_irma import constants
from magic_irma.bulk import IrmaBulk

def main():

    parser = argparse.ArgumentParser(description="Rapid tools for deal with CDC-IRMA batch analysis")
    parser.add_argument("--irma_cmd", default="IRMA")


    subparser = parser.add_subparsers(dest="command")

    bulk_parser = subparser.add_parser("bulker")

    bulk_parser.add_argument("-f", "--files", help="Archivos que se van a utilizar", nargs="+", required=True)
    bulk_parser.add_argument("-m", "--module", help="Modulo que debe usar IRMA", required=True)
    bulk_parser.add_argument("-r", "--reads", choices=[constants.SINGLE_READS, constants.PAIRED_READS], default=constants.SINGLE_READS,
                        help="Tipo de lecturas introducidas")
    bulk_parser.add_argument("--forward_pattern", help="Patrón de la lectura forward (solo para paired)",
                        default=constants.FORWARD_VALUE)
    bulk_parser.add_argument("--reverse_pattern", help="Patrón para la lectura reverse (solo en paried)",
                        default=constants.REVERSE_VALUE)
    

    args = parser.parse_args()

    if args.command == "bulker":
        bulker_main(args)
    

def bulker_main(args):

    #argparser

    # Declaramos el bulker:

    if not args.files:
        raise ValueError("Es necesario pasar archivos")
    
    irma_bulker = IrmaBulk(args.module, args.files, args.reads, args.irma_cmd, args.forward_pattern, args.reverse_pattern)

    irma_bulker.run()