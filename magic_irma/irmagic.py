import argparse
import os 
import magic_irma.constants as constants
from magic_irma.bulk import IrmaBulk
from magic_irma.renamer import FastaRenamer, IrmaRenamer, IrmaDirRenamer
import magic_irma.flu_tools as flu

def main():

    parser = argparse.ArgumentParser(description="Rapid tools for deal with CDC-IRMA batch analysis")
    
    subparser = parser.add_subparsers(dest="command")

    bulk_parser = subparser.add_parser("bulk", description="Ejecuta IRMA en bucle con el módulo especificado")

    bulk_parser.add_argument("-f", "--files", help="Archivos que se van a utilizar", nargs="+", required=True)
    bulk_parser.add_argument("-m", "--module", help="Modulo que debe usar IRMA", required=True)
    bulk_parser.add_argument("-r", "--reads", choices=[constants.SINGLE_READS, constants.PAIRED_READS], default=constants.SINGLE_READS,
                        help="Tipo de lecturas introducidas")
    bulk_parser.add_argument("--forward_pattern", help="Patrón de la lectura forward (solo para paired)",
                        default=constants.FORWARD_VALUE)
    bulk_parser.add_argument("--reverse_pattern", help="Patrón para la lectura reverse (solo en paried)",
                        default=constants.REVERSE_VALUE)
    bulk_parser.add_argument("--irma_cmd", default="IRMA", required=False)
    
    rename_parser = subparser.add_parser('rename')

    rename_parser.add_argument("-d", "--dirs", help="Directorios que se van a utilizar", nargs="+", required=True)
    rename_parser.add_argument("-t", "--types", nargs="+", choices=[constants.FASTA_FILE, constants.VCF_FILE, constants.BAM_FILE], required=True)
    rename_parser.add_argument("--out_dir", help="Directorio de salida")

    flu_parser = subparser.add_parser('flu', description="Tools for working with flu Fasta files from irma")

    flu_parser.add_argument("-f", "--files", help="Archivos que se van a utilizar", nargs="+", required=True)
    flu_parser.add_argument("-o", "--out_dir")
    flu_parser.add_argument("-m", "--mode", choices=["samples", "segments"], help= "Modo de agrupación")
    flu_parser.add_argument("--summary", action="store_true")


    args = parser.parse_args()

    if args.command == "bulk":
        bulker_main(args)
    elif args.command == "rename":
        rename_main(args)
    elif args.command == "flu":
        flu_main(args)

    
def flu_main(args):

    collection = flu.FluSegmentCollection.from_list_files(args.files)
    genome_collection = collection.return_genomes()


    if args.mode == "samples":
        segments = genome_collection.to_id_dict()
    else:
        segments = genome_collection.to_segment_dict()

    file_writer = flu.FastaFileWriter()
    flu_writer = flu.FluGenomeFileWriter(file_writer)
    flu_writer.write_fasta_from_dict(segments, args.out_dir)

    if args.summary:
        summary = flu.FluGenomeSummarizer.bulk_df_summarize(genome_collection)
        summary.to_csv('summary.csv')


def rename_main(args):

    renamer_dict = {
        constants.FASTA_FILE: FastaRenamer,
        constants.BAM_FILE: IrmaRenamer,
        constants.VCF_FILE: IrmaRenamer
    }

    dirs = [dir for dir in args.dirs if os.path.isdir(dir)]
    for dir in dirs:
        try:
            renamer = IrmaDirRenamer(dir, args.types, renamer_dict, args.out_dir)
            renamer.rename_files()
        except ValueError as e:
            print("No se ha podido ejecutar")
            print(e)

def bulker_main(args):

    #argparser

    # Declaramos el bulker:

    if not args.files:
        raise ValueError("Es necesario pasar archivos")
    
    irma_bulker = IrmaBulk(args.module, args.files, args.reads, args.irma_cmd, args.forward_pattern, args.reverse_pattern)

    irma_bulker.run()