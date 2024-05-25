#!/usr/bin/env python

import argparse
import subprocess
import os

IRMA_CMD = "IRMA"
FORWARD_VALUE = "_R1"
REVERSE_VALUE = "_R2"
SINGLE_READS = "single"
PAIRED_READS = "paired"

fastq_valid = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]


def return_files(files_list, valid_ext):
    """Return a list of files that have a valid extension"""
    return_list = [file for file in files_list if file.endswith(tuple(valid_ext))]
    return return_list

    
def pair_files(files_list, first_exp, second_exp):
    """Pair files with illumina scheme"""

    paired_list = []
    first_files = [file for file in files_list if first_exp in file]
    second_files = [file.replace(first_exp, second_exp) for file in first_files]

    for first, second in zip(first_files, second_files):
        if os.path.exists(first) and os.path.exists(second):
            paired_list.append((first, second))
    
    return paired_list

def delete_extension(filename, extensions):
    """Delete a list of extensions"""
    # Ordenamos para que elimine siempre la version larga primero (fastq.gz, por ejemplo)
    extensions.sort(reverse=True)
    for ext in extensions:
        filename = filename.replace(ext, "")
    return filename

def main():

    #argparser

    parser = argparse.ArgumentParser(description="Ejecuta CDC-IRMA en bucle con el módulo especificado")
    parser.add_argument("-m", "--module", help="Modulo que debe usar IRMA")
    parser.add_argument("-r", "--reads", choices=[SINGLE_READS, PAIRED_READS], default=SINGLE_READS,
                        help="Tipo de lecturas introducidas")
    parser.add_argument("--forward_pattern", help="Patrón de la lectura forward (solo para paired)",
                        default=FORWARD_VALUE)
    parser.add_argument("-reverse_pattern", help="Patrón para la lectura reverse (solo en paried)",
                        default=REVERSE_VALUE)
    
    args = parser.parse_args()

    # Listamos los archivos (añadir conf para subcarpeta)
    files = return_files(os.listdir(), fastq_valid)
    # Formateamos:
    # Si Paired:
    if args.reads == PAIRED_READS:
        files = pair_files(files, FORWARD_VALUE, REVERSE_VALUE)

    # Construimos el CMD.
    for file in files:
        cmd = [IRMA_CMD, args.module]
        # Si la lista viene en formato illumina cargamos las dos reads:
        if isinstance(file, tuple):
            cmd.extend(file)
            output_name = delete_extension(file[0], fastq_valid)
            output_name = output_name.replace(FORWARD_VALUE, "")
        else:
            cmd.append(file)
            output_name = delete_extension(file, fastq_valid)
        cmd.append(output_name)
        subprocess.run(cmd)

if __name__ == "__main__":
    main()