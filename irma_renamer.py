#!/usr/bin/env python
import os
import argparse
import shutil
from Bio import SeqIO


FASTA_FILE = "fasta"
VCF_FILE = "vcf"
BAM_FILE = "bam"

RENAME_PATTERNS = [
    FASTA_FILE,
    VCF_FILE,
    BAM_FILE
]


def extract_sample_id(string, sep="_", pos=0):
    """Devuelve el id de muestra de un esquema idmuestra_loquesea"""
    return string.split(sep)[pos]

# Todas las funciones de renombrado deben tener una interfaz común!
# SI no habrá que retomar la partia function
 
def general_rename(src_file, dest_file):
    print(f"Renombrando {src_file}")
    shutil.copyfile(src_file, dest_file)

def fasta_rename(src_file, dest_file):
    print(f"Renombrado {src_file}")
    with open(src_file, "r") as input_handle, open(dest_file, "w") as output_handle:
        name = os.path.splitext(os.path.basename(dest_file))[0]
        print(name)
        for record in SeqIO.parse(input_handle, "fasta"):
            # Eliminamos la descripción para mayor elegancia del fasta
            record.description = ""
            if record.id in name:
                record.id = name
            else:
                record.id = f"{name}_{record.id}"
            SeqIO.write(record, output_handle, "fasta")

def rebuild_filename(file, preffix, sep="_"):
    filename = os.path.basename(file)
    return f"{preffix}{sep}{filename}"

def create_dest_file(file, dest_folder, preffix, sep="_"):
    filename = rebuild_filename(file, preffix, sep)
    return os.path.join(dest_folder, filename)


def build_files_dict(dir, file_patterns):
    files_dict = {
        file_type: [
            file for file in os.listdir(dir) 
            if file.endswith(file_type)
            ] for file_type in file_patterns
            }
    
    return files_dict

rename_patterns = [FASTA_FILE, VCF_FILE, BAM_FILE]
rename_functions = {
    FASTA_FILE: fasta_rename,
}

def main():
    parser = argparse.ArgumentParser(description="Script para el renombrado flexible de salidas de CDC-RIMA")
    parser.add_argument("-o", "--output_dir", help="Carpeta de salida del renombrado. Por defecto, es el directorio actual", default=os.getcwd())
    parser.add_argument("-f", "--files", help="Tipos de archivo a renombrar, separados por coma. Ej: fasta,vcf", default="fasta")
    args = parser.parse_args()

    output_folder = args.output_dir
    dirs = [folder for folder in os.listdir('.') if os.path.isdir(folder) and folder != output_folder]
    rename_files = [file for file in args.files.lower().split(",") if file in RENAME_PATTERNS]

    if not rename_files:
        print("No se han indicado archivos válidos.")
        print(f"Archivos válidos: {', '.join(RENAME_PATTERNS)}")

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for dir in dirs:
        files_dict = build_files_dict(dir, rename_files)

        for file_type, files in files_dict.items():
            rename_function = rename_functions[file_type] if file_type in rename_functions.keys() else general_rename
            
            list(map(
                lambda file: rename_function(
                    os.path.join(dir,file), create_dest_file(file, output_folder, dir)
                    ), files
                    )
                )

if __name__ == "__main__":
    main()