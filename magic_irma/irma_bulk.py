#!/usr/bin/env python
import argparse
import subprocess
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

IRMA_CMD = "IRMA"
FORWARD_VALUE = "_R1"
REVERSE_VALUE = "_R2"
SINGLE_READS = "single"
PAIRED_READS = "paired"


class FastqFile:
    FASTQ_EXT = ".fastq"
    FQ_EXT = ".fq"
    GZ_EXT = ".gz"

    def __init__(self, filename):

        self._fastq_valid = (self.FASTQ_EXT, 
                             f"{self.FASTQ_EXT}{self.GZ_EXT}",
                             self.FQ_EXT,
                             f"{self.FQ_EXT}{self.GZ_EXT}")
        self.name = filename
        
    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, filename):
        if not filename.endswith(self._fastq_valid):
            raise ValueError(f"Filename do not have and appropiate extension. Extensions admitted:{self._fastq_valid}")
        
        self._name = os.path.basename(filename)
        self.path = os.path.dirname(filename)
        self.compressed = filename.endswith(self.GZ_EXT)

    def __str__(self) -> str:
        return self.name
    
    def __repr__(self) -> str:
        return self.name

class FastqSet:
    def __init__(self, file_tuple, paired, patterns) -> None:

        self.files = file_tuple
        self.paired = paired
        self.patterns = patterns
        check_file = file_tuple[0]
        self.compress = check_file.compressed
        if paired:
            common_name = self._extract_common_id(check_file.name, patterns[0], check_file.compressed)
        else:
            common_name = self._extract_common_id(check_file.name, compress=check_file.compressed)
        self.common_name = common_name

    def _extract_common_id(self, filename, pattern = None, compress = None):
        # Regex are not accept!
        name = os.path.splitext(filename)[0]

        if pattern:
            name = name.replace(pattern, "")
        
        # Si el archivo era un fq comprimido, tiene "dos extensiones"
        if compress:
            name = os.path.splitext(name)[0]
        return name
    
class FastqCollection:

    def __init__(self) -> None:
        self._collection = []

    @classmethod
    def from_list(cls, **kwargs):
        # Debe recibir los argumentos de append_files
        collection = cls()
        collection.append_files(**kwargs)
        return collection
    
    def append_files(self, files, paired: bool,
                     first_pattern: str = "_R1",
                     second_pattern: str = "_R2"):
        clean_files = [
            self._transform_file(file) for file
            in files if self._transform_file(file)
        ]
        if paired:
            files_tuple_list = self._pair_files(clean_files,
                                            first_pattern,
                                            second_pattern)
        else:
            # Convertimos en tuple para estandarizar tipos
            files_tuple_list = [(file,)for file in clean_files]

        for tuple in files_tuple_list:
            fastq_set = FastqSet(
                tuple,
                paired,
                (first_pattern, second_pattern)
            )
            self._collection.append(fastq_set)

    def _transform_file(self, file):
        try:
            fastq_file = FastqFile(file)
        except ValueError as e:
            print(f"Error: {file}: {e}")
            fastq_file = None
        return fastq_file
    
    def _pair_files(self, files_list, first_pattern, second_pattern):

        paired_list = []
        first_files = [file for file in files_list if first_pattern in file.name]
        second_files = [file for file in files_list if second_pattern in file.name]

        for first, second in zip(first_files, second_files):
                paired_list.append((first, second))
    
        return paired_list

    def iter_files(self):
        # Método que itera sobre las tuplas de archivos y las devuelve
        return iter(self._collection)

    def __repr__(self) -> str:
        return f"{self._collection}"

class IrmaBulk:

    SINGLE_READS = "single"
    PAIRED_READS = "paired"

    def __init__(self, 
                 module, file_list, reads_type, irma_cmd = "IRMA",
                 first_pattern = "_R1",
                 second_pattern = "_R2") -> None:
        
        self.irma_cmd = irma_cmd
        self.module = module
        self._cmd_list = []
        #Lógica de comprobación de módulos
        
        if reads_type == self.SINGLE_READS:
            paired = False
        elif reads_type == self.PAIRED_READS:
            paired = True
        else:
            raise ValueError(f"Reads type only addmit {self.SINGLE_READS} or {self.PAIRED_READS} values")
        
        fastq_collection = FastqCollection.from_list(files = file_list, paired=paired, first_pattern = first_pattern, 
                                                     second_pattern = second_pattern)
        self._extend_cmd_list(fastq_collection, module)


    def _extend_cmd_list(self, fastq_collection, module = None):
        for fastq_set in fastq_collection.iter_files():
            cmd = self._build_cmd(
                files = fastq_set.files,
                output= fastq_set.common_name,
                module= module
            )
            self._cmd_list.append(cmd)

    def _build_cmd(self, files, module = None , output = None):
        #FIles debe ser tuple (dar typing)
        if not module:
            module = self.module
        cmd = [self.irma_cmd, module]
        files = [file.name for file in files if file]
        cmd.extend(files)   
        if output:
            cmd.append(output)
        return cmd
    
    def add_execution(self, file_list, paired: bool, module=None, first_pattern = "_R1", second_pattern = "_R2"):
        fastq_collection = FastqCollection.from_list(files = file_list,paired = paired,
                                                     first_pattern = first_pattern, second_pattern = second_pattern)
        self._extend_cmd_list(fastq_collection, module)

    def run(self):
        for cmd in self._cmd_list:
            subprocess.run(cmd)



def main():

    #argparser

    parser = argparse.ArgumentParser(description="Ejecuta CDC-IRMA en bucle con el módulo especificado")
    parser.add_argument("-f", "--files", help="Archivos que se van a utilizar", nargs="+")
    parser.add_argument("-m", "--module", help="Modulo que debe usar IRMA")
    parser.add_argument("-r", "--reads", choices=[SINGLE_READS, PAIRED_READS], default=SINGLE_READS,
                        help="Tipo de lecturas introducidas")
    parser.add_argument("--forward_pattern", help="Patrón de la lectura forward (solo para paired)",
                        default=FORWARD_VALUE)
    parser.add_argument("--reverse_pattern", help="Patrón para la lectura reverse (solo en paried)",
                        default=REVERSE_VALUE)
    parser.add_argument("--irma_cmd", default=IRMA_CMD)
    
    args = parser.parse_args()

    # Declaramos el bulker:

    if not args.files:
        raise ValueError("Es necesario pasar archivos")
    
    irma_bulker = IrmaBulk(args.module, args.files, args.reads, args.irma_cmd, args.forward_pattern, args.reverse_pattern)

    irma_bulker.run()

if __name__ == "__main__":
    main()