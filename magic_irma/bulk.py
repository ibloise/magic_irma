#!/usr/bin/env python
import argparse
import subprocess
from files import FastqCollection

IRMA_CMD = "IRMA"
FORWARD_VALUE = "_R1"
REVERSE_VALUE = "_R2"
SINGLE_READS = "single"
PAIRED_READS = "paired"


# Hay que meterle un log a esto
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
        
        print(fastq_collection)
        self._extend_cmd_list(fastq_collection, module)
        print(self._cmd_list)


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
    parser.add_argument("-f", "--files", help="Archivos que se van a utilizar", nargs="+", required=True)
    parser.add_argument("-m", "--module", help="Modulo que debe usar IRMA", required=True)
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