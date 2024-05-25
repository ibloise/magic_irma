import argparse
import os
import shutil
from Bio import SeqIO
import files

class IrmaRenamer(files.File):
    def __init__(self, filename):
        super().__init__(filename)
        self._build_dest_name()

    def _build_dest_name(self):
        self.dest_file = f"{os.path.basename(self.path)}_{self.name}"
        print(self.dest_file)
    def rename(self):
        print(f"Renombrando {self.full_name}")
        shutil.copyfile(self.full_name, self.dest_file)


class FastaRenamer(IrmaRenamer):
    def rename(self):
        print(f"Renombrado {self.full_name}")
        with open(self.full_name, "r") as input_handle, open(self.dest_file, "w") as output_handle:
            name = os.path.splitext(os.path.basename(self.dest_file))[0]
            for record in SeqIO.parse(input_handle, files.FASTA_FILE):
                # Eliminamos la descripción para mayor elegancia del fasta
                record.description = ""
                if record.id in name:
                    record.id = name
                else:
                    record.id = f"{name}_{record.id}"
                SeqIO.write(record, output_handle, files.FASTA_FILE)


class IrmaDirRenamer:

    def __init__(self, dir, types, renamer_dict, out_dir=None,) -> None:
        self.src_dir = dir
        if not out_dir:
            self.out_dir = os.getcwd()
        else:
            self.out_dir = out_dir
        self.renamer_dict = renamer_dict
        self.types = types
        self._validate_dict()

    @property
    def src_dir(self):
        return self._src_dir

    @src_dir.setter
    def src_dir(self, dir):

        if not os.path.isdir(dir):
            raise ValueError("dir variable must be a directory")

        self._src_dir = dir

    def _validate_dict(self):
        for type in self.types:
            if type not in self.renamer_dict.keys():
                raise ValueError(f"Type {type} not present in renamer dict")
            
    def _bulk_renamer(self, files, renamer: IrmaRenamer):
        for file in files:
            renamer_ins = renamer(file)
            renamer_ins.rename()
    
    def rename_files(self):
        for type in self.types:
            files = [os.path.join(self.src_dir, file) for file in os.listdir(self.src_dir) if file.endswith(type)]
            renamer = self.renamer_dict[type]
            self._bulk_renamer(files, renamer)



def main():

    renamer_dict = {
        files.FASTA_FILE: FastaRenamer,
        files.BAM_FILE: IrmaRenamer,
        files.VCF_FILE: IrmaRenamer
    }

    parser = argparse.ArgumentParser(description="Renombra archivos de IRMA")
    parser.add_argument("-d", "--dirs", help="Directorios que se van a utilizar", nargs="+", required=True)
    parser.add_argument("-t", "--types", nargs="+", choices=[files.FASTA_FILE, files.VCF_FILE, files.BAM_FILE], required=True)
    parser.add_argument("--out_dir", help="Directorio de salida")

    args = parser.parse_args()
    dirs = [dir for dir in args.dirs if os.path.isdir(dir)]
    for dir in dirs:
        try:
            renamer = IrmaDirRenamer(dir, args.types, renamer_dict, args.out_dir)
            renamer.rename_files()
        except ValueError as e:
            print("No se ha podido ejecutar")
            print(e)

if __name__ == "__main__":
    main()

