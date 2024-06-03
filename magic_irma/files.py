import os
import magic_irma.constants as constants

class File:

    def __init__(self, filename, out_dir=None):

        self.name = filename

        self.out_dir = out_dir if out_dir else os.getcwd()

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, filename):
        self._name = os.path.basename(filename)
        self.ext = os.path.splitext(self._name)[1]
        self.path = os.path.dirname(filename)
        self.compressed = self._is_compressed(filename)
        self.full_name = os.path.join(self.path, self._name)

    #Abstract method...
    def _is_compressed(self, filename):
        return False

    def __str__(self) -> str:
        return self.name
    
    def __repr__(self) -> str:
        return self.name


class FastqFile(File):
    FASTQ_EXT = constants.FASTQ_EXT
    FQ_EXT = constants.FQ_EXT
    GZ_EXT = constants.GZ_EXT

    @File.name.setter
    def name(self, filename):
        self._fastq_valid = (
            self.FASTQ_EXT,
            self.FQ_EXT,
            self.GZ_EXT
        )
        if not filename.endswith(self._fastq_valid):
            raise ValueError(f"Filename do not have and appropiate extension. Extensions admitted:{self._fastq_valid}")
        
        self._name = os.path.basename(filename)
        self.path = os.path.dirname(filename)
        self.compressed = self._is_compressed(filename) # Setter!!

    def _is_compressed(self, filename):
        return filename.endswith(self.GZ_EXT)


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
    
    def __repr__(self) -> str:
        return f"{self.common_name}:{self.files}; {self.paired}"
    
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
