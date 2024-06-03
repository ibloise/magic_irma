import os
import re
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import argparse
import magic_irma.constants as constants


SEGMENTS = [
    constants.PB2,
    constants.PB1,
    constants.PA,
    constants.NP,
    constants.HA,
    constants.NA,
    constants.M,
    constants.NS
]

class FluSegment:

    A_TYPE = constants.A_TYPE
    B_TYPE = constants.B_TYPE

    SEGMENTS_METADATA = {
        constants.PB2: {
            "name":  "PB2"
        },
        constants.PB1: {
            "name": "PB1"
        },
        constants.PA: {
            "name": "PA"
        },
        constants.NP: {
            "name": "Nucleoprotein"
        },
        constants.HA: {
            "name": "Hemagglutinin"
        },
        constants.NA: {
            "name": "Neuraminidase"
        },
        constants.M: {
            "name": "Matrix"
        },
        constants.NS: {
            "name": "NS"
        }
    }

    def __init__(self, id, segment, sequence, flu_type: None, segment_type: None) -> None:
        self.id = id
        self.segment = segment
        self.type = flu_type

        if segment in [constants.NA, constants.HA] and flu_type == self.A_TYPE:
            self.segment_type = segment_type
        else:
            self.segment_type = None

        self.sequence = sequence # Property...
        self.length = len(self.sequence)

    def __repr__(self) -> str:
        return f"Flu Segment: {self.segment}-{self.SEGMENTS_METADATA[self.segment]['name']} - Length: {self.length}"

    @staticmethod
    def parse_id(record_id):
        pattern = r'(.*)_([AB])_([A-Za-z0-9]{2,3})_?([HN]\d)?'
        match = re.search(pattern, record_id)
        if match:
            id = match.group(1)
            type = match.group(2)
            segment = match.group(3)
            if match.group(4):
                segment_type = match.group(4)
            else:
                segment_type = None
            return id, type, segment, segment_type
        else:
            raise ValueError("Low information in segment type")

    @classmethod
    def from_record(cls, record):
        id, type, segment, segment_type = cls.parse_id(record.id)
        return cls(id, segment, sequence=record.seq, flu_type=type, segment_type=segment_type)

    @classmethod
    def read_file(cls, file):
        with open(file, "r") as f:
            records = list(SeqIO.parse(f, "fasta"))
        if len(records) > 1:
            raise ValueError("Only single-fasta segment are admitted")
        return cls.from_record(records[0])


class Collection: # Meter como abstract class
    CLASS_TYPE = None

    def __init__(self, *elements) -> None:
        self._collection = []
        self.extend(*elements)

    def append(self, element):
        if isinstance(element, self.CLASS_TYPE):
            self._collection.append(element)
    
    def extend(self, *elements):
        for element in elements:
            self.append(element)
    
    def __iter__(self):
        return iter(self._collection)
    
    def __repr__(self) -> str:
        return f"{self._collection}"
    

class FluSegmentCollection(Collection):

    CLASS_TYPE = FluSegment

    def return_genomes(self):
        return_dict = defaultdict(list)
        for segment in self._collection:
            return_dict[segment.id].append(segment)

        return_collection = FluGenomeCollection()
        for value in return_dict.values():
            flu_genome = FluGenome(*value)
            return_collection.append(flu_genome)
        return return_collection
    
    @classmethod
    def from_fasta(cls, fasta):

        with open(fasta, "r") as file:
            records = list(SeqIO.parse(file, "fasta"))
        
        instance = cls()
        for record in records:
            segment = FluSegment.from_record(record)
            instance.append(segment)
        
        return instance
    
    @classmethod
    def from_list_files(cls, fasta_list):
        ins_global = cls()
        for fasta in fasta_list:
            ins = cls.from_fasta(fasta)
            ins_global.extend(*ins)
        
        return ins_global

class FluGenome:
    A_TYPE = constants.A_TYPE
    B_TYPE = constants.B_TYPE

    def __init__(self, *segments) -> None:
        # Comprobamos que los segmentos tienen un id único:
        self._segments = FluSegmentCollection(*segments)
        self._validate_input()
        self._get_id()
        self._get_type()
        self._index_segments()
        self._get_subtype()

    def _validate_input(self):
        self._validate_id()
        self._validate_type()

    def _get_id(self):
        ids = {segment.id for segment in self._segments}
        self.id = list(ids)[0]

    def _get_type(self):
        types = {segment.type for segment in self._segments}
        self.type = list(types)[0]
    
    def _get_subtype(self):
        if self.type == self.A_TYPE:
            _ha = self._get_segment_type(constants.HA)
            _na = self._get_segment_type(constants.NA)
            self.subtype = f"H{_ha}N{_na}"
        else:
            self.subtype = None

    def _get_segment_type(self, segment):
        subtype = "x"
        if self._check_if_segment(segment):
            subtype = self.segments[segment].segment_type[1]
        return subtype

    def _check_if_segment(self, segment):
        return segment in self.segments

    def _validate_id(self):       
        id = {segment.id for segment in self._segments}
        if len(id) > 1:
            raise ValueError("FluGenome only admitted uniques list of segments")
        
    def _validate_type(self):
        type = {segment.type for segment in self._segments}

        if len(type) > 1:
            raise ValueError('Only one flu type are allowed')
        
    def __repr__(self) -> str:
        return f'Flu:\n {self.id}-Type-{self.type} {self.subtype} \nSegments:\n{self.segments.values()}'
    
    def _index_segments(self):
        self.segments = {segment.segment: segment for segment in self._segments}

    def get_segments(self):
        segments = [segment for segment in self.segments.values()]

        return FluSegmentCollection(*segments)


class FluGenomeSummarizer:
    ID_HEADER = "id"
    TYPE_HEADER = "type"
    SUBTYPE_HEADER = "subtype"
    HEADERS = [ID_HEADER, TYPE_HEADER, SUBTYPE_HEADER]
    SEGMENTS = SEGMENTS
    HEADERS.extend(SEGMENTS)

    @classmethod
    def summary(cls, flu_genome: FluGenome):
        sum_dict = {
            cls.ID_HEADER: flu_genome.id,
            cls.TYPE_HEADER: flu_genome.type,
            cls.SUBTYPE_HEADER: flu_genome.subtype
        }

        for segment in SEGMENTS:
            if segment in flu_genome.segments:
                sum_dict[segment] = flu_genome.segments[segment].length
            else:
                sum_dict[segment] = None
        return sum_dict
    
    @classmethod
    def summary_to_df(cls, flu_genome):
        summary = cls.summary(flu_genome)
        return pd.DataFrame([summary])
    
    @classmethod
    def bulk_summarize(cls, flu_genome_collection):
        return_list = []
        for genome in flu_genome_collection:
            summary = cls.summary(genome)
            return_list.append(summary)
        return return_list
    @classmethod
    def bulk_df_summarize(cls, flu_genome_collection):
        summary = cls.bulk_summarize(flu_genome_collection)
        return pd.DataFrame(summary)


class FluGenomeCollection(Collection):
    CLASS_TYPE = FluGenome

    def to_id_dict(self) -> dict:

        return_dict = {}
        for genome in self._collection:
            segments = genome.get_segments()
            return_dict[genome.id] = segments
        return return_dict

    def to_segment_dict(self):
        return_dict = {segment: FluSegmentCollection() for segment in SEGMENTS}

        for genome in self._collection:
            segments = genome.segments 
            for key, value in segments.items():
                return_dict[key].extend(value)
        return return_dict


class FluGenomeRecordCreator:
    @staticmethod
    def create_record(segment_collection):
        records = []
        for segment in segment_collection:
            id = f"{segment.id}_{segment.type}_{segment.segment}"
            if segment.segment_type:
                id = f"{id}_{segment.segment_type}"
            record = SeqRecord(segment.sequence,
                               id=id, description="")
            records.append(record)
        return records
    
    @staticmethod
    def create_record_dict(segment_dict):
        record_dict = {}
        for key, value in segment_dict.items():
            records = FluGenomeRecordCreator.create_record(value)
            record_dict[key] = records
        return record_dict

class FileWriterInterface:
    def write(self, data, dest_file):
        raise NotImplementedError

class FastaFileWriter(FileWriterInterface):
    def write(self, records, dest_file):
        with open(dest_file, 'w') as file:
            SeqIO.write(records, file, constants.FASTA_FILE)

class FluGenomeFileWriter:
    def __init__(self, file_writer: FileWriterInterface):
        self.file_writer = file_writer

    def write_fasta_from_dict(self, segment_dict, out_dir=""):
        record_dict = FluGenomeRecordCreator.create_record_dict(segment_dict)
        out_dir = out_dir if out_dir else os.getcwd()
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        for key, value in record_dict.items():
            dest_file = os.path.join(out_dir, f"{key}{constants.FASTA_EXT}")
            self.file_writer.write(value, dest_file)