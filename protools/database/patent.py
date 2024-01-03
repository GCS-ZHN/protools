from enum import Enum
from collections import namedtuple, OrderedDict
from typing import Union
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils import IUPACData

SequenceListElementTuple = namedtuple('SequenceListElementTuple', ['id', 'is_mandatory'])


class SequenceListElement(Enum):
    """
    Patent sequence list element.
    doc: https://www.wipo.int/export/sites/www/standards/en/pdf/03-25-01.pdf
    """
    APPLICANT_NAMES = SequenceListElementTuple(110, True)
    INVENTION_TITLE = SequenceListElementTuple(120, True)
    FILE_REFERENCE = SequenceListElementTuple(130, True)
    CURRENT_PATENT_APPLICATION = SequenceListElementTuple(140, True)
    CURRENT_FILING_DATE = SequenceListElementTuple(141, True)
    EARLIER_PATENT_APPLICATION = SequenceListElementTuple(150, True)
    EARLIER_FILING_DATE = SequenceListElementTuple(151, True)
    NUMBER_OF_SEQ_IDS = SequenceListElementTuple(160, True)
    SOFTWARE = SequenceListElementTuple(170, False)
    SEQ_ID = SequenceListElementTuple(210, True)
    SEQ_LENGTH = SequenceListElementTuple(211, True)
    SEQ_MOLTYPE = SequenceListElementTuple(212, True)
    SEQ_ORGANISM = SequenceListElementTuple(213, True)
    SEQ_FEATURE = SequenceListElementTuple(220, True)
    SEQ_NAME = SequenceListElementTuple(221, True)
    SEQ_LOCATION = SequenceListElementTuple(222, True)
    OTHER_INFORMATION = SequenceListElementTuple(223, True)
    PUBLICATION = SequenceListElementTuple(300, False)
    AUTHORS = SequenceListElementTuple(301, False)
    TITLE = SequenceListElementTuple(302, False)
    JOURNAL = SequenceListElementTuple(303, False)
    VOLUME = SequenceListElementTuple(304, False)
    ISSUE = SequenceListElementTuple(305, False)
    PAGES = SequenceListElementTuple(306, False)
    DATE = SequenceListElementTuple(307, False)
    DATABASE_ACCESSION = SequenceListElementTuple(308, False)
    DATABASE_ENTRY_DATE = SequenceListElementTuple(309, False)
    DOCUMENT_NUMBER = SequenceListElementTuple(310, False)
    FILING_DATE = SequenceListElementTuple(311, False)
    PUBLICATION_DATE = SequenceListElementTuple(312, False)
    RELEVANT_RESIDUES = SequenceListElementTuple(313, False)
    SEQUENCE = SequenceListElementTuple(400, True)
    UNDEFINED = SequenceListElementTuple(-1, False)

    @staticmethod
    def is_id_style(line: str) -> bool:
        return line.startswith('<') and line.endswith('>') and line[1:-1].isdigit()

    @staticmethod
    def get_element_by_id(id: Union[int, str]) -> 'SequenceListElement':
        for element in SequenceListElement:
            if isinstance(id, str) and f'<{element.value.id}>' == id:
                return element
            elif isinstance(id, int) and element.value.id == id:
                return element
        return SequenceListElement.UNDEFINED


class PatentSequenceList(OrderedDict):
    """
    Patent sequence list object.
    """
    @classmethod
    def from_file(cls, file_path) -> 'PatentSequenceList':
        res = cls()
        with open(file_path, 'r') as f:
            element = SequenceListElement.UNDEFINED
            current_record = None
            element_hook = None
            line_reader = iter(map(lambda x: x.strip().split(), f))

            for lines in line_reader:
                if not lines:
                    continue
                if SequenceListElement.is_id_style(lines[0]):
                    element = SequenceListElement.get_element_by_id(lines[0])
                    print(element)
                    if element != SequenceListElement.UNDEFINED:
                        if element_hook:
                            element_hook()
                            element_hook = None
                    if element == SequenceListElement.SEQ_ID:
                        current_record = SeqRecord(
                            Seq(''),
                            id=lines[1],
                            name='',
                            description='',
                            annotations={},
                            features=[])
                        res[int(lines[1])] = current_record
                    elif element in (
                        SequenceListElement.SEQ_LENGTH,
                        SequenceListElement.SEQ_MOLTYPE,
                        SequenceListElement.SEQ_ORGANISM
                    ):
                        current_record.annotations[element.name.lower()] = lines[1]
                    elif element == SequenceListElement.SEQUENCE:
                        
                        setattr(current_record, 'seq_raw_data', [])
                        def element_hook():
                            seq_data = getattr(current_record, 'seq_raw_data')
                            current_record.seq = Seq(''.join(seq_data))

                else:
                    if element == SequenceListElement.SEQUENCE:
                        seq_data = getattr(current_record, 'seq_raw_data')
                        if current_record.annotations['seq_moltype'] in ('RNA', 'DNA'):
                             # TODO: add support for RNA and DNA
                            pass
                        if current_record.annotations['seq_moltype'] == 'PRT':
                            seq_data.extend(
                                map(lambda x: IUPACData.protein_letters_3to1[x], lines))
                        next(line_reader)
        return res
