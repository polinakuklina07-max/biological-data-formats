from base import SequenceReader

class FastaReader(SequenceReader):
    """Класс для чтения и анализа файлов формата FASTA."""

    def read(self):
        if self._file is None or self._file.closed:
            self._file = open(str(self.filename), 'r', encoding='utf-8')
        self._file.seek(0)
        
        current_ID = None
        current_seq = []
        
        for line in self._file:
            line = line.strip()
            if line == "":
                continue
                
            if line.startswith(">"):
                if current_ID is not None:
                    full_sequence = "".join(current_seq)
                    yield {"ID": current_ID, "sequence": full_sequence}
                current_ID = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
                
        if current_ID is not None:
            full_sequence = "".join(current_seq)
            yield {"ID": current_ID, "sequence": full_sequence}

    def analyze_line(self, line):
        """Реализация обязательного абстрактного метода"""
        pass

    def get_sequence_ID(self):
        ID = []
        for part in self.read():
            ID.append(part["ID"])
        return ID

    def validate_sequence(self, sequence):
        return len(sequence) > 0

    def get_sequence_count(self):
        count = 0
        for i in self.read():
            count += 1
        return count

    def get_average_length(self):
        total_length = 0
        count = 0
        for item in self.read():
            total_length += len(item["sequence"])
            count += 1
        if count > 0:
            return total_length / count
        else:
            return 0.0
