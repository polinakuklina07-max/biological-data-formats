from base import SequenceReader
class FastaReader(SequenceReader):
    """Класс для чтения и анализа файлов формата FASTA."""

    def read(self):
        self.close()
        self._file = open(str(self.filename))
        current_ID = None
        current_seq = []

        for line in self._file:
            line = line.strip()
            if line == "":
                continue
            if line.startswith(">"):
                if current_ID is not None:
                    full_sequence = ""
                    for chunk in current_seq:
                        full_sequence += chunk
                    yield {"ID": current_ID, "sequence": full_sequence}
                    current_ID = line[1:]
                    current_seq = []
            else:
                current_seq.append(line)
        if current_ID is not None:
            full_sequence = ""
            for chunk in current_seq:
                full_sequence += chunk
            yield {"ID": current_ID, "sequence": full_sequence}
        self.close()

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
        for lenght in self.read():
            total_length += len(lenght["sequence"])
            count += 1
        if count > 0:
            return total_length / count
        else:
            return 0.0
        
