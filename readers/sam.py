import pandas as pd
import re
from base import GenomicDataReader
class SamReader(GenomicDataReader):
    def __init__(self, filepath):
        super().__init__(filepath)
        
    def parse_headers(self):
        """Собирает информацию из строк, начинающихся с @"""
        headers= {'PG': [], 'RG': []}

        with open(self.filename, 'r') as f:
            for line in f:
                if not line.startswith('@'):
                    break

                parts=line.strip().split('\t')
                record_type=parts[0][1:]

                if record_type in ['PG', 'RG']:
                    info = {}
                    for field in parts [1:]:
                        if ':' in field:
                            key, value = field.split (':', 1)
                            info[key]= value 
                    headers[record_type].append(info)
        return headers

    def read_alignments(self):
         """Генератор, который выдает по одному выравниванию"""
         with open(self.filename, 'r') as f:
              for line in f:
                   if line.startswith('@'):
                       continue

                   parts=line.strip().split('\t')
                   if len(parts)<11:
                       continue
                   yield {'qname': parts[0],
                          'rname': parts[2],
                          'pos': int(parts[3])}

    def get_chromosome_stats(self):
        """Статистика: количество выравниваний на хромосому с помощью Pandas"""
        data = list(self.read_alignments())

        if not data:
            return pd.DataFrame(columns=['rname', 'count'])

        df = pd.DataFrame(data)
        stats = df.groupby('rname').size().reset_index(name='count')
        return stats

    def get_region_intersect(self, chrom, start, end):
        """Поиск пересечений"""
        hits = []
        for alignment in self.read_alignments():
            if alignment['rname'] == chrom and start <= alignment['pos'] <= end:
                hits.append(alignment)
        return hits

    def read(self):
        return self.read_alignments()

    def analyze_line(self, line):
        """Разбирает одну строку SAM."""
        parts = line.strip().split('\t')
        if len(parts) >= 11:
            return {'qname': parts[0], 'rname': parts[2], 'pos': int(parts[3])}
        return None

    def get_chromosomes(self):
        """Обязательный метод из GenomicDataReader"""
        chroms = set()
        for aln in self.read_alignments():
            chroms.add(aln['rname'])
        return list(chroms)

    def get_refrence_genome(self):
        """Возвращает имя референсного генома."""
        with open(self.filename, 'r') as f:
            for line in f:
                if line.startswith('@SQ'):
                    parts = line.strip().split('\t')
                    for field in parts:
                        if field.startswith('SN:'):
                            return field.split(':')[1]
                if not line.startswith('@'):
                    break
        return "Unknown"

    def validate_coordinate(self, chrom, pos):
        """Проверяет корректность."""
        valid_chroms = self.get_chromosomes()
        return chrom in valid_chroms and isinstance(pos, int) and pos > 0

