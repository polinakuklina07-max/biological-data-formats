from abc import ABC, abstractmethod
from typing import Generator, Any

class Reader(ABC):
    """Абстрактный класс для чтения биологических файлов."""

    def __init__(self, filename: str):
        self.filename = filename
        self._file = None

    @abstractmethod
    def read(self):
        """Абстрактный метод генератор."""
        pass

    def close(self):
        """Закрывает открытый файл."""
        if self.__file and not self.__file.closed:
            self._file.close()

    @abstractmethod
    def analyze_line(self, line):
        """Метод для анализа и обработки отдельной строки файла."""
        pass
class SequenceReader(Reader, ABC):
    """Абстрактный класс для последовательностей FASTA и FASTQ."""

    @abstractmethod
    def get_sequence_ID(self):
        """Возвращает список ID всех последовательностей."""
        pass

    def validate_sequence(self, sequence):
        """Проверка последовательности на корректность."""
        pass

class GenomicDataReader(Reader, ABC):
    """Абстрактный класс для геномных данных формата SAM и VCF."""

    @abstractmethod
    def get_chromosomes(self):
        """Возвращает список хромосом."""
        pass

    @abstractmethod
    def get_refrence_genome(self):
        """Возвращает имя референсного генома."""
        pass

    @abstractmethod
    def validate_coordinate(self, chrom, pos):
        """Проверяет корректность координаты"""
        pass