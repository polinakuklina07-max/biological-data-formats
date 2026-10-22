"""Чтение и анализ файлов VCF (Variant Call Format)."""
from pathlib import Path
from collections import defaultdict
from typing import Iterator
import re
import pandas as pd
from base import GenomicDataReader
class VcfReader(GenomicDataReader):
    """Reader для текстовых файлов VCF."""
    REQUIRED_COLUMNS = [
        "#CHROM", "POS", "ID", "REF",
        "ALT", "QUAL", "FILTER", "INFO"
    ]
    def __init__(self, filename: str):
        super().__init__(filename)
        self.filename = str(filename)
    def get_headers(self) -> dict:
        """Возвращает метаданные и названия столбцов VCF."""
        headers = {}
        found_columns = False
        with open(self.filename, "r", encoding="utf-8") as file:
            for line in file:
                line = line.rstrip("\r\n")
                if line.startswith("##"):
                    key, separator, value = line[2:].partition("=")
                    if separator:
                        headers.setdefault(key, []).append(value)
                elif line.startswith("#CHROM"):
                    columns = line.split("\t")
                    self._validate_columns(columns)
                    headers["columns"] = columns
                    found_columns = True
                    break
                elif line.strip() and not line.startswith("#"):
                    raise ValueError(
                        "VCF: записи обнаружены до заголовка #CHROM"
                    )
        if not found_columns:
            raise ValueError("VCF: отсутствует заголовок #CHROM")
        return headers
    @classmethod
    def _validate_columns(cls, columns: list) -> None:
        """Проверяет структуру заголовка."""
        if columns[:8] != cls.REQUIRED_COLUMNS:
            raise ValueError("VCF: неверные обязательные столбцы")
        if len(columns) > 8:
            if columns[8] != "FORMAT":
                raise ValueError("VCF: девятый столбец должен быть FORMAT")
            if len(columns) == 9:
                raise ValueError("VCF: FORMAT указан без образцов")
        if len(columns) != len(set(columns)):
            raise ValueError("VCF: повторяющиеся названия столбцов")
    def get_header_group(self, group: str) -> list:
        """Возвращает элементы определённой группы заголовков."""
        return self.get_headers().get(group, [])
    def get_refrence_genome(self) -> str:
        """Референсный геном. Имя сохранено для совместимости с base.py."""
        values = self.get_header_group("reference")
        return values[0] if values else "Не указан"
    def get_reference_genome(self) -> str:
        """Альтернативное имя метода с правильным написанием."""
        return self.get_refrence_genome()
    def _get_contig_lengths(self) -> dict:
        """Получает объявленные длины хромосом из ##contig."""
        lengths = {}
        for record in self.get_header_group("contig"):
            match_id = re.search(r"(?:^|,)ID=([^,>]+)", record)
            match_len = re.search(r"(?:^|,)length=(\d+)", record)
            if match_id and match_len:
                lengths[match_id.group(1)] = int(match_len.group(1))
        return lengths
    def analyze_line(self, line: str) -> dict:
        """Проверяет строку и возвращает одну запись VCF."""
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) < 8:
            raise ValueError("VCF: запись содержит меньше 8 полей")
        chrom, raw_pos, identifier, ref, alt, qual, flt, info = parts[:8]
        if not chrom or chrom == ".":
            raise ValueError("VCF: отсутствует имя хромосомы")
        try:
            pos = int(raw_pos)
        except ValueError as error:
            raise ValueError("VCF: POS должен быть целым числом") from error
        if pos < 1:
            raise ValueError("VCF: POS должен быть положительным")
        if not re.fullmatch(r"[ACGTNacgtn]+", ref):
            raise ValueError("VCF: некорректное значение REF")
        if not alt:
            raise ValueError("VCF: отсутствует ALT")
        for allele in alt.split(","):
            valid_sequence = re.fullmatch(r"[ACGTNacgtn]+", allele)
            valid_symbolic = re.fullmatch(r"<[A-Za-z0-9_:.-]+>", allele)
            valid_breakend = "[" in allele or "]" in allele
            if not (
                valid_sequence
                or valid_symbolic
                or valid_breakend
                or allele in {".", "*"}
            ):
                raise ValueError(f"VCF: некорректное ALT: {allele}")
        if qual != ".":
            try:
                quality = float(qual)
            except ValueError as error:
                raise ValueError("VCF: некорректное QUAL") from error
            if not (0 <= quality < float("inf")):
                raise ValueError("VCF: QUAL должен быть конечным и >= 0")
        if len(parts) > 8:
            if len(parts) < 10:
                raise ValueError("VCF: FORMAT указан без образцов")
            if not parts[8] or parts[8] == ".":
                raise ValueError("VCF: некорректное FORMAT")
        return {
            "chrom": chrom,
            "pos": pos,
            "id": identifier,
            "ref": ref,
            "alt": alt,
            "qual": qual,
            "filter": flt,
            "info": info,
            "format": parts[8] if len(parts) > 8 else "",
            "samples": parts[9:]
        }
    def read(self) -> Iterator[dict]:
        """Последовательно выдаёт варианты, используя генератор."""
        headers = self.get_headers()
        expected_columns = len(headers["columns"])
        with open(self.filename, "r", encoding="utf-8") as file:
            for line_number, line in enumerate(file, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                try:
                    if len(line.rstrip("\r\n").split("\t")) != expected_columns:
                        raise ValueError(
                            "число полей не совпадает с заголовком"
                        )
                    yield self.analyze_line(line)
                except ValueError as error:
                    raise ValueError(
                        f"VCF: ошибка в строке {line_number}: {error}"
                    ) from error
    def count_variants(self) -> int:
        """Подсчитывает число вариантов."""
        return sum(1 for _ in self.read())
    def get_chromosomes(self) -> list:
        """Возвращает отсортированный список хромосом с вариантами."""
        return sorted({variant["chrom"] for variant in self.read()})
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """Проверяет координату по хромосоме и её известной длине."""
        if not isinstance(chrom, str):
            return False
        if type(pos) is not int or pos < 1:
            return False
        lengths = self._get_contig_lengths()
        if chrom in lengths:
            return pos <= lengths[chrom]
        return chrom in self.get_chromosomes()
    @staticmethod
    def _variant_interval(variant: dict) -> tuple:
        """Возвращает покрываемый вариантом интервал 1-based."""
        start = variant["pos"]
        end = start + len(variant["ref"]) - 1
        for entry in variant["info"].split(";"):
            if entry.startswith("END="):
                try:
                    end = int(entry.split("=", 1)[1])
                except ValueError as error:
                    raise ValueError("VCF: некорректное INFO/END") from error
                if end < start:
                    raise ValueError("VCF: END меньше POS")
                break
        return start, end
    def get_statistics(self, region_size: int = 1_000_000) -> pd.DataFrame:
        """Статистика числа вариантов в непересекающихся регионах."""
        if type(region_size) is not int or region_size < 1:
            raise ValueError("Размер региона должен быть положительным целым")
        counts = defaultdict(int)
        for variant in self.read():
            chrom = variant["chrom"]
            pos = variant["pos"]
            start = ((pos - 1) // region_size) * region_size + 1
            counts[(chrom, start)] += 1
        rows = [
            {
                "chrom": chrom,
                "start": start,
                "end": start + region_size - 1,
                "count": count
            }
            for (chrom, start), count in counts.items()
        ]
        columns = ["chrom", "start", "end", "count"]
        result = pd.DataFrame(rows, columns=columns)
        if result.empty:
            return result
        return result.sort_values(
            ["chrom", "start"]
        ).reset_index(drop=True)
    def get_region_intersect(
        self, chrom: str, start: int, end: int
    ) -> Iterator[dict]:
        """Выдаёт варианты, пересекающие интервал [start, end]."""
        if not isinstance(chrom, str) or not chrom:
            raise ValueError("Необходимо указать хромосому")
        if type(start) is not int or type(end) is not int:
            raise ValueError("Координаты должны быть целыми числами")
        if start < 1 or end < start:
            raise ValueError("Некорректные координаты региона")
        for variant in self.read():
            if variant["chrom"] != chrom:
                continue
            variant_start, variant_end = self._variant_interval(variant)
            if variant_start <= end and variant_end >= start:
                yield variant
def main():
    """Демонстрационная программа для VCF."""
    path = input("Путь к VCF-файлу: ").strip().strip('"').strip("'")
    if not path:
        print("Ошибка: путь к файлу не указан")
        return
    if not Path(path).is_file():
        print("Ошибка: файл не найден")
        return
    reader = VcfReader(path)
    try:
        print("\n=== ЗАГОЛОВКИ ===")
        for group, values in reader.get_headers().items():
            print(f"{group}:")
            for value in values:
                print(" ", value)
        print("\n=== КОЛИЧЕСТВО ВАРИАНТОВ ===")
        print(reader.count_variants())
        print("\n=== ХРОМОСОМЫ ===")
        print(reader.get_chromosomes())
        print("\n=== РЕФЕРЕНСНЫЙ ГЕНОМ ===")
        print(reader.get_refrence_genome())
        print("\n=== СТАТИСТИКА ПО РЕГИОНАМ ===")
        statistics = reader.get_statistics()
        if statistics.empty:
            print("Варианты отсутствуют")
        else:
            print(statistics.to_string(index=False))
        print("\n=== ПОИСК ПО ИНТЕРВАЛУ ===")
        chrom = input("Хромосома (Enter — пропустить): ").strip()
        if chrom:
            start = int(input("Начало: "))
            end = int(input("Конец: "))
            matches = list(reader.get_region_intersect(chrom, start, end))
            for variant in matches:
                print(variant)
            print(f"Найдено вариантов: {len(matches)}")
        answer = input("\nПоказать все варианты? (да/нет): ").lower().strip()
        if answer == "да":
            for variant in reader.read():
                print(variant)
    except (OSError, UnicodeError, ValueError) as error:
        print(f"Ошибка: {error}")
if __name__ == "__main__":
    main()
