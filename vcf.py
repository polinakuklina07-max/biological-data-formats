
from base import GenomicDataReader
import pandas as pd

class VcfReader(GenomicDataReader):
    """Класс для чтения и анализа VCF-файлов."""

    def get_headers(self):
        """Возвращает заголовки VCF."""
        headers = {}
        with open(self.filename, "r", encoding="utf-8") as file:
            for line in file:
                if line.startswith("##"):
                    key, sep, value = line[2:].strip().partition("=")
                    if sep:
                        if key not in headers:
                            headers[key] = []
                        headers[key].append(value)
                elif line.startswith("#CHROM"):
                    columns = line.strip().split("\t")
                    if columns[:8] != [
                        "#CHROM", "POS", "ID", "REF",
                        "ALT", "QUAL", "FILTER", "INFO"
                    ]:
                        raise ValueError("Неправильные столбцы VCF")
                    if len(columns) > 8 and (columns[8] != "FORMAT" or len(columns) < 10):
                        raise ValueError("Неправильные столбцы с образцами")
                    headers["columns"] = columns
                    break
                elif line.strip() and not line.startswith("#"):
                    raise ValueError("Данные появились до заголовка #CHROM")
        if "columns" not in headers:
            raise ValueError("Нет заголовка #CHROM")
        return headers

    def get_header_group(self, group):
        """Возвращает указанную группу заголовков."""
        return self.get_headers().get(group, [])

    def read(self):
        """Читает варианты по одному с помощью генератора."""
        columns = self.get_headers()["columns"]
        with open(self.filename, "r", encoding="utf-8") as file:
            for line in file:
                if not line.strip() or line.startswith("#"):
                    continue
                if len(line.strip().split("\t")) != len(columns):
                    raise ValueError("Количество полей не совпадает с заголовком")
                yield self.analyze_line(line)

    def analyze_line(self, line):
        """Превращает строку VCF в словарь."""
        parts = line.strip().split("\t")
        if len(parts) < 8:
            raise ValueError("В строке VCF меньше восьми полей")
        try:
            pos = int(parts[1])
        except ValueError:
            raise ValueError("Позиция должна быть целым числом")
        if pos < 1:
            raise ValueError("Позиция должна быть больше нуля")
        if not parts[0] or parts[0] == ".":
            raise ValueError("Не указана хромосома")
        if not parts[3] or parts[3] == ".":
            raise ValueError("Не указан REF")
        if not parts[4]:
            raise ValueError("Не указан ALT")
        return {
            "chrom": parts[0],
            "pos": pos,
            "id": parts[2],
            "ref": parts[3],
            "alt": parts[4],
            "qual": parts[5],
            "filter": parts[6],
            "info": parts[7],
            "format": parts[8] if len(parts) > 8 else "",
            "samples": parts[9:]
        }

    def count_variants(self):
        """Считает количество вариантов."""
        count = 0
        for variant in self.read():
            count += 1
        return count

    def get_chromosomes(self):
        """Возвращает список хромосом с вариантами."""
        chromosomes = set()
        for variant in self.read():
            chromosomes.add(variant["chrom"])
        return sorted(chromosomes)

    def get_refrence_genome(self):
        """Возвращает название референсного генома."""
        reference = self.get_header_group("reference")
        if reference:
            return reference[0]
        return "Не указан"

    def validate_coordinate(self, chrom, pos):
        """Проверяет хромосому и положительную позицию."""
        if type(pos) is not int or pos < 1:
            return False
        return chrom in self.get_chromosomes()

    def get_statistics(self, region_size=1000000):
        """Считает варианты по регионам и возвращает таблицу pandas."""
        if type(region_size) is not int or region_size < 1:
            raise ValueError("Размер региона должен быть положительным целым числом")
        counts = {}
        for variant in self.read():
            chrom = variant["chrom"]
            pos = variant["pos"]
            start = ((pos - 1) // region_size) * region_size + 1
            region = (chrom, start)
            if region not in counts:
                counts[region] = 0
            counts[region] += 1
        rows = []
        for (chrom, start), count in counts.items():
            rows.append({
                "chrom": chrom,
                "start": start,
                "end": start + region_size - 1,
                "count": count
            })
        table = pd.DataFrame(rows, columns=["chrom", "start", "end", "count"])
        return table.sort_values(["chrom", "start"])

    def get_region_intersect(self, chrom, start, end):
        """Выдаёт варианты, пересекающие указанный интервал."""
        if type(start) is not int or type(end) is not int or start < 1 or end < start:
            raise ValueError("Неправильные координаты региона")
        for variant in self.read():
            if variant["chrom"] != chrom:
                continue
            variant_start = variant["pos"]
            variant_end = variant_start + len(variant["ref"]) - 1
            for field in variant["info"].split(";"):
                if field.startswith("END="):
                    variant_end = int(field.split("=", 1)[1])
                    if variant_end < variant_start:
                        raise ValueError("END меньше позиции варианта")
                    break
            if variant_start <= end and variant_end >= start:
                yield variant

if __name__ == "__main__":
    path = input("Путь к VCF-файлу: ").strip().strip('"').strip("'")
    reader = VcfReader(path)
    try:
        print("\nЗАГОЛОВКИ:")
        for group, values in reader.get_headers().items():
            print(group + ":", values)
        print("\nКОЛИЧЕСТВО ВАРИАНТОВ:")
        print(reader.count_variants())
        print("\nХРОМОСОМЫ:")
        print(reader.get_chromosomes())
        print("\nРЕФЕРЕНСНЫЙ ГЕНОМ:")
        print(reader.get_refrence_genome())
        print("\nСТАТИСТИКА ПО РЕГИОНАМ:")
        print(reader.get_statistics().to_string(index=False))
        print("\nПОИСК ПО ГЕНОМНОМУ ИНТЕРВАЛУ")
        chrom = input("Хромосома (Enter — пропустить): ").strip()
        if chrom:
            start = int(input("Начало региона: "))
            end = int(input("Конец региона: "))
            count = 0
            for variant in reader.get_region_intersect(chrom, start, end):
                print(variant)
                count += 1
            print("Найдено вариантов:", count)
        answer = input("\nПоказать все варианты? (да/нет): ").strip().lower()
        if answer == "да":
            for variant in reader.read():
                print(variant)
    except (OSError, UnicodeError, ValueError) as error:
        print("Ошибка:", error)
