
from base import GenomicDataReader
import pandas as pd


class VcfReader(GenomicDataReader):
    """Класс для чтения и анализа VCF-файлов."""

    # Чтение вариантов с помощью генератора
    def read(self):
        found_header = False

        with open(self.filename, "r", encoding="utf-8") as file:
            for line in file:

                if not line.strip() or line.startswith("##"):
                    continue

                if line.startswith("#CHROM"):
                    columns = line.strip().split("\t")

                    if columns[:8] != [
                        "#CHROM", "POS", "ID", "REF",
                        "ALT", "QUAL", "FILTER", "INFO"
                    ]:
                        raise ValueError("Неправильные столбцы VCF")

                    found_header = True
                    continue

                if line.startswith("#"):
                    continue

                if not found_header:
                    raise ValueError("Нет заголовка #CHROM")

                yield self.analyze_line(line)

        if not found_header:
            raise ValueError("Нет заголовка #CHROM")

    # Разбор одной строки VCF
    def analyze_line(self, line):
        parts = line.strip().split("\t")

        if len(parts) < 8:
            raise ValueError("В строке VCF меньше восьми полей")

        pos = int(parts[1])

        if pos < 1:
            raise ValueError("Позиция должна быть больше нуля")

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

    # Получение заголовков
    def get_headers(self):
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
                    headers["columns"] = line.strip().split("\t")
                    break

                elif line.strip() and not line.startswith("#"):
                    break

        if "columns" not in headers:
            raise ValueError("Нет заголовка #CHROM")

        return headers

    # Получение отдельной группы заголовков
    def get_header_group(self, group):
        return self.get_headers().get(group, [])

    # Подсчёт вариантов
    def count_variants(self):
        count = 0

        for variant in self.read():
            count += 1

        return count

    # Получение списка хромосом
    def get_chromosomes(self):
        chromosomes = set()

        for variant in self.read():
            chromosomes.add(variant["chrom"])

        return sorted(chromosomes)

    # Получение референсного генома
    def get_refrence_genome(self):
        reference = self.get_header_group("reference")

        if reference:
            return reference[0]

        return "Не указан"

    # Проверка координат
    def validate_coordinate(self, chrom, pos):
        if type(pos) is not int or pos < 1:
            return False

        return chrom in self.get_chromosomes()

    # Статистика вариантов по геномным регионам
    def get_statistics(self, region_size=1000000):
        if region_size < 1:
            raise ValueError("Размер региона должен быть больше нуля")

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

        table = pd.DataFrame(
            rows,
            columns=["chrom", "start", "end", "count"]
        )

        return table.sort_values(["chrom", "start"])

    # Поиск вариантов в геномном интервале
    def get_region_intersect(self, chrom, start, end):
        if start < 1 or end < start:
            raise ValueError("Неправильные координаты региона")

        for variant in self.read():

            if variant["chrom"] != chrom:
                continue

            pos = variant["pos"]
            variant_end = pos + len(variant["ref"]) - 1

            # Учитываем поле END, если оно есть
            for field in variant["info"].split(";"):
                if field.startswith("END="):
                    variant_end = int(field.split("=", 1)[1])
                    break

            if pos <= end and variant_end >= start:
                yield variant


# Запуск программы для выбранного файла
if __name__ == "__main__":

    path = input("Путь к VCF-файлу: ").strip().strip('"').strip("'")

    reader = VcfReader(path)

    try:
        print("\nЗАГОЛОВКИ:")

        for group, values in reader.get_headers().items():
            print(group + ":")

            for value in values:
                print("  ", value)

        print("\nКОЛИЧЕСТВО ВАРИАНТОВ:")
        print(reader.count_variants())

        print("\nХРОМОСОМЫ:")
        print(reader.get_chromosomes())

        print("\nРЕФЕРЕНСНЫЙ ГЕНОМ:")
        print(reader.get_refrence_genome())

        print("\nСТАТИСТИКА ПО РЕГИОНАМ:")
        print(reader.get_statistics().to_string(index=False))

        print("\nПОИСК ПО ГЕНОМНОМУ ИНТЕРВАЛУ")

        chrom = input(
            "Хромосома (Enter — пропустить): "
        ).strip()

        if chrom:
            start = int(input("Начало региона: "))
            end = int(input("Конец региона: "))

            count = 0

            for variant in reader.get_region_intersect(
                chrom, start, end
            ):
                print(variant)
                count += 1

            print("Найдено вариантов:", count)

        answer = input(
            "\nПоказать все варианты? (да/нет): "
        ).strip().lower()

        if answer == "да":
            for variant in reader.read():
                print(variant)

    except (OSError, UnicodeError, ValueError) as error:
        print("Ошибка:", error)
