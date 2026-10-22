import os
from readers.fasta import FastaReader
from readers.fastq import FastqReader
from readers.sam import SamReader

def run_fasta_demo():
    fasta_path = os.path.join("A0A6I4IPK0.fasta")

    if not os.path.exists(fasta_path):
        print(f"Файла не существует")
        return

    print("Тестирование FastaReader")
    reader = FastaReader(fasta_path)

    count = reader.get_sequence_count()
    avg_len = reader.get_average_length()
    id = reader.get_sequence_ID()

    print(f"Количество последовательностей: {count}")
    print(f"Средняя длина последовательности: {avg_len:.2f}")
    print(f"Список найденных ID: {id}")

def run_fastq_demo():
    print("Тестирование FastqReader")
    fastq_path = os.path.join("test_data", "sample.fastq")
    if not os.path.exists(fastq_path):
        print(f"Файла не существует:{fastq_path}")
        return
    reader = FastqReader(fastq_path)
    count = reader.count_sequences()
    avg_len = reader.get_average_length()
    print(f"Количество ридов:{count}")
    print(f"Средняя длина:{avg_len:.1f} bp")
    print("Строим графики (matplotlib)...")
    reader.plot_quality(save_path="fastq_report.png")
    print("Успех! Графики сохранены в файл 'fastq_report.png'\n")
    

def run_sam_demo():
    print("\n--- Тестирование SamReader ---")
    sam_path = os.path.join("test_data", "test.sam")

    if not os.path.exists(sam_path):
        print(f"Файла SAM не существует по пути: {sam_path}")
        return

    reader = SamReader(sam_path)
    print(f"Заголовки PG/RG: {reader.parse_headers()}")
    print("\nСтатистика по хромосомам:")
    print(reader.get_chromosome_stats())
    print("\nПоиск в регионе chr1:100-160:")
    print(reader.get_region_intersect("chr1", 100, 160))

if __name__ == "__main__":
    run_fasta_demo()
    run_fastq_demo()
    run_sam_demo()

