import os
from fasta import FastaReader
def run_fasta_demo():
    fasta_path = os.path.join("A0A6I4IPK0.fasta")

    if not os.path.exists(fasta_path):
        print(f"Файла не существует")
        return

    print(f"Содержимое файла: {fasta_path}")
    with open(fasta_path, 'r') as f:
        print(f.read())

    print()

    print("Тестирование FastaReader")
    reader = FastaReader(fasta_path)

    count = reader.get_sequence_count()
    avg_len = reader.get_average_length()
    id = reader.get_sequence_ID()

    print(f"Количество последовательностей: {count}")
    print(f"Средняя длина последовательности: {avg_len:.2f}")
    print(f"Список найденных ID: {id}")

if __name__ == "__main__":
    run_fasta_demo()

