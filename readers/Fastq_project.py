import matplotlib.pyplot as plt
import seaborn as sns
from typing import Iterator, Dict, List

from base import SequenceReader

class FastqReader(SequenceReader):
    """
    Парсер и анализатор файлов формата FASTQ.
    Наследуется от SequenceReader.
    Использует генераторы для эффективного чтения больших файлов.
    """

    def read(self) -> Iterator[dict]:
        """
        Генератор для построчного чтения FASTQ файла.
        
        Yields:
            dict: Словарь с ключами 'header', 'sequence' и 'quality'.
        """
        with open(self.filename, 'r') as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                sequence = f.readline().strip()
                f.readline()  # Пропускаем строку "+"
                quality = f.readline().strip()
                
                yield {
                    "header": header,
                    "sequence": sequence,
                    "quality": quality
                }

    def analyze_line(self, line: str):
        """
        Реализация абстрактного метода из базового класса.
        Так как мы читаем файл целиком через генератор read(), 
        этот метод нам не нужен, но мы обязаны его создать.
        """
        pass

    def get_sequence_ID(self) -> List[str]:
        """
        Реализация абстрактного метода из базового класса.
        Возвращает список всех ID (заголовков) последовательностей.
        """
        return [record["header"] for record in self.read()]

    def validate_sequence(self, sequence: str) -> bool:
        """
        Переопределяем метод из базового класса.
        Проверяет, что последовательность состоит только из символов A, T, G, C, N.
        """
        valid_chars = set("ATGCatgcNn")
        return all(char in valid_chars for char in sequence)

    def count_sequences(self) -> int:
        """
        Возвращает общее количество последовательностей в файле.
        """
        return sum(1 for _ in self.read())

    def get_average_length(self) -> float:
        """
        Вычисляет среднюю длину последовательностей.
        """
        total_len, count = 0, 0
        for record in self.read():
            total_len += len(record["sequence"])
            count += 1
        return total_len / count if count else 0.0

    def plot_quality(self, save_path: str = None):
        """
        Строит графики качества (аналог FastQC).
        
        Args:
            save_path (str, optional): Путь для сохранения картинки. 
                                       Если None, график выводится на экран.
        """
        pos_qualities: Dict[int, List[int]] = {}
        pos_content = {}
        lengths = []

        for record in self.read():
            seq = record["sequence"]
            qual = record["quality"]
            if len(seq) != len(qual):
                print(
                    f"Пропускаем некорректную запись {record['header']}: "
                    f"длина последовательности {len(seq)}, "
                    f"длина качества {len(qual)}"
                )
                continue
            lengths.append(len(seq))
            # Перевод ASCII в Phred score (Phred+33)
            for pos in range(len(seq)):
                pos_qualities.setdefault(pos, []).append(ord(qual[pos]) - 33)
                if pos not in pos_content:
                    pos_content[pos] = {'A':0, 'T':0, 'G':0, 'C':0}
                base = seq[pos].upper()
                if base in pos_content[pos]:
                    pos_content[pos][base] += 1

        if not lengths:
            print("Нет данных для построения графиков.")
            return

        # Создаем фигуру с тремя графиками
        fig, axes = plt.subplots(1, 3, figsize=(20, 5))

        # 1. Распределение длин
        sns.histplot(lengths, kde=True, ax=axes[0], color="skyblue")
        axes[0].set_title("Sequence Length Distribution")
        axes[0].set_xlabel("Length (bp)")
        axes[0].set_ylabel("Count")

        # 2. Качество по позициям
        positions = sorted(pos_qualities.keys())
        step = max(1, len(positions) // 40)
        sampled_pos = positions[::step]
        data_to_plot = [pos_qualities[p] for p in sampled_pos]

        sns.boxplot(data=data_to_plot, ax=axes[1], color="lightgreen")
        axes[1].set_title("Per Base Sequence Quality")
        axes[1].set_xlabel("Position in read (bp)")
        axes[1].set_ylabel("Phred Score")
        axes[1].set_xticks(range(len(sampled_pos)))
        axes[1].set_xticklabels(sampled_pos, rotation=90, fontsize=8)

        # 3. Состав нуклеотидов
        a_vals, t_vals, g_vals, c_vals = [], [], [], []
        for p in positions:
            total = sum(pos_content[p].values())

            a_vals.append(pos_content[p]['A']/total*100)
            t_vals.append(pos_content[p]['T']/total*100)
            g_vals.append(pos_content[p]['G']/total*100)
            c_vals.append(pos_content[p]['C']/total*100)
            
        axes[2].plot(positions, a_vals, label='A', color='green')
        axes[2].plot(positions, t_vals, label='T', color='red')
        axes[2].plot(positions, g_vals, label='G', color='blue')
        axes[2].plot(positions, c_vals, label='C', color='orange')
                          
        axes[2].set_title("Per Base Sequence Content")
        axes[2].set_xlabel("Position in read (bp)")
        axes[2].set_ylabel("Percentage(%)")
        axes[2].legend()

        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300)
            print(f"Графики сохранены в {save_path}")
        else:
            plt.show()
