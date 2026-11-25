import gzip
import logging
from typing import Iterator, Dict, List, Optional

import matplotlib.pyplot as plt
import seaborn as sns

from base import SequenceReader

logger = logging.getLogger(__name__)


class FastqReader(SequenceReader):
    """
    Парсер и анализатор файлов формата FASTQ.
    Поддерживает .fastq / .fq, а также сжатые .fastq.gz / .fq.gz.
    """


    def _open(self):
        """Открывает файл: gzip, если расширение .gz, иначе обычный текст."""
        if self.filename.lower().endswith(".gz"):
            return gzip.open(self.filename, "rt", encoding="utf-8", errors="replace")
        return open(self.filename, "r", encoding="utf-8", errors="replace")

    def read(self) -> Iterator[dict]:
        """
        Генератор для построчного чтения FASTQ файла.
        Yields dict с ключами 'header', 'sequence', 'quality'.
        """
        with self._open() as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                sequence = f.readline().strip()
                plus = f.readline()          # строка "+"
                quality = f.readline().strip()

                # неполная запись в конце файла — прекращаем
                if not quality and not sequence:
                    break

                yield {
                    "header": header,
                    "sequence": sequence,
                    "quality": quality,
                }

    def analyze_line(self, line: str):
        """Абстрактный метод из базового класса — не используется."""
        pass


    def get_sequence_ID(self) -> List[str]:
        return [record["header"] for record in self.read()]

    def validate_sequence(self, sequence: str) -> bool:
        valid_chars = set("ATGCatgcNn")
        return all(char in valid_chars for char in sequence)

    def count_sequences(self) -> int:
        return sum(1 for _ in self.read())

    def get_average_length(self) -> float:
        total_len, count = 0, 0
        for record in self.read():
            total_len += len(record["sequence"])
            count += 1
        return total_len / count if count else 0.0


    def compute_stats(self) -> dict:
        """
        Считает всю статистику за один проход по файлу.
        Возвращает dict:
          count, avg_length, min_length, max_length,
          lengths, pos_qualities, pos_content,
          avg_quality, gc_percent, warnings
        """
        lengths: List[int] = []
        pos_qualities: Dict[int, List[int]] = {}
        pos_content: Dict[int, Dict[str, int]] = {}
        warnings: List[str] = []
        total_qual = 0
        total_bases = 0
        gc_bases = 0

        for record in self.read():
            seq = record["sequence"]
            qual = record["quality"]

            if len(seq) != len(qual):
                warnings.append(
                    f"{record['header']}: длина seq={len(seq)} != длина qual={len(qual)} — запись пропущена"
                )
                continue

            lengths.append(len(seq))

            for pos, (base, qchar) in enumerate(zip(seq, qual)):
                phred = ord(qchar) - 33
                if phred < 0 or phred > 93:
                    warnings.append(
                        f"{record['header']}: символ {qchar!r} (pos {pos}) вне диапазона Phred+33 — пропущен"
                    )
                    continue

                pos_qualities.setdefault(pos, []).append(phred)
                total_qual += phred
                total_bases += 1

                b = base.upper()
                slot = pos_content.setdefault(pos, {"A": 0, "T": 0, "G": 0, "C": 0})
                if b in slot:
                    slot[b] += 1
                    if b in ("G", "C"):
                        gc_bases += 1

        avg_length = sum(lengths) / len(lengths) if lengths else 0.0
        min_length = min(lengths) if lengths else 0
        max_length = max(lengths) if lengths else 0
        avg_quality = (total_qual / total_bases) if total_bases else 0.0
        gc_percent = (gc_bases / total_bases * 100) if total_bases else 0.0

        return {
            "count": len(lengths),
            "avg_length": avg_length,
            "min_length": min_length,
            "max_length": max_length,
            "lengths": lengths,
            "pos_qualities": pos_qualities,
            "pos_content": pos_content,
            "avg_quality": avg_quality,
            "gc_percent": gc_percent,
            "warnings": warnings,
        }


    def build_figure(self, stats: Optional[dict] = None):
        """
        Строит matplotlib.Figure с тремя графиками.
        Возвращает fig (или None, если данных нет).
        """
        if stats is None:
            stats = self.compute_stats()

        lengths = stats["lengths"]
        if not lengths:
            logger.warning("Нет данных для построения графиков.")
            return None

        pos_qualities = stats["pos_qualities"]
        pos_content = stats["pos_content"]

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
        data_to_plot = [pos_qualities[p] for p in sampled_pos if pos_qualities[p]]

        sns.boxplot(data=data_to_plot, ax=axes[1], color="lightgreen")
        axes[1].set_title("Per Base Sequence Quality")
        axes[1].set_xlabel("Position in read (bp)")
        axes[1].set_ylabel("Phred Score")
        axes[1].set_xticks(range(len(sampled_pos)))
        axes[1].set_xticklabels(sampled_pos, rotation=90, fontsize=8)

        # 3. Состав нуклеотидов
        a_vals, t_vals, g_vals, c_vals = [], [], [], []
        for p in positions:
            total = sum(pos_content[p].values()) or 1
            a_vals.append(pos_content[p]["A"] / total * 100)
            t_vals.append(pos_content[p]["T"] / total * 100)
            g_vals.append(pos_content[p]["G"] / total * 100)
            c_vals.append(pos_content[p]["C"] / total * 100)

        axes[2].plot(positions, a_vals, label="A", color="green")
        axes[2].plot(positions, t_vals, label="T", color="red")
        axes[2].plot(positions, g_vals, label="G", color="blue")
        axes[2].plot(positions, c_vals, label="C", color="orange")
        axes[2].set_title("Per Base Sequence Content")
        axes[2].set_xlabel("Position in read (bp)")
        axes[2].set_ylabel("Percentage(%)")
        axes[2].legend()

        fig.tight_layout()
        return fig


    def plot_quality(self, save_path: str = None):

        fig = self.build_figure()
        if fig is None:
            return None

        if save_path:
            fig.savefig(save_path, dpi=300)
            logger.info("Графики сохранены в %s", save_path)
            plt.close(fig)
            return None

        return fig