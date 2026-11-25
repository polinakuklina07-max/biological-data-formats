"""
Tkinter GUI для анализа FASTQ-файлов.
"""

import os
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import tkinter.font as tkfont

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk
)
from matplotlib.figure import Figure

try:
    from tkinterdnd2 import TkinterDnD, DND_FILES
    HAS_DND = True
except ImportError:
    HAS_DND = False

from readers.Fastq_project import FastqReader

FASTQ_FILETYPES = [
    ("FASTQ files", "*.fastq *.fq *.fastq.gz *.fq.gz"),
    ("All files", "*.*"),
]
FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


class FastqApp:
    """Главное окно приложения."""

    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Урезанный FastQC")
        self.root.geometry("1200x700")
        self.root.minsize(900, 600)

        self.default_font = tkfont.Font(family="Segoe UI", size=12) 
        self.bold_font = tkfont.Font(family="Segoe UI", size=12, weight="bold")
        self.mono_font = tkfont.Font(family="Consolas", size=11)    

        self.current_path: str | None = None

        # Применяе шрифт ко всем стандартным элементам
        self.root.option_add("*Font", self.default_font)
        self.root.option_add("*TButton.Font", self.default_font)
        self.root.option_add("*TCheckbutton.Font", self.default_font)
        self.root.option_add("*TNotebook.Tab.Font", self.bold_font)

        self.current_path: str | None = None
        self._canvas = None       # FigureCanvasTkAgg
        self._toolbar = None      # NavigationToolbar2Tk
        self._fig = None          # текущая Figure

        self._build_ui()
        self._setup_dnd()


    def _build_ui(self):
        # Панель файла 
        top = ttk.Frame(self.root, padding=(10, 8))
        top.pack(side=tk.TOP, fill=tk.X)

        self.btn_open = ttk.Button(top, text="Открыть FASTQ",
                                   command=self.open_dialog)
        self.btn_open.grid(row=0, column=0, padx=(0, 8), sticky="w")

        self.path_var = tk.StringVar(value="(файл не выбран)")
        
        entry = ttk.Entry(top, textvariable=self.path_var, state="readonly", 
                          font=self.default_font)
        entry.grid(row=0, column=1, sticky="ew")

        self.btn_run = ttk.Button(top, text="Анализировать",
                                  command=self.run_analysis, state="disabled")
        self.btn_run.grid(row=1, column=0, padx=(0, 8), pady=(6, 0), sticky="w")

        self.show_warnings = tk.BooleanVar(value=False)
        ttk.Checkbutton(top, text="Показывать предупреждения",
                        variable=self.show_warnings,
                        command=self._refresh_summary)\
            .grid(row=1, column=1, pady=(6, 0), sticky="w")

        top.columnconfigure(1, weight=1)

        # Вкладки
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=(0, 6))

        # "Сводка"
        self.tab_summary = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_summary, text="Сводка")

        self.summary_text = tk.Text(self.tab_summary, wrap="word",
                                    font=self.mono_font, state="disabled")

        
        scroll = ttk.Scrollbar(self.tab_summary, orient="vertical",
                               command=self.summary_text.yview)
        self.summary_text.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.summary_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # "Графики"
        self.tab_plots = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_plots, text="Графики")

        self.plot_placeholder = ttk.Label(
            self.tab_plots,
            text="Выберите FASTQ-файл и нажмите «Анализировать»",
            anchor="center",
            font=self.bold_font
        )
        self.plot_placeholder.pack(fill=tk.BOTH, expand=True)

        # Статус-бар
        self.status_var = tk.StringVar(value="Готово")
        
        status = ttk.Label(self.root, textvariable=self.status_var,
                           relief="sunken", anchor="w", padding=(6, 2), 
                           font=self.default_font)
        status.pack(side=tk.BOTTOM, fill=tk.X)

        self._last_stats = None

    def _setup_dnd(self):
        if not HAS_DND:
            return
        try:
            self.root.drop_target_register(DND_FILES)
            self.root.dnd_bind("<<Drop>>", self._on_drop)
        except Exception as e:
            print("DND недоступен:", e)

    def _on_drop(self, event):
        raw = event.data.strip()
        if raw.startswith("{") and raw.endswith("}"):
            raw = raw[1:-1]
        path = raw
        if os.path.isfile(path) and path.lower().endswith(FASTQ_EXTS):
            self._set_path(path)
            self.run_analysis()
        else:
            messagebox.showwarning("Неверный файл",
                                   "Перетащите файл .fastq / .fq / .fastq.gz / .fq.gz")

# Actions
    def open_dialog(self):
        paths = filedialog.askopenfilenames(
            title="Выберите FASTQ-файл",
            filetypes=FASTQ_FILETYPES,
        )
        if not paths:
            return
        if len(paths) > 1:
            messagebox.showinfo(
                "Выбрано несколько файлов",
                "Анализируется только первый выбранный файл.",
            )
        self._set_path(paths[0])

    def _set_path(self, path: str):
        if not path.lower().endswith(FASTQ_EXTS):
            messagebox.showerror("Неверный формат",
                                 "Поддерживаются только .fastq")
            return
        self.current_path = path
        self.path_var.set(path)
        self.btn_run.configure(state="normal")
        self.status_var.set(f"Выбран файл: {os.path.basename(path)}")

    def run_analysis(self):
        if not self.current_path:
            return
        self.btn_open.configure(state="disabled")
        self.btn_run.configure(state="disabled")
        self.status_var.set("Читаю файл...")

        # анализ в отдельном потоке, чтобы окно не висело
        threading.Thread(target=self._analyze_worker,
                         args=(self.current_path,), daemon=True).start()

    def _analyze_worker(self, path: str):
        try:
            reader = FastqReader(path)
            stats = reader.compute_stats()
            fig = reader.build_figure(stats) if stats["count"] else None
            self.root.after(0, self._on_analyzed, path, stats, fig)
        except Exception as e:
            self.root.after(0, self._on_error, str(e))

    def _on_analyzed(self, path, stats, fig):
        self._last_stats = stats
        self._render_summary(path, stats)
        self._render_figure(fig)
        self.status_var.set(
            f"Готово: {stats['count']} ридов, средняя длина {stats['avg_length']:.1f} bp"
        )
        self.btn_open.configure(state="normal")
        self.btn_run.configure(state="normal")

    def _on_error(self, msg: str):
        self.status_var.set(f"Ошибка: {msg}")
        messagebox.showerror("Ошибка анализа", msg)
        self.btn_open.configure(state="normal")
        self.btn_run.configure(state="normal")


    def _render_summary(self, path, stats):
        self.summary_text.configure(state="normal")
        self.summary_text.delete("1.0", tk.END)

        lines = [
            f"Файл:              {path}",
            f"Размер файла:      {os.path.getsize(path) / 1024:.1f} КБ",
            f"Число ридов:       {stats['count']}",
            f"Средняя длина:     {stats['avg_length']:.1f} bp",
            f"Мин. длина:        {stats['min_length']} bp",
            f"Макс. длина:       {stats['max_length']} bp",
            f"Среднее качество:  {stats['avg_quality']:.1f} (Phred)",
            f"GC-состав:         {stats['gc_percent']:.1f} %",
        ]
        self.summary_text.insert(tk.END, "\n".join(lines))

        if self.show_warnings.get() and stats["warnings"]:
            self.summary_text.insert(tk.END,
                "\n\n⚠ Предупреждения:\n" +
                "\n".join(f"  - {w}" for w in stats["warnings"][:50])
            )
            if len(stats["warnings"]) > 50:
                self.summary_text.insert(
                    tk.END,
                    f"\n  ... и ещё {len(stats['warnings']) - 50} шт."
                )

        self.summary_text.configure(state="disabled")
        self.notebook.select(self.tab_summary)

    def _refresh_summary(self):
        # перерисовать вкладку Сводка, если пользователь переключил чекбокс
        if self._last_stats and self.current_path:
            self._render_summary(self.current_path, self._last_stats)

    def _render_figure(self, fig):
        # очистить место заглушки
        self.plot_placeholder.pack_forget()

        if self._canvas is not None:
            self._canvas.get_tk_widget().destroy()
            self._canvas = None
        if self._toolbar is not None:
            self._toolbar.destroy()
            self._toolbar = None

        if fig is None:
            self.plot_placeholder.pack(fill=tk.BOTH, expand=True)
            return

        self._fig = fig
        self._canvas = FigureCanvasTkAgg(fig, master=self.tab_plots)
        self._canvas.draw()
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self._toolbar = NavigationToolbar2Tk(self._canvas, self.tab_plots)
        self._toolbar.update()
        self.notebook.select(self.tab_plots)


def launch_gui():
    """Запуск GUI."""
    if HAS_DND:
        root = TkinterDnD.Tk()
    else:
        root = tk.Tk()
    FastqApp(root)
    root.mainloop()