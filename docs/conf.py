
import os
import sys

sys.path.insert(
    0, os.path.abspath("..")
)

project = "Biological Data Formats"
copyright = "2026, Biological Data Formats Team"
author = "Biological Data Formats Team"
release = "1.0.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = []

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

language = "ru"

html_theme = "alabaster"
html_static_path = []
