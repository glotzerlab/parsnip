"""Configuration file for the Sphinx documentation builder."""

import datetime

CURRENT_YEAR = datetime.date.today().year
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "parsnip"
copyright = f"2024-{ CURRENT_YEAR } The Regents of the University of Michigan"
author = "Jen Bradley"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

add_module_names = False # Simplify paths to parsnip.oo

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "autodocsumm",
]

templates_path = ["_templates"]
exclude_patterns = ["build", "Thumbs.db", ".DS_Store"]
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

autodoc_default_options = {
    "inherited-members": True,
    "show-inheritance": True,
    "autosummary": True,
}
autodoc_typehints = "description"

pygments_style = "friendly"
pygments_dark_style = "native"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
source_repository = "https://github.com/glotzerlab/parsnip/"
html_theme = "furo"
html_static_path = ["_static"]
html_theme_options = {
    "sidebar_hide_name": True,
    "light_logo": "parsnip_header_dark.svg",
    "dark_logo": "parsnip_header_light.svg",
    "dark_css_variables": {
        "color-brand-primary": "#4AA092",
        "color-brand-content": "#5187b2",
    },
    "light_css_variables": {
        "color-brand-primary": "#005A50",
        "color-brand-content": "#406a8c",
    },
    "source_edit_link": "https://github.com/glotzerlab/parsnip/edit/main/doc/source/{filename}",
    "source_view_link": "https://github.com/glotzerlab/parsnip",
}

html_favicon = "_static/parsnip_logo_favicon.svg"
