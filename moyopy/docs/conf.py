# -- Project information -----------------------------------------------------
from moyopy import __version__ as version

project = "moyopy"
copyright = "2023, Kohei Shinohara"
repository_url = "https://github.com/spglib/moyo"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
    # "nbsphinx",
    "myst_parser",
    "autoapi.extension",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build"]

# The suffix(es) of source filenames.
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

# -----------------------------------------------------------------------------
# napoleon
# -----------------------------------------------------------------------------

# napoleon_type_aliases = {}
napoleon_use_rtype = True
napoleon_use_ivar = True

# -----------------------------------------------------------------------------
# bibtex
# -----------------------------------------------------------------------------

# https://pypi.org/project/sphinxcontrib-bibtex/
bibtex_bibfiles = ["references.bib"]

# -----------------------------------------------------------------------------
# intersphinx
# -----------------------------------------------------------------------------
intersphinx_mapping = {
    "spglib": ("https://spglib.readthedocs.io/en/stable/", None),
}

# -----------------------------------------------------------------------------
# MyST
# -----------------------------------------------------------------------------
myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "tasklist",
]
myst_dmath_double_inline = True
myst_heading_anchors = 3

# -----------------------------------------------------------------------------
# autoapi
# -----------------------------------------------------------------------------
autoapi_dirs = [
    "../python/moyopy",
]
# Drop `private-members` (default-on) so leading-underscore helpers like
# `_species_from_numbers` don't appear in the rendered API reference.
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
    "imported-members",
]
# We hand-curate the API Reference into grouped sections in `api.md` and use
# autoapi's directives there. Stop autoapi from generating its own per-module
# stub pages so we don't end up with two competing entry points.
autoapi_add_toctree_entry = False
autoapi_generate_api_docs = False

# -----------------------------------------------------------------------------
# sphinx-book-theme
# -----------------------------------------------------------------------------
html_theme = "sphinx_book_theme"
html_title = project + " " + version
html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
    "navigation_with_keys": True,
    "globaltoc_includehidden": "true",
    "show_toc_level": 2,
}

# hide sphinx footer
html_show_sphinx = False
html_show_sourcelink = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
