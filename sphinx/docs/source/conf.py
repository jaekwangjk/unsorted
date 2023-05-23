# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os


# -- Project information -----------------------------------------------------

project = 'statsticsGBM '
copyright = '2023, Jaekwang Kim'
author = 'Jaekwang Kim'

# The full version, including alpha/beta/rc tags
release = '1.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinxcontrib.matlab', 'sphinx.ext.autodoc',
'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.imgmath', 
    'sphinx.ext.todo'
]


this_dir = os.path.dirname(os.path.abspath(__file__))

#matlab_src_dir = os.path.abspath(os.path.join(this_dir, '../../'))
matlab_src_dir = os.path.abspath(os.path.join(this_dir, '../../statisticsGBM/'))

templates_path = ['_templates']
exclude_patterns = []
html_theme = 'sphinx_book_theme'
html_static_path = ['_static']