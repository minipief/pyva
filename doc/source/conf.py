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
import sys
#sys.path.insert(0, os.path.abspath('C:/Users/alexander/Documents/python_VA_lib/VAScript'))

# -- Set this value to your path of the pyva package to allow access to python code ------
#sys.path.insert(0, os.path.abspath('C:/Users/alexander/Documents/git_python'))
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------

project = 'pyva - python toolbox for vibroacoustics'
copyright = '2025, Dr. Alexander Peiffer'
author = 'Dr. Alexander Peiffer'

# The full version, including alpha/beta/rc tags
release = '1.3.0'

html_baseurl = 'https://pyva.eu'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.imgmath',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx_sitemap'
]

html_baseurl = 'https://pyva.eu/'

# Allow for autosummery files
autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_theme = 'sphinx_rtd_theme'
html_logo = 'pyva_logo.png'
html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'style_nav_header_background': 0x15b9b4 # 29808b
#    'stickysidebar': True 
}

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#
html_favicon = 'favicon.ico'


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']