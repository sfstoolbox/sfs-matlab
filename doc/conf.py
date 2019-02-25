# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import sphinx_rtd_theme
# Allow local modules and extensions
sys.path.insert(0, os.path.abspath('.'))  # for acronyms.py
from acronyms import acronyms  # |HRTF| etc.


# -- GENERAL -------------------------------------------------------------

project = 'SFS Toolbox'
copyright = '2010-2019 SFS Toolbox Developers'
author = 'SFS Toolbox Developers'

needs_sphinx = '1.3'  # minimal sphinx version
extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.mathjax',
        'sphinx.ext.viewcode',
]
master_doc = 'index'
source_suffix = '.rst'
exclude_patterns = ['_build']
# The full version, including alpha/beta/rc tags.
#release = version
try:
    release = subprocess.check_output(
            ('git', 'describe', '--tags', '--always'))
    release = release.decode().strip()
except Exception:
    release = '<unknown>'

# Append acronyms at the end of every page
rst_epilog = acronyms


# -- FIGURES AND CODE ----------------------------------------------------
numfig = True
# Code syntax highlighting style
pygments_style = 'trac'


# -- HTML ----------------------------------------------------------------

def setup(app):
    """Include custom theme files to sphinx HTML header"""
    app.add_stylesheet('css/abbr.css')
    app.add_stylesheet('css/title.css')

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']
templates_path = ['_template']
html_title = "SFS Toolbox"
html_short_title = ""
htmlhelp_basename = 'sfs-matlab'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    'preamble': '',
    'figure_align': 'htbp',
}
# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
  (master_doc, 'sfs-matlab.tex', u'SFS Toolbox - Matlab Documentation',
   u'SFS Toolbox Developers', 'manual'),
]
# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None
