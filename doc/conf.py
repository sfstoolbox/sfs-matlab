# -*- coding: utf-8 -*-

import sys
import os
import shlex
import sphinx_rtd_theme
#import subprocess

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('./_include'))

from acronyms import rst_epilog # This includes things like |HRTF| etc.
import version

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '1.3'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
#extensions = ['sphinx.ext.autodoc','nbsphinx','sphinx.ext.mathjax']
extensions = [
	'sphinx.ext.autodoc',
        'mathjax',
	'sphinx.ext.viewcode',
	'matplotlib.sphinxext.only_directives',
	'matplotlib.sphinxext.plot_directive',
]

# Enable numbering of figures and tables
numfig = True

# Plot settings ofr matplot
plot_include_source = True
plot_html_show_source_link = False
plot_html_show_formats = False
plot_formats = ['png']
plot_rcparams = {'figure.figsize' : [8, 4.5] }

# Add any paths that contain templates here, relative to this directory.
templates_path = ['./_theme/sfs/static/']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = '.txt'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'contents'

# General information about the project. (substitutions)
project = 'SFS Toolbox - Documentation'
copyright = '2016, SFS Toolbox Team'
author = 'SFS Toolbox Team'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = version.get_version()

# The full version, including alpha/beta/rc tags.
release = version

html_context = {'versions': [('1.0', '1.0/'), ('latest', 'latest/')],
                'downloads': [('PDF', '/sfs-toolbox-documentation.pdf')],
                'home_url': 'http://sfstoolbox.org/doc/'}

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = True

# The name of the Pygments (syntax highlighting) style to use.
#pygments_style = 'sphinx'
pygments_style = 'trac'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sfs"

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = ["./_theme",sphinx_rtd_theme.get_html_theme_path()]

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "SFS Toolbox"

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = ""

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = False

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Language to be used for generating the HTML full-text search index.
# Sphinx supports the following languages:
#   'da', 'de', 'en', 'es', 'fi', 'fr', 'hu', 'it', 'ja'
#   'nl', 'no', 'pt', 'ro', 'ru', 'sv', 'tr'
#html_search_language = 'en'

# A dictionary with options for the search language support, empty by default.
# Now only 'ja' uses this config value
#html_search_options = {'type': 'default'}

# The name of a javascript file (relative to the configuration directory) that
# implements a search results scorer. If empty, the default will be used.
#html_search_scorer = 'scorer.js'

# Output file base name for HTML help builder.
htmlhelp_basename = 'sfs-doc'

# Appended to every page
rst_epilog = rst_epilog + """
.. |SFS Toolbox|     replace:: SFS Tooblox
"""

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',

# Latex figure (float) alignment
#'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
  (master_doc, 'sfs-toolbox-documentation.tex', u'SFS Toolbox - Documentation',
   u'SFS Toolbox team', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


