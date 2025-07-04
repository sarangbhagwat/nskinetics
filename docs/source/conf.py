# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'NSKinetics'
copyright = '2025, Sarang S. Bhagwat'
author = 'Sarang S. Bhagwat'

release = ''
version = ''

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']

#
html_theme_options = {
    "logo" : {
        'image_light': '_static/images/logo/logo_nskinetics.png',
        'image_dark': '_static/images/logo/logo_nskinetics_dark.png'
    },
    "show_toc_level": 2,
#     "announcement": (
#         "<p> ..."
#         "<a href='https://...'>...</a></p>"
    ),
#     "external_links": [
#       {"name": "...", "url": "https://..."},

  ]
}

# -- Options for EPUB output
epub_show_urls = 'footnote'

# Autosummary settings
autosummary_generate = True

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

# Napoleon settings
# napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

