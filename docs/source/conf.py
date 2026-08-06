# Configuration file for the Sphinx documentation builder.
import nskinetics

project = 'NSKinetics'
copyright = '2025-2026, Sarang S. Bhagwat'
author = 'Sarang S. Bhagwat'

version = nskinetics.__version__
release = nskinetics.__version__

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx_multitoc_numbering',
    'sphinx_autodoc_typehints',
    'sphinx_design',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_css_files = ['css/custom.css']
html_theme_options = {
    'logo': {
        'image_light': '_static/images/logo/logo_nskinetics_light_white-circle.png',
        'image_dark': '_static/images/logo/logo_nskinetics_dark_white-circle.png',
    },
    'show_toc_level': 2,
}

epub_show_urls = 'footnote'

autosummary_generate = True
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}
# RTD may not install the heavy runtime stack; mock it for autodoc imports.
# numpy/pandas are cheap wheels and are used at import time (e.g. np.inf default
# args), so they are NOT mocked.
autodoc_mock_imports = ['tellurium', 'biosteam', 'thermosteam', 'numba']

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_init_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True
