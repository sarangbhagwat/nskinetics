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
# conf.py imports nskinetics above, so any environment that can build these
# docs already has the full runtime stack importable (locally and on RTD,
# which installs requirements.txt plus the package itself) — autodoc
# therefore documents real modules and real docstrings, notably the
# processes system factory (a biosteam.SystemFactory instance whose full
# __doc__ exists only with real biosteam). The conditional list below is a
# fallback that mocks only modules genuinely missing; it is inert when
# everything is installed. numpy/pandas are cheap wheels used at import
# time (e.g. np.inf default args), so they are never mocked.
import importlib.util
autodoc_mock_imports = [
    m for m in ('tellurium', 'biosteam', 'thermosteam', 'numba')
    if importlib.util.find_spec(m) is None
]

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_init_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True
