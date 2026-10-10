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
    # Search as an icon button (Ctrl+K still opens it) instead of the
    # 175-px search field: with five header links the field made the links
    # wrap onto two lines, doubling the header height, at ~960-1150 px.
    'navbar_persistent': ['search-button'],
    'icon_links': [
        {
            'name': 'GitHub',
            'url': 'https://github.com/sarangbhagwat/nskinetics',
            'icon': 'fa-brands fa-github',
            'type': 'fontawesome',
        },
    ],
}

# Case-study pages carry dense 1600-px publication figures and wide numeric
# tables. On those pages the "On this page" contents move from the right
# (secondary) sidebar into the left sidebar, under the section navigation,
# so the article gets the full content width (about 920 px instead of
# 645 px in a 1280-px window). Every other page keeps the default layout.
# The secondary-sidebar overrides name each page explicitly: pydata warns
# when a page matches two wildcard patterns ('**' and 'case_studies/**').
import pathlib
_DOCS_SOURCE = pathlib.Path(__file__).resolve().parent
_WIDE_PAGES = sorted(
    p.relative_to(_DOCS_SOURCE).with_suffix('').as_posix()
    for p in (_DOCS_SOURCE / 'case_studies').rglob('*.rst')
)
html_sidebars = {
    'case_studies/**': ['sidebar-nav-bs.html', 'page-toc.html', 'sourcelink.html'],
}
html_theme_options['secondary_sidebar_items'] = {
    '**': ['page-toc', 'edit-this-page', 'sourcelink'],
    **{page: [] for page in _WIDE_PAGES},
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
# The second group are the optional dependencies of the isobutanol
# biorefinery case study (the ``isobutanol`` extra) and of the shipped
# model; Read the Docs does not install them, so they are mocked there and
# the case-study API page (API/biorefineries.rst) can still be autodoc'd.
# Importing nskinetics.biorefineries.isobutanol is build-free; never call
# its load() from this file or from a docs page.
import importlib.util
autodoc_mock_imports = [
    m for m in ('tellurium', 'biosteam', 'thermosteam', 'numba',
                'biorefineries', 'hensmith', 'contourplots', 'optuna',
                'SALib', 'openpyxl', 'sklearn', 'simplesbml', 'torch')
    if importlib.util.find_spec(m) is None
]

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_init_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True


# Typography: keep a number and its percent sign on one line. The pages
# write "12.60 %" with a space, which narrow table cells otherwise break
# before the "%". Text inside literals, code blocks, math and raw HTML is
# left untouched.
import re
from docutils import nodes as _nodes

_NUMBER_SPACE_PERCENT = re.compile(r'(?<=\d) (?=%)')
_NO_TYPOGRAPHY = (_nodes.literal, _nodes.literal_block, _nodes.doctest_block,
                  _nodes.math, _nodes.math_block, _nodes.raw, _nodes.comment)


def _join_numbers_to_percent_signs(app, doctree):
    for text in list(doctree.findall(_nodes.Text)):
        if not _NUMBER_SPACE_PERCENT.search(text):
            continue
        parent = text.parent
        node = parent
        while node is not None and not isinstance(node, _NO_TYPOGRAPHY):
            node = node.parent
        if node is None and parent is not None:
            parent.replace(text, _nodes.Text(
                _NUMBER_SPACE_PERCENT.sub('\u00a0', str(text))))


# Keep a numeric range, and a time unit right after it, on one line:
# narrow table cells otherwise break "hours (~7\u201311 h)" after the en dash
# ("~7\u2013" / "11 h)"). Each range is wrapped in an inline node with the
# ``ibo-nobr`` class (``white-space: nowrap`` in custom.css), so the text,
# the search index and copy-paste are unchanged. Same exclusions as above.
_NUMBER_RANGE = re.compile(
    r'~?\d[\d.,]*\u2013\d[\d.,]*(?:[ \u00a0](?:h|min|s)\b)?')


def _keep_number_ranges_together(app, doctree):
    for text in list(doctree.findall(_nodes.Text)):
        s = str(text)
        if not _NUMBER_RANGE.search(s):
            continue
        parent = text.parent
        if parent is None or 'ibo-nobr' in parent.get('classes', ()):
            continue
        node = parent
        while node is not None and not isinstance(node, _NO_TYPOGRAPHY):
            node = node.parent
        if node is not None:
            continue
        pieces, pos = [], 0
        for m in _NUMBER_RANGE.finditer(s):
            if m.start() > pos:
                pieces.append(_nodes.Text(s[pos:m.start()]))
            pieces.append(_nodes.inline(m.group(0), m.group(0),
                                        classes=['ibo-nobr']))
            pos = m.end()
        if pos < len(s):
            pieces.append(_nodes.Text(s[pos:]))
        parent.replace(text, pieces)


def setup(app):
    app.connect('doctree-read', _join_numbers_to_percent_signs)
    app.connect('doctree-read', _keep_number_ranges_together)
