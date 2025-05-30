# Sphinx Configuration for CLI Enhanced Validation Engine API Documentation
# Based on Phase 5 research findings and industry best practices

import os
import sys
from pathlib import Path

# Path setup
project_root = Path(__file__).parent.parent.parent  # Go up to VCF_Agent root
sys.path.insert(0, str(project_root / "src"))
sys.path.insert(0, str(project_root / "scripts"))

# Project information
project = 'CLI Enhanced Validation Engine'
copyright = '2025, VCF Agent Development Team'
author = 'VCF Agent Development Team'
version = '1.0.0'
release = '1.0.0'

# General configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx_autodoc_typehints',
    'autoapi.extension',
]

# Templates path
templates_path = ['_templates']

# List of patterns to exclude when looking for source files
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Language for content autogeneration
language = 'en'

# -- Options for HTML output -------------------------------------------------

# HTML theme
html_theme = 'sphinx_rtd_theme'

# Theme options
html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#2980B9',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Custom sidebar templates
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
        'donate.html',
    ]
}

# Static files path
html_static_path = ['_static']

# Custom CSS files
html_css_files = [
    'custom.css',
]

# -- Extension configuration -------------------------------------------------

# AutoAPI configuration
autoapi_dirs = [
    str(project_root / "scripts"),  # For cli_enhanced_validation.py
]
autoapi_type = 'python'
autoapi_template_dir = '_templates/autoapi'
autoapi_options = [
    'members',
    'undoc-members',
    'show-inheritance',
    'show-module-summary',
]
autoapi_python_class_content = 'both'
autoapi_member_order = 'groupwise'
autoapi_root = 'api'
autoapi_keep_files = True

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
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

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'special-members': '__init__',
    'exclude-members': '__weakref__'
}
autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented'
autodoc_preserve_defaults = True

# Autosummary settings
autosummary_generate = True
autosummary_generate_overwrite = True

# Intersphinx configuration for external documentation
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'click': ('https://click.palletsprojects.com/en/stable/', None),
}

# Coverage settings
coverage_show_missing_items = True 