# Documentation

The ready-to-read documentation is available at https://decneo.readthedocs.io/.
The documentation of our software is built with [**Sphinx**](https://www.sphinx-doc.org/ "Sphinx") at 
[**ReadTheDocs.org**](https://readthedocs.org/).

- [Documentation Updates](#documentation-updates)
- [Documentation Development](#documentation-development)

## Documentation Updates

Any changes to source code and docstrings
are automatically reflected at [**ReadTheDocs.org**](https://readthedocs.org/) 
(a new version of documentation is built silently). 

## Documentation Development

For development and testing of the documentation locally (on the development machine) 
install [**Sphinx**](https://www.sphinx-doc.org/ "Sphinx") by:

	$ pip install --user sphinx
	$ pip install sphinx-press-theme

To compile **html** version of the documentation go to **docs/** directory and run:

	$ sphinx-build -E -a -b html ./docs/source ./docs/build

We are utilizing a 3rd party **Sphinx** extension [**sphinxcontrib-images**](https://github.com/sphinx-contrib/images 
"GitHub repository") extension, allowing to display documentation images in a organized way.
