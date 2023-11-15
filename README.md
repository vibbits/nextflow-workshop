# Nextflow workshop

This repository contains the data, scripts, documentation and relevant information for the website containing the training materials. If you would like to get started with the materials, please visit the website: [vibbits-nextflow-workshop.rtfd.io](https://vibbits-nextflow-workshop.rtfd.io).

## Website

The website was generated with Sphinx read-the-docs with Furo template.

## Development

### Dependencies

```
python>=3.8
sphinx-autobuild
sphinx-copybutton
sphinx-inline-tabs
myst-parser
furo
```

```
conda create -n workshop_docs python>=3.11 -c conda-forge sphinx-autobuild sphinx-copybutton sphinx-inline-tabs myst-parser furo

conda activate workshop_docs
```

Local test build with:

```
sphinx-autobuild docs docs/_build/html
```

Then go to: http://127.0.0.1:8000/index.html

[![Documentation Status](https://readthedocs.org/projects/rtd-bioinformatics/badge/?version=latest)](https://rtd-bioinformatics.readthedocs.io/en/latest/?badge=latest)
