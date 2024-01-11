To build the documentation you will need the following packages:

- sphinx
- sphinx-apidoc
- autodoc

### Steps
1. First, generate the autdoc sources for the API reference from the root using

        sphinx-apidoc -f -o docs asimtools

2. Then go into the docs directory and run

        make html

This will generate the docs in the `_build` directory


