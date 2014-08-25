Virtual Machine IPython Notebooks
=================================

Data and notebooks for the following courses:
- Climatologia e Meteorologia
- Oceanografia Física Descritiva
- Oceanografia Física Dinâmica
- Ondas e Marés

This [virtual machine](https://susestudio.com/a/YfJVDT/python4oceanographers--2)
open IPython notebooks as soon as it starts.  To get some help navigating them
watch [this](http://www.youtube.com/watch?v=QSmmt2hDIwY).

Some of the notebooks will not work because iris/cartopy are not python3
compatible yet.

To be able to run those we need to crate an iris virtual environment first.

mkvirtualenv iris
workon iris
pip install -r iris_env.txt
deactivate
