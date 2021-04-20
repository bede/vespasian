#!/bin/bash
python setup.py sdist
twine upload dist/* --skip-existing