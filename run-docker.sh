#!/bin/bash -ie

pushd data
./download.sh
popd

output_folder='output_notebooks'
mkdir $output_folder

for pyfile in jovian-01-upstream-qc.py jovian-02a-downstream-qc.py jovian-02b-downstream-integrated.py jovian-03-DE.py; do
    jupytext --to notebook $pyfile
    ipynbfile=${pyfile/.py/.ipynb}
    papermill $ipynbfile $output_folder/${ipynbfile/.ipynb/-output.ipynb}
done
