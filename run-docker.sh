#!/bin/bash -ie

pushd example-data
./download.sh
popd

output_folder='output_notebooks'
input_folder='input_notebooks'
reports_folder='reports'

mkdir $output_folder $input_folder $reports_folder

for pyfile in jovian-01-upstream-qc.py jovian-02a-downstream-qc.py jovian-02b-downstream-integrated.py jovian-03-DE.py; do
    ipynbfile=${pyfile/.py/.ipynb}
    jupytext --to notebook -o $input_folder/$ipynbfile $pyfile
    papermill -f parameters.yaml $input_folder/$ipynbfile $output_folder/${ipynbfile/.ipynb/-output.ipynb}
    jupyter nbconvert --output ${ipynbfile/.ipynb/.html} --output-dir $reports_folder --to html $output_folder/$ipynbfile
    jupyter nbconvert --output ${ipynbfile/.ipynb/.pdf} --output-dir $reports_folder --to pdf $output_folder/$ipynbfile
done
