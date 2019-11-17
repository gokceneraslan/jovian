#!/bin/bash -ie

pushd data
./download.sh
popd

output_folder='output_notebooks'
mkdir $output_folder

for notebook in jovian-01-upstream-qc.ipynb jovian-02a-downstream-qc.ipynb jovian-02b-downstream-integrated.ipynb jovian-03-DE.ipynb; do
    papermill $notebook $output_folder/${notebook/.ipynb/-output.ipynb}
done
