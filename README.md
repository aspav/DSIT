# DSIT: Data Science & Information Technologies

[![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/vagmcs/prml/tree/master/)
[![License: LGPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/python-3.8-3670A0?style=flat&logo=python&logoColor=ffdd54)]()
[![Jupyter](https://img.shields.io/badge/Jupyter-white?style=flat&logo=Jupyter)](https://jupyter.org/try)
[![Matlab](https://img.shields.io/badge/MATLAB-9.4--2018a-red?style=flat&logo=)]()

This repository contains the projects that were developed in the scope of the courses of the MSc program [Data Science and Information Technologies-Bioinformatics & Biomedical Data Science](http://dsit.di.uoa.gr/).

## ÂŠī¸ License

This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions; See the [GNU General Public License v3 for more details](http://www.gnu.org/licenses/gpl-3.0.en.html).

## đ AIMB

##### âī¸ Algorithms in Molecular Biology

This folder contains the code and reports for the assignments of the course **Algorithms in Molecular Biology** as well as the data that are required to run the notebooks.

##### đ data đī¸

##### đ notebooks

* đž [frequent_words.ipynb](https://nbviewer.org/github/aspav/DSIT/blob/main/AIMB/notebooks/frequent_words.ipynb)

* đž [motif_search_algorithms.ipynb](https://nbviewer.org/github/aspav/DSIT/blob/main/AIMB/notebooks/motif_search_algorithms.ipynb)

* đž cutoff_heuristic (to be uploaded)

* đž miRNA_BUFET (to be uploaded)

## đ AISB

##### âī¸ Algorithms in Structural Bioinformatics

This project contains four notebooks that were developed in the scope of the course **Algorithms in Structural Bioinformatics**.

##### đ data đī¸

##### đ notebooks

* đž [Implementation of c-RMSD & d-RMSD](https://nbviewer.org/github/aspav/DSIT/blob/main/AISB/notebooks/cRMSD_dRMSD.ipynb)

* đž [Implementation of a simplified Zucker minimization algorithm](https://nbviewer.org/github/aspav/DSIT/blob/main/AISB/notebooks/zucker_minimization_algorithm.ipynb)

* đž [Construction of Cayley-Menger matrix B-rank(B) computation-perturbation](https://nbviewer.org/github/aspav/DSIT/blob/main/AISB/notebooks/distance_geometry.ipynb)

* đž [Playing with pdb files](https://nbviewer.org/github/aspav/DSIT/blob/main/AISB/notebooks/playing_with_pdb_files.ipynb)

##### đ reports

* đ playing_with_pdb_files.pdf

## đ DSITNeuro

##### âī¸ Application of Data Science and Information Technologies in Neurosciences

This folder contains the code (implemented in `matlab`) and the reports (*in pdf format*) for the following assignments of the course:

##### đ **IntegrateFireNeuron**: Integrate and Fire Neuron Model f-I curve

* đž `IntegrateFireNeuron.m`

* đ IntegrateFireNeuron.pdf

##### đ **HodgkinHuxley**: The Hodgkin-Huxley Model as an Oscillator

* đž `HodgkinHuxley.m`

* đ HodgkinHuxley.pdf

## 

## đ MLCB

### âī¸ Machine Learning in Computational Biology

##### đŋ Supervised Machine Learning

* đ output

* đž `ex_3.py`

* đž `ex_4.py`

* đž `ex_5.py`

* đ Supervised_Learning.pdf 

##### đŋ Unsupervised Machine Learning

* đ datasets đī¸

* đž [Unsupervised_Learning.ipynb](https://nbviewer.org/github/aspav/DSIT/blob/main/MLCB/Unsupervised%20Learning/Unsupervised_Learning.ipynb)

* đ Unsupervised_Learning.pdf

## âī¸ HOW TO

### đ Regarding the `python` files and notebooks, the next steps should be followed to run them:

**âī¸ī¸ NOTE: I used Python 3.8 for running the following commands.** 

1. Install the `virtualenv` package:

```python
pip install virtualenv
```

2. Create a virtual environment

```python
virtualenv venv
```

3. Activate the virtual environment

```python
. vnenv/bin/activate
```

4. Install the required packages, provided in the `requirements.txt` file:

```python
pip install -r requirements.txt
```

5. Running the files.
* To run a `jupyter notebook` type:
  
  ```python
  jupyter notebook
  ```

* To run the `.py` files type:
  
  ```python
  python3 *.py
  ```

### đ§Ž Regarding the `matlab` files, you should download MATLAB 9.4 2018a.
