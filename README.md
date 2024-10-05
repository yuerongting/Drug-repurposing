The code for the paper "Repurposing Drugs for Infectious Diseases by Graph Convolutional Network wish Sensitivity-based Graph Reduction".

**Dependencies**
Pandas version: 2.2.2
Matplotlib version: 3.9.0
Scikit-learn version: sklearn 2.2.2
NumPy version: 1.23.5
PyTorch version: 2.3.0
Torch Geometric version: 2.3.0



**The main algorithm of the sensitivity analysis presented in our paper is implemented below:**
For COVID-19: "sensitivity_covid.py";
For Zika virus: "sensitivity_Zika_GPU.py";


**For Zika virus drug prediction:**
"hetero_Zika_retrain_new_grid_search.py" for the retrained model with grid search;
"hetero_Zika_retrain_same_setting.py" for the retrained model with the same settings before graph reduction.

**For COVID-19 drug prediction:**
"hetero_COVID_retrain_new_grid.py" for the retrained model with grid search;
"hetero_COVID_retrain_same.py" for the retrained model with the same settings before graph reduction.


All data is in the folder "Data".
