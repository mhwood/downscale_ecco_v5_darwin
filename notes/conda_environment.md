## Python Packages
To use the python scripts provided in this directory, several packages are required. 

These steps were successfully tested on 10 Jan 2023 on a MacBook.

It is recommended that a fresh conda environment is established to include these packages. To create a new environment, use
```
conda create --name mitgcm
conda activate mitgcm
```
Then, with the package activated, first install the following packages:
```
conda install numpy
conda install matplotlib
conda install scipy
conda install pyproj
conda install netcdf4
conda install -c conda-forge xesmf
```
Next, install the required [simplegrid](https://github.com/nasa/simplegrid) package, which is not available by `pip` or `conda install`. Instead, it must be cloned and then installed locally: 
```
git clone https://github.com/nasa/simplegrid.git
cd simplegrid
pip install .
```
Next, we will install the [ecco_v4_py](https://github.com/ECCO-GROUP/ECCOv4-py). While a `pip` install of `ecco-v4-py` is available (`pip install ecco-v4-py`) I have run into issues with it. The `conda` implementation also seems not to work. Instead, i find it easiest manually generate a package in my conda environment site-packages as follows:
```
git clone https://github.com/ECCO-GROUP/ECCOv4-py.git
mkdir [conda dir]/envs/mitgcm2/lib/python3.10/site-packages/ecco_v4_py
cp ECCOv4-py/ecco_v4_py/* [conda dir]/envs/mitgcm2/lib/python3.10/site-packages/ecco_v4_p
```
This package will in turn require a few more dependencies, installed as follows:
```
conda install -c conda-forge xgcm
conda install -c conda-forge xmitgcm
conda install -c conda-forge pyresample
conda install -c conda-forge cartopy
```

Finally, add the MITgcm utils to your environment (from your local clone of MITgcm):
```
cd /path/to/MITgcm/utils/python/MITgcmutils
python setup.py install
```
Alternatively, you can copy this into your environment site-packages as was done above for `ecco_v4_py`.
