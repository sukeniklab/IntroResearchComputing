## Getting Nardini parameters for multiple sequences at once

### 1. Set up your directory that you're going to work in 
```bash
mkdir batch_nardini
```
*** examples in this directory assume the path /home/{username}/batch_nardini ***

---

### 2. Create your conda environment
```bash
cd batch_nardini 
conda env create -f environment.yml
conda install -c conda-forge gxx_linux-64=12 gcc_linux-64=12 cmake ninja -y
CC=$(dirname $(which python))/x86_64-conda-linux-gnu-cc 
CXX=$(dirname $(which python))/x86_64-conda-linux-gnu-c++ 
conda deactivate
conda activate nardini_features
pip install localcider
pip install git+https://github.com/idptools/sparrow.git
```
You should be able to do `$CXX --version` and get 12.4.0, if not talk to claude about updating g++

### 3. Generate your job list 
Place any fasta that you want to analyze in fastas/, and then run: 
```bash
python generate_job_list.py fastas/homo_sapiens_short.fasta job_list.txt
```
*** Note, OrangeGrid has a max job number of 15,000 (from my experience). You may need to break up very large jobs. *** 

### 4. Submit jobs 
```bash
mkdir logs
chmod +x run_calc_features_single.sh
condor_submit sub_nardini_calc.sub
```
This will generate one csv per sequence in `output_csvs`

### 4. Pull sequence-level csvs to master csv
```bash
 python combine_csvs.py output_csvs/ combined_sequence_features.csv
```


