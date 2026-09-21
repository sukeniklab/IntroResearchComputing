
## Embedding sequences with ESM

### 1. Set up your directory that you're going to work in 
```bash
mkdir esm_embed
```
*** examples in this directory assume the path /home/{username}/esm_embed ***

---


### 2. Create and activate the conda environment
```bash
cd esm_embed 
conda create -n esm2 -c pytorch -c conda-forge python=3.9 pytorch
conda activate esm2
pip install fair-esm
git clone https://github.com/facebookresearch/esm.git
```
This will pull a directory from Meta that contains the scripts necessary to run ESM . There is also a .yml file if this command does not work.

---


### 3. Prep fasta files
```bash
mkdir fastas
mkdir fastas/to_run/
```
Put a fasta file (any number of sequences) in fastas

```bash
cd fastas
python split_fasta.py
```
split_fasta.py will break your fasta into 50 chunks (could be even more chunks) that are written to_run/

---


### 4. Submit the job 
```bash
mkdir logs
chmod +x run_esm2_to_embed.sh 
```
You should check to make sure all paths look good in both run_esm2_to_embed.sh and submit_esm_job.sub.  If paths look good, run: 

```bash
condor_submit submit_esm_job.sub
```

---

### 5. Output 

Embeddings will be saved in output_pt/{fasta_name}, with a .pt file saved per sequence for further analysis.