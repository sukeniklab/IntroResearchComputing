## Running multi-protein simulations with CALVADOS

### 1. Set up your directory that you're going to work in 
```bash
mkdir CALVADOS
```
*** examples in this directory assume the path /home/{username}/CALVADOS ***

---

### 2. Create your conda environments
First make the modeller environment to generate multiple starting structures:
```bash
cd CALVADOS/init_AD_struc 
conda env create -f environment.yml
```

Then, build the CALVADOS package (this may take a sec):
```bash
cd ../src/ 
conda env create -f environment.yml
```

---

### 3. Generate multiple starting conformations
Best practices usually run simulations in triplicate. Place as many AlphaFold structures in this directory, and submit the jobs. 

```bash
cd ../init_AD_struc
mkdir logs
conda activate modeller_env
condor_submit run_modeller.sub
```

This will generate 3 (modify in generate_initStructure.py, line 28) initial conformations and move them to ../src/starting

---

### 4. Start calvados simulations 
All simulations are run in the /src/ folder.

```bash
cd ../src
mkdir logs
conda activate CALVADOS3
```
Check to make sure the starting structures are present: 

```bash
ls starting/
```
make sure that the domains for your protein are included in single_IDR.py: 
```bash
vim single_IDR.py
```
modify domain_dict = {...}, ~line 28 to add in your sequence and change any temperature, pH, or ionic strength. To close vim, pres esc and type:
 
```bash
:wq 
```

If it all looks good, submit the job 
```bash
condor_submit submit_calvados_job.sub
```

---

### 5. Output 
Simulations will be in CALVADOS3COM_2.0_MD_gpu_trial{trial}_{sequence_name} located in the src directory.