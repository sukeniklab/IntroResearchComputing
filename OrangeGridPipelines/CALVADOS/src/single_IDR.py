########################################################################################################################
#                                   1. import modules                                                                  #
########################################################################################################################
from utils import *
from rawdata import *
import os
import sys 
from protein_repo import get_ssdomains
import ast
import subprocess
import time
from datetime import timedelta, datetime

########################################################################################################################
#                2. you only need to change parameters below to customize your protein simulation                      #
########################################################################################################################
batch_size = 4
trial = sys.argv[1] 
HOME = os.path.expanduser("~")
cwd = f"{HOME}/CALVADOS/src/"


# GPU BATCH SETTINGS - NEW
check_gpu_memory = True  # Set to False to skip GPU memory check

        
record = sys.argv[2] #"UBQLN2_RTL8"
domain_dict = { ##
    # [[1,233]] this is one domain boundary
    # [[230,240],[250,351]]  # a nested list excludes a small loop (241-249) within domain

    "UBQLN2_RTL8": [[[32,102],[162,260],[300,381]]],
    "UBQLN2_RTL8_modB":  [[[32,102],[162,260],[300,381]]],
    "UBQLN2_UBL-RTL8_bound": [[[32,102],[300,381]],[162,260]],
    "UBQLN2_STI1-RTL8_bound": [[32,102],[[162,260],[300,381]]],
    "RTL8": [[1,82]],
    "RTL8C_FL": [[25,113]],
    "MMACHC": [[1,233]],
    }

temp = 298.15  # unit: K
pH = 7.2
ionic =  0.1 # unit: M
chain_breaks = None   

with open(f"{HOME}/CALVADOS/src/starting/{record}-seq.txt", "r") as file:
    fasta_sequence = ""
    for line in file:
        if not line.startswith(">"):
            fasta_sequence += line.strip()
fasta = list(fasta_sequence)


replicas = 20  # nums of simulation replica for your protein

# is it a IDP or MDP?
isIDP = False
if not isIDP:
    domain_boundaries = {record: domain_dict[record]}
    path2pdb = f'{HOME}/CALVADOS/src/starting/{record}_trial{trial}.pdb'

# customize your desired simulation time, unit: ns;
customized_simulation_time = 100
if customized_simulation_time == None:
    nframes = 300
    discard_first_nframes = 10
else:
    interval = 0.01
purpose = "CALVADOS3"
gpu = True
gpu_id = 0

########################################################################################################################
#          3. no need to change parameters below unless you know exactly what you are doing                            #
########################################################################################################################
initial_type = "C3"
CoarseGrained = "COM"
k_restraint = 700
dataset_replica = 1
cutoff = 2.0
cycle = 0
eps_factor = 0.2
Usecheckpoint = False
slab = False

########################################################################################################################
#                                  LOGGING AND GPU FUNCTIONS - NEW                                                     #
########################################################################################################################
class Logger:
    """Class to handle both console and file logging"""
    def __init__(self, log_file):
        self.log_file = log_file
        self.terminal = sys.stdout
        
    def write(self, message):
        """Write to both console and file"""
        self.terminal.write(message)
        with open(self.log_file, 'a') as f:
            f.write(message)
    
    def flush(self):
        """Flush both outputs"""
        self.terminal.flush()

def log_print(message, log_file=None):
    """Print to console and optionally to log file"""
    print(message)
    if log_file:
        with open(log_file, 'a') as f:
            f.write(message + '\n')

def get_gpu_memory():
    try:
        result = subprocess.run(['nvidia-smi', '--query-gpu=memory.total,memory.used,memory.free', 
                                '--format=csv,noheader,nounits'], 
                               capture_output=True, text=True)
        total, used, free = map(int, result.stdout.strip().split(','))
        utilization = (used/total)*100
        return {
            'total': total,
            'used': used,
            'free': free,
            'utilization': utilization
        }
    except Exception as e:
        return None

def format_gpu_info(gpu_info):
    """Format GPU info for display"""
    if gpu_info is None:
        return "GPU Memory: Could not query"
    return (f"GPU Memory: {gpu_info['used']}MB used / {gpu_info['total']}MB total "
            f"({gpu_info['free']}MB free, {gpu_info['utilization']:.1f}% utilized)")

def format_time(seconds):
    """Format seconds into human-readable time"""
    return str(timedelta(seconds=int(seconds)))

########################################################################################################################
#                                            4. submit simulations                                                     #
########################################################################################################################
# Start overall timer and setup logging
overall_start_time = time.time()
start_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

dataset = f"{purpose}{CoarseGrained}_{cutoff}_MD_gpu_trial{trial}_{record}"
print(dataset)

if not os.path.isdir(f"{cwd}/{dataset}"):
    os.system(f"mkdir -p {cwd}/{dataset}")

# Create log file
log_file = f"{cwd}/{dataset}/simulation_log_trial{trial}_{record}.txt"
log_print(f"{'='*70}", log_file)
log_print(f"CALVADOS SIMULATION LOG", log_file)
log_print(f"{'='*70}", log_file)
log_print(f"Start time: {start_datetime}", log_file)
log_print(f"Record: {record}", log_file)
log_print(f"Trial: {trial}", log_file)
log_print(f"Dataset: {dataset}", log_file)
log_print(f"Replicas: {replicas}", log_file)
log_print(f"Batch size: {batch_size}", log_file)
log_print(f"Simulation time: {customized_simulation_time} ns", log_file)
log_print(f"Temperature: {temp} K", log_file)
log_print(f"pH: {pH}", log_file)
log_print(f"Ionic strength: {ionic} M", log_file)
log_print(f"Sequence length: {len(fasta)} residues", log_file)
log_print(f"{'='*70}\n", log_file)

os.system(f"cp {cwd}/residues_pub.csv {cwd}/{dataset}")
create_parameters(cwd, dataset, cycle, initial_type)

name = record.split("@")[0]
L = int(np.ceil((len(fasta) - 1) * 0.38 + 4))
fdomains = f'{cwd}/{dataset}/domain_simple.yaml'
print(fdomains)

if not isIDP:
    yaml.dump(domain_boundaries, open(fdomains, 'w'))
    input_pae = ""
    use_pdb = True
    use_hnetwork = True
    use_ssdomains = True
    if customized_simulation_time == None:
        domain_len = 0
        for domain in get_ssdomains(name, fdomains, output=False):
            domain_len += len(domain)
        N_res = int(len(fasta) - domain_len)
        N_save = 3000 if N_res < 100 else int(np.ceil(3e-4 * N_res ** 2) * 1000)
        N_steps = (nframes + discard_first_nframes) * int(N_save)
    else:
        N_steps = int(customized_simulation_time * 1E5)
        total_frames = int(customized_simulation_time / interval)
        discard_first_nframes = int(total_frames * 0.05)
        nframes = int(total_frames - discard_first_nframes)
        N_save = int(interval * 1E5)
else:
    input_pae = None
    path2pdb = ""
    use_pdb = False
    use_hnetwork = False
    use_ssdomains = False
    if customized_simulation_time == None:
        N_save = 3000 if len(fasta) < 100 else int(np.ceil(3e-4 * len(fasta) ** 2) * 1000)
        N_steps = (nframes + discard_first_nframes) * int(N_save)
    else:
        N_steps = int(customized_simulation_time * 1E5)
        total_frames = int(customized_simulation_time / interval)
        discard_first_nframes = int(total_frames * 0.05)
        nframes = int(total_frames - discard_first_nframes)
        N_save = int(interval * 1E5)

if not os.path.isdir(f"{cwd}/{dataset}/{record}/{cycle}"):
    os.system(f"mkdir -p {cwd}/{dataset}/{record}/{cycle}")

replicas_list4MD = list(range(replicas))
config_sim_data = dict(cwd=cwd, name=name, dataset=dataset, temp=temp, ionic=ionic, cycle=cycle, pH=pH,
   replicas_list4MD=replicas_list4MD, cutoff=cutoff, L=L, wfreq=int(N_save), slab=slab,
   use_pdb=use_pdb, path2pdb=path2pdb, use_hnetwork=use_hnetwork, fdomains=fdomains,
   use_ssdomains=use_ssdomains, input_pae=input_pae, k_restraint=k_restraint, record=record,
   gpu_id=gpu_id, Threads=1, overwrite=True, N_res=len(fasta),
   CoarseGrained=CoarseGrained, isIDP=isIDP, Usecheckpoint=Usecheckpoint, eps_factor=eps_factor,
   initial_type=initial_type, seq=fasta, steps=N_steps, gpu=gpu, replicas=replicas,
   discard_first_nframes=discard_first_nframes,validate=False, nframes=nframes,chain_breaks=chain_breaks)

# Check GPU memory before starting
if check_gpu_memory and gpu:
    log_print("\n=== Initial GPU Memory Status ===", log_file)
    initial_gpu = get_gpu_memory()
    log_print(format_gpu_info(initial_gpu), log_file)
    log_print(f"Running {replicas} replicas in batches of {batch_size}", log_file)
    log_print(f"Total batches: {int(np.ceil(replicas / batch_size))}\n", log_file)

# Initialize Ray with GPU support
if gpu:
    ray.init(num_cpus=batch_size, num_gpus=1, include_dashboard=False, _temp_dir=f"/tmp/")
else:
    ray.init(num_cpus=len(replicas_list4MD), include_dashboard=False, _temp_dir=f"/tmp/")

# Pre-create and save all config files
log_print("Creating configuration files...", log_file)
config_start = time.time()
configs = []
for replica in replicas_list4MD:
    config_sim_data["replica"] = replica
    config_path = f"{cwd}/{dataset}/{record}/{cycle}/config_{replica}.yaml"
    yaml.dump(config_sim_data, open(config_path, 'w'))
    configs.append(yaml.safe_load(open(config_path, 'r')))
config_time = time.time() - config_start
log_print(f"Config files created in {config_time:.2f} seconds\n", log_file)

# Run simulations in batches with timing
log_print("="*70, log_file)
log_print("Starting batched simulation runs...", log_file)
log_print("="*70, log_file)
simulation_start_time = time.time()
batch_times = []
gpu_stats = []

for i in range(0, len(replicas_list4MD), batch_size):
    batch_replicas = replicas_list4MD[i:i+batch_size]
    batch_configs = configs[i:i+batch_size]
    batch_num = i//batch_size + 1
    total_batches = int(np.ceil(len(replicas_list4MD) / batch_size))
    
    log_print(f"\n{'='*70}", log_file)
    log_print(f"BATCH {batch_num}/{total_batches}: Running replicas {batch_replicas}", log_file)
    log_print(f"Batch start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", log_file)
    log_print(f"{'='*70}", log_file)
    
    batch_start = time.time()
    ray.get([simulate_simple.remote(config) for config in batch_configs])
    batch_elapsed = time.time() - batch_start
    batch_times.append(batch_elapsed)
    
    log_print(f"\nBatch {batch_num} completed in {format_time(batch_elapsed)}", log_file)
    log_print(f"Average time per replica in this batch: {batch_elapsed/len(batch_replicas):.2f} seconds", log_file)
    
    if check_gpu_memory and gpu:
        gpu_info = get_gpu_memory()
        gpu_stats.append(gpu_info)
        log_print(f"\nGPU status after batch {batch_num}:", log_file)
        log_print(format_gpu_info(gpu_info), log_file)
    
    # Show progress
    completed_replicas = min((batch_num) * batch_size, replicas)
    avg_batch_time = np.mean(batch_times)
    remaining_batches = total_batches - batch_num
    estimated_remaining = avg_batch_time * remaining_batches
    
    log_print(f"\nProgress: {completed_replicas}/{replicas} replicas completed", log_file)
    if remaining_batches > 0:
        log_print(f"Estimated time remaining: {format_time(estimated_remaining)}", log_file)

total_simulation_time = time.time() - simulation_start_time
log_print("\n" + "="*70, log_file)
log_print("All simulations complete!", log_file)
log_print("="*70, log_file)
log_print(f"Total simulation time: {format_time(total_simulation_time)}", log_file)
log_print(f"Average time per batch: {format_time(np.mean(batch_times))}", log_file)
log_print(f"Average time per replica: {total_simulation_time/replicas:.2f} seconds", log_file)
log_print(f"Fastest batch: {format_time(min(batch_times))}", log_file)
log_print(f"Slowest batch: {format_time(max(batch_times))}", log_file)
log_print("="*70 + "\n", log_file)

# Write detailed batch statistics
log_print("\n" + "="*70, log_file)
log_print("DETAILED BATCH STATISTICS", log_file)
log_print("="*70, log_file)
for i, (batch_time, replicas_batch) in enumerate(zip(batch_times, 
                                                       [replicas_list4MD[j:j+batch_size] 
                                                        for j in range(0, len(replicas_list4MD), batch_size)])):
    log_print(f"\nBatch {i+1}:", log_file)
    log_print(f"  Replicas: {replicas_batch}", log_file)
    log_print(f"  Time: {format_time(batch_time)}", log_file)
    log_print(f"  Time per replica: {batch_time/len(replicas_batch):.2f} seconds", log_file)
    if i < len(gpu_stats) and gpu_stats[i]:
        log_print(f"  GPU Memory Used: {gpu_stats[i]['used']}MB ({gpu_stats[i]['utilization']:.1f}%)", log_file)
log_print("="*70 + "\n", log_file)

# Merging trajectories
log_print("Merging trajectories...", log_file)
merge_start = time.time()
config_merge_data = dict(cwd=cwd, name=name, dataset=dataset, temp=temp, ionic=ionic,
                         cycle=cycle, pH=pH, replicas_list4MD=replicas_list4MD, cutoff=cutoff, L=L,
                         wfreq=int(N_save), slab=slab, use_pdb=use_pdb, path2pdb=path2pdb,
                         use_hnetwork=use_hnetwork, fdomains=fdomains, use_ssdomains=use_ssdomains,
                         input_pae=input_pae, k_restraint=k_restraint, record=record, gpu_id=gpu_id, Threads=1,
                         overwrite=True, N_res=len(fasta), CoarseGrained=CoarseGrained, isIDP=isIDP,
                         Usecheckpoint=Usecheckpoint, eps_factor=eps_factor, initial_type=initial_type, seq=fasta,
                         steps=N_steps, gpu=gpu, replicas=replicas, discard_first_nframes=discard_first_nframes,
                         validate=False, nframes=nframes)
centerDCD_simple(config_merge_data)
merge_time = time.time() - merge_start
log_print(f"Trajectories merged in {format_time(merge_time)}\n", log_file)

# Calculate simulated Rg
log_print("Calculating Rg...", log_file)
rg_start = time.time()
df = load_parameters(cwd, dataset, cycle, initial_type).set_index("three")
t = md.load_dcd(f"{cwd}/{dataset}/{record}/{cycle}/{record}.dcd",
                f"{cwd}/{dataset}/{record}/{cycle}/{record}.pdb")
residues = [res.name for res in t.top.atoms]
masses = df.loc[residues,'MW'].values
masses[0] += 2
masses[-1] += 16

# Optimized Rg calculation
masses_normalized = masses / masses.sum()
cm = np.einsum('ijk,j->ik', t.xyz, masses) / masses.sum()
si_squared = np.sum((t.xyz - cm[:, np.newaxis, :])**2, axis=2)
rgarray = np.sqrt(np.einsum('ij,j->i', si_squared, masses_normalized))
np.save(f"{cwd}/{dataset}/{record}/{cycle}/Rg_traj.npy", rgarray)
rg_time = time.time() - rg_start
log_print(f"Rg calculation completed in {rg_time:.2f} seconds\n", log_file)

# Final summary
total_time = time.time() - overall_start_time
end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

log_print("="*70, log_file)
log_print("FINAL TIMING SUMMARY", log_file)
log_print("="*70, log_file)
log_print(f"Configuration setup:    {format_time(config_time)}", log_file)
log_print(f"Simulations:            {format_time(total_simulation_time)}", log_file)
log_print(f"Trajectory merging:     {format_time(merge_time)}", log_file)
log_print(f"Rg calculation:         {format_time(rg_time)}", log_file)
log_print(f"{'─'*70}", log_file)
log_print(f"TOTAL RUNTIME:          {format_time(total_time)}", log_file)
log_print(f"\nEnd time: {end_datetime}", log_file)
log_print("="*70, log_file)

print(f"\nLog file saved to: {log_file}")
