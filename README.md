# Welcome to a brief overview of research computing 

This repository contains working pipelines and a detailed [Wiki](https://github.com/sukeniklab/IntroResearchComputing/wiki) for those looking to start using Syracuse University's HTC and HPC resources. 

If you're not sure whether you even need any of this, start with the [Home page](../../wiki/Home), it has a short checklist for "when should I start using computing resources?"

# Where to start
If you are looking to get started, look at different examples to figure out what suits your needs and pull one of the `.sub` files.

# Example pipelines

The [`OrangeGridPipelines/`](OrangeGridPipelines) folder has complete, working `.sub` and `.sh` files you can copy directly and adjust for your own data. 
  - **[`esm_embed/`](OrangeGridPipelines/esm_embed)** — generating ESM2 protein embeddings, one job per FASTA file
  - **[`CALVADOS/`](OrangeGridPipelines/CALVADOS)** — sweeping multiple variables (e.g. replicates × constructs) from a file


# Getting help 
If your job isn't getting picked up, try  `condor_q -better-analyze <jobid>` to see if any nodes match your specs, or check the `.log` file if everything instantly dies. For any other issues, send an email to Research Computing at **researchcomputing@syr.edu**. 
