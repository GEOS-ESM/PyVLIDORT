#!/usr/bin/env python 

"""
  This script is expected to be run from the AeroApps/install/bin/missions/SBG folder.
  It assumes that the template file sbg_vlidort_run.j is available in the folder.

  The script asks the user to provide:
    - an experiment name (exp_name)
    - the group id (group_id), i.e., the NCCS sponsor code to be used in SLURM.
 
  It will then create an experiment directory that has a self-contained and ready
  to use SLURM script sbg_vlidort_run.j.
"""

from pathlib import Path
import sys
import os
import subprocess
import shutil


def print_message():
    mssg = """
    ---------------------------------------------------------------------------------
    This setup script creates a self-contained experiment directory to run
    sampler.

    The script is interactive and asks the user to provide:

        - an experiment name (exp_name)
        - an experiment directory (experiment_directory)
        - the group id (group_id), i.e., the NCCS sponsor code to be used in SLURM.
 
    It will then create an experiment directory that has a self-contained and ready
    to use SLURM script sbg_vlidort_run.j.
    ---------------------------------------------------------------------------------
    """
    print(mssg)

def search_reaplace_in_file(loc_filename: Path, 
                            target_dir: Path, 
                            dict_words: dict) -> None:
    """
    Take a file template to search and replace collection of words.
    The new file (with the same name) will be created in the target directory.

    Parameters
    ----------
    loc_filename : str
       Local template file name.
    target_dir : Path
       Target directory where the new file will be created.
    dict_words : dict
       Dictionary where the keys are old words and the corresponding values
       are the new words.
    """
    new_filename = target_dir / loc_filename.name
    shutil.copy(loc_filename, new_filename)

    try:
        with open(new_filename, 'r') as fid:
            file_content = fid.read()

        for key in dict_words:
            file_content = file_content.replace(key, dict_words[key])
            print(f"Successfully replaced '{key}' with '{dict_words[key]}' in '{new_filename}'.")

        with open(new_filename, 'w') as file:
            file.write(file_content)

    except FileNotFoundError:
        print(f"Error: File '{new_filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def create_experiment_directory():

    # Get the current directory
    # Will be in the form FULL_PATH/AeroApps/install/bin/missions/SBG
    current_directory = Path.cwd()

    # Determine the source code main directory
    # Will be FULL_PATH/AeroApps
    source_directory = current_directory.parent.parent.parent.parent 

    reference_directory = source_directory.parent

    # Obtain the experiment name
    experiment_name = input("Provide the experiment name (one word):  ")
    experiment_name = experiment_name.strip()

    if not experiment_name:
        print("You need to provide and experiment name")
        sys.exit()

    if len(experiment_name.split()) > 1:
        print(f"The experiment name ({experiment_name}) should be in one word.")
        sys.exit()

    # Create the experiment directory
    experiment_directory = reference_directory / experiment_name
    print()
    exp_dir = input(f"Provide the experiment directory you want to use (default: {experiment_directory}): ")
    exp_dir = exp_dir.strip()

    if not exp_dir:
        pass
    else:
        experiment_directory = exp_dir

    print(f"The following experiment directory will be created: \n\n\t {experiment_directory}")
    print()

    experiment_directory.mkdir(parents=True, exist_ok=True)

    # Copy the config yaml files to the experiment directory
    configs = ["sbg_vlidort.yaml","sbg.yaml","ssd650.yaml","sbg_files.yaml","m2_aop.yaml"]
    for cf in configs:
        config_filepath = current_directory / cf
        shutil.copy(config_filepath, experiment_directory / config_filepath.name)

    # Do you want the benchmark codes
    benchmark = input("Is this a benchmark experiment (y/n,default: n):  ")
    benchmark = benchmark.strip()

    if not benchmark:
        print("Defaulting to not providing benchmarking scripts")
        benchmark = 'n'

    if len(benchmark) > 1:
        print(f"The benchmark response should be y or n")
        sys.exit()

    if benchmark == 'y':
        # Copy script to the experiment directory.
        scripts = ["benchmarking/sbg_vlidort_pyexample.py","benchmarking/lamb_vlidort_pyexample.py",
                   "benchmarking/sbg_files.yaml","benchmarking/lamb_files.yaml","benchmarking/m2_aop.yaml"]
    for sc in scripts:
        config_filepath = current_directory / sc
        shutil.copy(config_filepath, experiment_directory / config_filepath.name)

    # Get the sponsor code id

    result = subprocess.run(["groups"], shell=True, capture_output=True, text=True)

    groups = result.stdout.strip().split()

    print(f"The list of available group ids is: \n\n\t {result.stdout.strip()}")
    print()
    my_group = input(f"Provide the group id do you want to use (default: {groups[0]}): ")
    my_group = my_group.strip()

    if not my_group:
        my_group = groups[0]

    print()
    print(f"Your group is is: {my_group}")
    print()

    if my_group not in groups:
        print(f"You selected and invalid group.")
        print(f"The group {groups[0]} will be used.")
        print(f"You can change the group id in the SLURM script available in the experiment directory")
        print()

    if benchmark == 'y':
        loc_filename = Path("benchmarking/sbg_vlidort_run.j")
    else:
        loc_filename = Path("sbg_vlidort_run.j")
    target_dir = experiment_directory
    dict_words = {"@SRCDIR": str(source_directory), "@GROUPID": my_group}
    search_reaplace_in_file(loc_filename, target_dir, dict_words)

    if benchmark == 'y':
        loc_filename = Path("benchmarking/lamb_vlidort_run.j")
        target_dir = experiment_directory
        dict_words = {"@SRCDIR": str(source_directory), "@GROUPID": my_group}
        search_reaplace_in_file(loc_filename, target_dir, dict_words)


    print()
    print("-"*70)
    print(f"The experiment directory was created: \n\n\t {experiment_directory}")
    print()
    print(f"Go to the folder and if necessary edit the file {loc_filename}.")
    print()
    print("From the experiment directory, issue the command: ")
    print(f"   sbatch {loc_filename}")
    print("-"*70)
    print()

if __name__ == "__main__":
    print_message()
    create_experiment_directory()

