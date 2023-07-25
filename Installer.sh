#!/bin/bash
PYTHON_VERSION="3"

echo "Welcome to the interactive setup script for the RNA editing analysis pipeline."
echo "This script will guide you through the setup process."
echo "Please make sure that you have the following dependencies installed:"
echo "  - Python $PYTHON_VERSION"
echo "  - virtualenv"
echo "  - ViennaRNA (tested with version 2.4.2)"
echo "  - UMI-tools (tested with version 1.1.0)"
echo "  - STAR aligner (tested with version 2.6.1d)"
echo "  - BWA aligner (tested with version 0.7.17)"
echo "  - SAMtools (tested with version 1.9)"
echo "  - Cutadapt (tested with version 2.10)"
echo "If you have not installed these dependencies, please do so before continuing."
echo "Press enter to continue."
read

echo "This tool will setup a Python virtual environment and install all required packages through PIP."
echo "Press enter to continue."
read

# Check if virtualenv is installed; if not, install it
if ! command -v virtualenv &>/dev/null; then
    echo "Installing virtualenv..."
    python$PYTHON_VERSION -m pip install virtualenv
fi

# Create a virtual environment using virtualenv
virtualenv -p python$PYTHON_VERSION venv

# Install required packages using pip
venv/bin/pip install -r biopython==1.79 cutadapt==4.0 cycler==0.11.0 Cython==3.0.0 dnaio==0.8.1 fonttools==4.32.0 isal==0.11.1 kiwisolver==1.4.2 matplotlib==3.5.1 numpy==1.21.6 packaging==21.3 Pillow==9.1.0 pyparsing==3.0.8 pysam==0.15.4 python-dateutil==2.8.2 setuptools==41.2.0 six==1.16.0 typing-extensions==4.2.0 xopen==1.5.0 

# Function to prompt the user with path completion support
prompt_with_completion() {
    local prompt_message="$1"
    local user_input

    # Read user input with readline's -e option to enable path completion
    read -e -p "$prompt_message" -i "$user_input" user_input

    # Output the user's selected option
    echo "$user_input"
}

# Function to check if a path exists
check_path_exists() {
    local path="$1"

    if [! -e "$path" ]; then
        echo "File does not exist: $path"
    fi
}

path_to_ViennaRNA=$(prompt_with_completion "Enter path to ViennaRNA executable: ")
check_path_exists "$path_to_ViennaRNA"
path_to_umi_tools=$(prompt_with_completion "Enter path to UMI Tools executable: ")
check_path_exists "$path_to_umi_tools"
path_to_star=$(prompt_with_completion "Enter path to STAR  aligner executable: ")
check_path_exists "$path_to_star"
path_to_bwa=$(prompt_with_completion "Enter path to BWA aligner executable: ")
check_path_exists "$path_to_bwa"
path_to_samtools=$(prompt_with_completion "Enter path to SAMtools executable: ")
check_path_exists "$path_to_samtools"
path_to_cutadapt=$(prompt_with_completion "Enter path to Cutadapt executable: ")
check_path_exists "$path_to_cutadapt"

# Create a config file
echo "Creating config file..."
echo "path_to_ViennaRNA = \"$path_to_ViennaRNA\"" > config.txt
echo "path_to_umi_tools = \"$path_to_umi_tools\"" >> config.txt
echo "path_to_star = \"$path_to_star\"" >> config.txt
echo "path_to_bwa = \"$path_to_bwa\"" >> config.txt
echo "path_to_samtools = \"$path_to_samtools\"" >> config.txt
echo "path_to_cutadapt = \"$path_to_cutadapt\"" >> config.txt

# Indicate that the setup is complete
echo "Setup process complete."
