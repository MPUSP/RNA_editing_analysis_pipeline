#!/bin/bash
PYTHON_VERSION="3"

echo "Welcome to the interactive setup script for the RNA editing analysis pipeline."
echo "This script will guide you through the setup process."
echo "Please make sure that you have the following dependencies installed:"
echo "  - Python $PYTHON_VERSION"
echo "  - virtualenv"
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
venv/bin/pip install biopython==1.79 cycler==0.11.0 Cython==3.0.0 dnaio==0.8.1 fonttools==4.32.0 isal==0.11.1 kiwisolver==1.4.2 matplotlib numpy packaging==21.3 Pillow==9.1.0 pyparsing==3.0.8 python-dateutil==2.8.2 setuptools==41.2.0 six==1.16.0 typing-extensions==4.2.0 xopen==1.5.0 
venv/bin/pip install pysam==0.21.0

# Function to prompt the user with path completion support
prompt_with_completion() {
    local prompt_message="$1"
    local user_input=""

    while true; do
        # Read user input with readline's -e option to enable path completion
        read -e -p "$prompt_message" user_input 

        # If the user input is empty, prompt again
        if [ -z "$user_input" ]; then
            echo "Please enter a value." >&2
        # check if the path exists
        elif [ ! -e "$user_input" ]; then
            echo "Path does not exist: $user_input">&2
        else
            break
        fi
    done
    
    # Output the user's selected option
    echo "$user_input" "-> path confirmed." >&2
    echo "$user_input"
}

path_to_umi_tools=$(prompt_with_completion "Enter path to UMI Tools executable: ")
path_to_star=$(prompt_with_completion "Enter path to STAR  aligner executable: ")
path_to_bwa=$(prompt_with_completion "Enter path to BWA aligner executable: ")
path_to_samtools=$(prompt_with_completion "Enter path to SAMtools executable: ")
path_to_cutadapt=$(prompt_with_completion "Enter path to Cutadapt executable: ")

# Create a config file
echo "Creating config file..."
echo "path_to_umi_tools = \"$path_to_umi_tools\"" > config.txt
echo "path_to_star = \"$path_to_star\"" >> config.txt
echo "path_to_bwa = \"$path_to_bwa\"" >> config.txt
echo "path_to_samtools = \"$path_to_samtools\"" >> config.txt
echo "path_to_cutadapt = \"$path_to_cutadapt\"" >> config.txt

# Indicate that the setup is complete
echo "Setup process complete."
