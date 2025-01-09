#!/bin/bash

# Name of the environment
ENV_NAME="phgtools"

# Download the .yml file from the PHG repository at the same folder as the script is located
wget https://github.com/maize-genetics/phg_v2/raw/refs/heads/main/src/main/resources/phg_environment.yml

echo "The .yml file has been downloaded"

# Change the name of the file
mv phg_environment.yml phgtools_environment.yml
sed -i 's/name: .*/name: phgtools/' phgtools_environment.yml

echo "The .yml file has been renamed"

# Add the pygenometracks package to the .yml file
echo "  - pandas" >> phgtools_environment.yml
echo "  - pip" >> phgtools_environment.yml
echo "  - pip:" >> phgtools_environment.yml
echo "      - pygenometracks" >> phgtools_environment.yml
echo "      - jupyter" >> phgtools_environment.yml

# Check if the environment already exists and remove it if it does
if conda env list | grep -q $ENV_NAME; then
    echo "The environment $ENV_NAME already exists. Removing it..."
    conda env remove -n $ENV_NAME
fi

# Install the environment
conda env create -f phgtools_environment.yml

# Initialize Conda for bash
conda init bash

# Reload the shell to apply changes
#exec $SHELL
#echo "The shell has been reloaded"
# Activate the environment
conda activate phgtools

echo "The environment has been set up and pygenometracks has been installed."
echo "The environment name is: $ENV_NAME"
echo "To activate the environment, run the following command: conda activate $ENV_NAME"

# Remove the .yml file
rm phgtools_environment.yml
echo "The .yml file has been removed"

# Check if the pwd is included in the PYTHONPATH
if [[ ":$PYTHONPATH:" != *":$PWD:"* ]]; then
    echo "export PYTHONPATH=$PWD:\$PYTHONPATH" >> ~/.bashrc
    echo "The current directory has been added to the PYTHONPATH"
fi
source ~/.bashrc

