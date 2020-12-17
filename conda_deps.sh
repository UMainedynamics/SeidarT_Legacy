#!/bin/bash

SEIDART_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" # current directory
ARCH=$(uname -m)    # machine hardware
KERN=$(uname -s)    # kernel name
NODE=$(uname -n)    # node name
TMP="/tmp"          # temporary directory for install file
EXE="conda-install.sh" # install file name
TMP_EXE="$TMP/$EXE" # install file loc/name
MACOS_EXE="Miniconda3-latest-MacOSX-x86_64.sh"
LINUX_EXE="Miniconda3-latest-Linux-x86_64.sh"
ARM_EXE="Berryconda3-2.0.0-Linux-armv7l.sh"
X86_BASE_URL="https://repo.anaconda.com/miniconda/"
ARM_BASE_URL="https://github.com/jjhelmus/berryconda/releases/download/v2.0.0/"
if [[ "$ARCH" == "armv"* ]]; then RELEASE='berryconda3'; else RELEASE='miniconda3'; fi
# conda install location:
PREFIX="$HOME/$RELEASE"         # $HOME/miniconda3 is default location
FULL="$HOME/anaconda3"          # full release install location
BERRYCONDA="$HOME/berryconda3"  # berryconda install location
MINICONDA="$HOME/miniconda3"    # miniconda install location

echo "Please follow instructions in prompts."
read -rp $'Press Enter to continue...\n'

# first we have to test if there is an existing anaconda installation
# the simplest case, that the conda command works:
echo "Looking for conda installation in $HOME..."
command -v conda >/dev/null 2>&1 &&
conda activate >/dev/null 2>&1 &&
CONDA_EXISTS=1

if [ -z ${CONDA_EXISTS+x} ]; then
  # if conda command doesn't exist,
  if [ -f "$MINICONDA/bin/conda" ]; then
    # now we look in the default install location
    . $PREFIX/etc/profile.d/conda.sh &&
    conda activate &&
    CONDA_EXISTS=1
  elif [ -f "$BERRYCONDA/bin/conda" ]; then
    # look for a berryconda release
    . $BERRYCONDA/etc/profile.d/conda.sh &&
    conda activate &&
    CONDA_EXISTS=1
  elif [ -f "$FULL/bin/conda" ]; then
    # finally, look for a full release
    . $FULL/etc/profile.d/conda.sh &&
    PREFIX=$FULL &&
    conda activate &&
    CONDA_EXISTS=1
  elif grep -Fxq "etc/profile.d/conda.sh" "$HOME/.bashrc"; then
    `grep etc/profile.d/conda.sh $HOME/.bashrc` &&
    conda activate &&
    CONDA_EXISTS=1
  else
    conda="$PREFIX/bin/conda"
  fi
else
  PREFIX="$(cd $(dirname $(which conda))/../; pwd)"
fi

if [ -z ${CONDA_EXISTS+x} ]; then
  echo "Cannot find conda installation; will try installing $RELEASE."
  # get ready to install anaconda or berryconda
  echo "Found $KERN environment on $ARCH."
  echo "Install location: $PREFIX"
  echo "Ready to download $RELEASE"
  echo "The download could be as large as 90 MB."
  read -n1 -rp $'Press any key to continue or Ctrl+C to exit...\n\n'

  if [ ! -z ${PYTHONPATH+x} ]; then
    # conda does not like $PYTHONPATH, and $PYTHONPATH is deprecated,
    # so we can get away with disabling it during installation.
    # because it is sourced, it will come back when the user opens a new shell
    # and conda will complain about it directly to the user.
    unset $PYTHONPATH
  fi

  if [[ "$ARCH" == "armv"* ]]; then
    # installing on ARM architecture (RPi or similar)
    rpi="rpi"

    if [[ "$NODE" == "raspberryshake" ]]; then
      # warn the user about installing on a Shake
      echo '---------------------------------------------------------------'
      echo "WARNING: You are installing this on the Raspberry Shake itself."
      echo "Although this is possible, it is not tested or supported."
      echo "Raspberry Shake S.A. is not liable for strange Shake behavior"
      echo "if you choose to do this! Proceed at your own risk."
      read -n1 -rp $'Press any key to continue or Ctrl+C to exit...\n'
    fi

    wget "$ARM_BASE_URL$ARM_EXE" -O "$TMP_EXE" && DL=1

  else
    if [[ "$KERN" == "Linux" ]]; then
      CONDA_INSTALLER=$LINUX_EXE
      wget "$X86_BASE_URL$CONDA_INSTALLER" -O "$TMP_EXE" && DL=1

    elif [[ "$KERN" == "Darwin" ]]; then
      CONDA_INSTALLER=$MACOS_EXE
      curl "$X86_BASE_URL$CONDA_INSTALLER" -o "$TMP_EXE" && DL=1

    else
      echo "ERROR: Script does not support this OS."
      echo "Please install Anaconda 3 by hand from the following link:"
      echo "https://www.anaconda.com/distribution/#download-section"
      exit 1
    fi

  fi

  if [ ! -z ${DL+x} ]; then
    chmod +x "$TMP_EXE"
    echo "Installing $RELEASE..."
    cd "$TMP" && ./$EXE -b -p $PREFIX
    echo "Cleaning up temporary files..."
    rm "$TMP_EXE"
    echo "Updating base conda environment..."
    $conda update conda -y
  else
    echo "Something went wrong downloading $RELEASE. Check the error and try again."
    exit 2
  fi

else
    PREVIOUS_CONDA=1
    echo "Anaconda installation found at $PREFIX"
    echo "conda executable: $(which conda)"
fi

COMMENT="# added by SeidarT/conda installer"
SOURCELINE=". $PREFIX/etc/profile.d/conda.sh"

if grep -Fxq "$SOURCELINE" "$HOME/.bashrc"; then
  echo "Source line already exists in $HOME/.bashrc"
  SOURCED=1
else
  echo "----------------------------------------------"
  echo "The script will now append a sourcing line to your ~/.bashrc file in order to"
  echo 'make activating conda easier (just type "conda activate" into a terminal).'
  echo 'This is highly recommended to make it easier to access conda environments in the future.'
  echo "This line is: $SOURCELINE"
  read -rp $'Press Enter to continue, or type no and press Enter to prevent this...\n' key

  if [[ "$key" == "no" ]]; then
    echo "Not appending sourcing line to bashrc."
    echo "You can add it later by adding the following line to the bottom of ~/.bashrc:"
    echo $SOURCELINE
  else
    echo "Appending sourcing line to bashrc..."
    echo $COMMENT >> $HOME/.bashrc
    echo $SOURCELINE >> $HOME/.bashrc
    SOURCED=1
  fi
fi
echo "Sourcing..."
$SOURCELINE
echo "Activating conda..."
conda activate && CONDA_EXISTS=1

if [ -z ${CONDA_EXISTS+x} ]; then
  echo "ERROR: Anaconda install failed. Check the error output and try again."
  exit 2
fi

if [[ "$ARCH" == "armv"* ]]; then
  # may not need to differentiate depending on testing, otherwise arm python=3.6.6
  ENV_INSTALL="conda create -n SeidarT python=3 pip git dos2unix ghostscript imagemagick numpy matplotlib scipy pyevtk vtk -y"
else
  ENV_INSTALL="conda create -n SeidarT python=3 pip git dos2unix ghostscript imagemagick numpy matplotlib scipy pyevtk vtk -y"
fi

# check for conda forge channel; if it's not there add it
if [ ! -f $HOME/.condarc ]; then
  echo "No $HOME/.condarc file exists. Creating..."
  echo $'channels:\n  -\n   defaults\n  -\n   rpi\n  -\n   conda-forge\n' > $HOME/.condarc
fi
if [[ "$ARCH" == "armv"* ]]; then
  cat $HOME/.condarc | grep "rpi" >/dev/null && echo "Found rpi channel in $HOME/.condarc" ||
  (echo "Appending rpi to conda channels..." &&
  conda config --append channels rpi)
fi
cat $HOME/.condarc | grep "conda-forge" >/dev/null && echo "Found conda-forge channel in $HOME/.condarc"  ||
(echo "Appending conda-forge to conda channels..." &&
conda config --append channels conda-forge)

if [ -d $PREFIX/envs/SeidarT ]; then
  echo "Another SeidarT conda environment already exists at $PREFIX/envs/SeidarT" &&
  echo "Do you want to use it, or remove it and install a new one?"
  read -rp $'Press Enter to use it, or type yes and press Enter to reinstall:\n' REINSTALL
  
  if [[ "$REINSTALL" == "yes" ]]; then
    echo "Removing old environment..."
    rm -r $PREFIX/envs/SeidarT
    echo "Reinstalling SeidarT conda environment..." &&
    $ENV_INSTALL
  fi
else
  echo "Creating and installing SeidarT conda environment..." &&
  $ENV_INSTALL
fi

if [ -d $PREFIX/envs/SeidarT ]; then
  conda env config vars set PATH="$PATH:$SEIDART_DIR/bin" -n SeidarT
  echo "Activating SeidarT environment..." &&
  conda activate SeidarT && echo "Success: SeidarT environment activated." &&
  echo "Installing SeidarT..." &&
  pip install mplstereonet && SUCCESS=1
else
  echo "ERROR: SeidarT failed to install."
fi

if [ ! -z ${SUCCESS+x} ]; then
  echo "SeidarT has installed successfully!"

  # from https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
  echo "Appending $SEIDART_DIR/bin to PATH in SeidarT environment..."
  mkdir -p $CONDA_PREFIX/etc/conda/activate.d
  mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
  touch $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
  touch $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
  echo 'export PATH="$PATH:'"$SEIDART_DIR/bin"'"' >> "$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh"
  echo "unset PATH" >> "$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh"

  if [ -z ${PREVIOUS_CONDA+x} ]; then
    if [ -z ${SOURCED+x} ]; then
      echo 'You will need to tell your shell where to find conda by entering ". ~/'"$RELEASE"'/etc/profile.d/conda.sh"'
      THEN='then '
    fi
    echo 'You can'$THEN' enter the command "conda activate SeidarT" to activate the SeidarT conda environment'
  else
    if [ -z ${SOURCED+x} ]; then
      echo 'You need to re-source your shell before using conda. To do this, type "source ~/.bashrc" or just open a new terminal window.'
      THEN='then '
    fi
    echo 'You can'$THEN' enter the SeidarT conda environment by typing "conda activate SeidarT"'
  fi
  echo 'and then run any of the SeidarT scripts as stated in the manual'
  exit 0
else
  echo "---------------------------------"
  echo "Something went wrong."
  echo "Check the error output and try again."
  exit 2
fi
