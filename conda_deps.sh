#!/bin/bash

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" # current directory
arch=$(uname -m)    # machine hardware
os=$(uname -s)      # kernel name
node=$(uname -n)    # node name
gnu=$(uname -a | grep GNU)  # is this GNU?
tmp="/tmp"          # temporary directory for install file
exe="conda-install.sh" # install file name
tmp_exe="$tmp/$exe" # install file loc/name
conda="conda"       # anaconda executable or alias
macos_exe="Miniconda3-latest-MacOSX-x86_64.sh"
linux_exe="Miniconda3-latest-Linux-x86_64.sh"
arm_exe="Berryconda3-2.0.0-Linux-armv7l.sh"
x86_base_url="https://repo.anaconda.com/miniconda/"
arm_base_url="https://github.com/jjhelmus/berryconda/releases/download/v2.0.0/"
if [[ "$arch" == "armv"* ]]; then release='berryconda3'; else release='miniconda3'; fi
# conda install location:
prefix="$HOME/$release"         # $HOME/miniconda3 is default location
full="$HOME/anaconda3"          # full release install location
berryconda="$HOME/berryconda3"  # berryconda install location
miniconda="$HOME/miniconda3"    # miniconda install location

echo "Please follow instructions in prompts."
read -rp $'Press Enter to continue...\n'

# first we have to test if there is an existing anaconda installation
# the simplest case, that the conda command works:
echo "Looking for conda installation in $HOME..."
command -v conda >/dev/null 2>&1 &&
conda activate >/dev/null 2>&1 &&
conda_exists=1

if [ -z ${conda_exists+x} ]; then
  # if conda command doesn't exist,
  if [ -f "$miniconda/bin/conda" ]; then
    # now we look in the default install location
    . $prefix/etc/profile.d/conda.sh &&
    conda activate &&
    conda_exists=1
  elif [ -f "$berryconda/bin/conda" ]; then
    # look for a berryconda release
    . $berryconda/etc/profile.d/conda.sh &&
    conda activate &&
    conda_exists=1
  elif [ -f "$full/bin/conda" ]; then
    # finally, look for a full release
    . $full/etc/profile.d/conda.sh &&
    prefix=$full &&
    conda activate &&
    conda_exists=1
  elif grep -Fxq "etc/profile.d/conda.sh" "$HOME/.bashrc"; then
    `grep etc/profile.d/conda.sh $HOME/.bashrc` &&
    conda activate &&
    conda_exists=1
  else
    conda="$prefix/bin/conda"
  fi
else
  prefix="$(cd $(dirname $(which conda))/../; pwd)"
fi

if [ -z ${conda_exists+x} ]; then
  echo "Cannot find conda installation; will try installing $release."
  # get ready to install anaconda or berryconda
  echo "Found $os environment on $arch."
  echo "Install location: $prefix"
  echo "Ready to download $release"
  echo "The download could be as large as 90 MB."
  read -n1 -rp $'Press any key to continue or Ctrl+C to exit...\n\n'

  if [ ! -z ${PYTHONPATH+x} ]; then
    # conda does not like $PYTHONPATH, and $PYTHONPATH is deprecated,
    # so we can get away with disabling it during installation.
    # because it is sourced, it will come back when the user opens a new shell
    # and conda will complain about it directly to the user.
    unset $PYTHONPATH
  fi

  if [[ "$arch" == "armv"* ]]; then
    # installing on ARM architecture (RPi or similar)
    rpi="rpi"

    if [[ "$node" == "raspberryshake" ]]; then
      # warn the user about installing on a Shake
      echo '---------------------------------------------------------------'
      echo "WARNING: You are installing this on the Raspberry Shake itself."
      echo "Although this is possible, it is not tested or supported."
      echo "Raspberry Shake S.A. is not liable for strange Shake behavior"
      echo "if you choose to do this! Proceed at your own risk."
      read -n1 -rp $'Press any key to continue or Ctrl+C to exit...\n'
    fi

    wget "$arm_base_url$arm_exe" -O "$tmp_exe" && dl=1

  else
    if [[ "$os" == "Linux" ]]; then
      conda_installer=$linux_exe
      wget "$x86_base_url$conda_installer" -O "$tmp_exe" && dl=1

    elif [[ "$os" == "Darwin" ]]; then
      conda_installer=$macos_exe
      curl "$x86_base_url$conda_installer" -o "$tmp_exe" && dl=1

    else
      echo "ERROR: Script does not support this OS."
      echo "Please install Anaconda 3 by hand from the following link:"
      echo "https://www.anaconda.com/distribution/#download-section"
      exit 1
    fi

  fi

  if [ ! -z ${dl+x} ]; then
    chmod +x "$tmp_exe"
    echo "Installing $release..."
    cd "$tmp" && ./$exe -b -p $prefix
    echo "Cleaning up temporary files..."
    rm "$tmp_exe"
    echo "Updating base conda environment..."
    $conda update conda -y
  else
    echo "Something went wrong downloading $release. Check the error and try again."
    exit 2
  fi

else
    previous_conda=1
    echo "Anaconda installation found at $prefix"
    echo "conda executable: $(which conda)"
fi

comment="# added by SeidarT/conda installer"
sourceline=". $prefix/etc/profile.d/conda.sh"

if grep -Fxq "$sourceline" "$HOME/.bashrc"; then
  echo "Source line already exists in $HOME/.bashrc"
  sourced=1
else
  echo "----------------------------------------------"
  echo "The script will now append a sourcing line to your ~/.bashrc file in order to"
  echo 'make activating conda easier (just type "conda activate" into a terminal).'
  echo 'This is highly recommended to make it easier to access conda environments in the future.'
  echo "This line is: $sourceline"
  read -rp $'Press Enter to continue, or type no and press Enter to prevent this...\n' key

  if [[ "$key" == "no" ]]; then
    echo "Not appending sourcing line to bashrc."
    echo "You can add it later by adding the following line to the bottom of ~/.bashrc:"
    echo $sourceline
  else
    echo "Appending sourcing line to bashrc..."
    echo $comment >> $HOME/.bashrc
    echo $sourceline >> $HOME/.bashrc
    sourced=1
  fi
fi
echo "Sourcing..."
$sourceline
echo "Activating conda..."
conda activate && conda_exists=1

if [ -z ${conda_exists+x} ]; then
  echo "ERROR: Anaconda install failed. Check the error output and try again."
  exit 2
fi

if [[ "$arch" == "armv"* ]]; then
  # may not need to differentiate depending on testing, otherwise arm python=3.6.6
  env_install="conda create -n SeidarT python=3 pip git dos2unix ghostscript imagemagick numpy matplotlib scipy pyevtk vtk -y"
else
  env_install="conda create -n SeidarT python=3 pip git dos2unix ghostscript imagemagick numpy matplotlib scipy pyevtk vtk -y"
fi

# check for conda forge channel; if it's not there add it
if [ ! -f $HOME/.condarc ]; then
  echo "No $HOME/.condarc file exists. Creating..."
  echo $'channels:\n  -\n   defaults\n  -\n   rpi\n  -\n   conda-forge\n' > $HOME/.condarc
fi
if [[ "$arch" == "armv"* ]]; then
  cat $HOME/.condarc | grep "rpi" >/dev/null && echo "Found rpi channel in $HOME/.condarc" ||
  (echo "Appending rpi to conda channels..." &&
  conda config --append channels rpi)
fi
cat $HOME/.condarc | grep "conda-forge" >/dev/null && echo "Found conda-forge channel in $HOME/.condarc"  ||
(echo "Appending conda-forge to conda channels..." &&
conda config --append channels conda-forge)

if [ -d $prefix/envs/SeidarT ]; then
  echo "Another SeidarT conda environment already exists at $prefix/envs/SeidarT" &&
  echo "Do you want to use it, or remove it and install a new one?"
  read -rp $'Press Enter to use it, or type yes and press Enter to reinstall:\n' reinstall
  
  if [[ "$reinstall" == "yes" ]]; then
    echo "Removing old environment..."
    rm -r $prefix/envs/SeidarT
    echo "Reinstalling SeidarT conda environment..." &&
    $env_install
  fi
else
  echo "Creating and installing SeidarT conda environment..." &&
  $env_install
fi

if [ -d $prefix/envs/SeidarT ]; then
  echo "Activating SeidarT environment..." &&
  conda activate SeidarT && echo "Success: SeidarT environment activated." &&
  echo "Installing SeidarT..." &&
  pip install mplstereonet && success=1
else
  echo "ERROR: SeidarT failed to install."
fi

if [ ! -z ${success+x} ]; then
  echo "SeidarT has installed successfully!"

  if [ -z ${previous_conda+x} ]; then
    if [ -z ${sourced+x} ]; then
      echo 'You will need to tell your shell where to find conda by entering ". ~/'"$release"'/etc/profile.d/conda.sh"'
      then='then '
    fi
    echo 'You can'$then' enter the command "conda activate SeidarT" to activate the SeidarT conda environment'
  else
    if [ -z ${sourced+x} ]; then
      echo 'You need to re-source your shell before using conda. To do this, type "source ~/.bashrc" or just open a new terminal window.'
      then='then '
    fi
    echo 'You can'$then' enter the SeidarT conda environment by typing "conda activate SeidarT"'
  fi
  echo 'and then run any of the SeidarT scripts as stated in the manual'
  exit 0
else
  echo "---------------------------------"
  echo "Something went wrong."
  echo "Check the error output and try again."
  exit 2
fi
