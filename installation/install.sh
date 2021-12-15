#!/bin/bash

fatalerror () {
    echo "===================== FATAL ERROR ======================" >&2
    echo "An error occurred during installation!!" >&2
    echo "$1" >&2
    echo "========================================================" >&2
    exit 2
}

usage () {
  echo
  echo "=========================== Manual ============================"
  echo
  echo "Start installation:"
  echo
  echo "    bash install.sh or zsh install.sh"
  echo
  echo "Run installation with certain python version. python3.8 is an example here:"
  echo
  echo "    bash install.sh -p python3.8"
  echo
  exit
}

if [[ $1 == "-h" ]]; then
  usage
fi

PYTHON=$(which python3)
PROJECT_PATH=$(cd ../; pwd) || fatalerror "Unable to find project folder path"
VENV_FILE='.env'

echo
echo "--------------------------------------------------------"
echo "Creating the virtual environment named as '$VENV_FILE'"
echo "--------------------------------------------------------"

$PYTHON -m venv $PROJECT_PATH/$VENV_FILE || fatalerror "Unable to create virtual environment"

echo
echo "--------------------------------------------------------"
echo "Activating the virtual environment"
echo "--------------------------------------------------------"

echo source $PROJECT_PATH/$VENV_FILE/bin/activate
source $PROJECT_PATH/$VENV_FILE/bin/activate
# pip install --upgrade pip
$PROJECT_PATH/$VENV_FILE/bin/python3 -m pip install --upgrade pip

echo
echo "--------------------------------------------------------"
echo "Installing the required Python packages through pip"
echo "--------------------------------------------------------"
# pip install -r requirements.txt || fatalerror "Unable to install Python requirements."

echo
echo "--------------------------------------------------------"
echo "== Installation completed =="
echo "Do not forget activating the virtual environment before running the project."
echo "You can activate the virtual environment as follows:  source $VIRTUAL_ENV/bin/activate"
echo "--------------------------------------------------------"
echo
