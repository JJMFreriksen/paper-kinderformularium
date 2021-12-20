# Scripts - Reviewing Scientific Literature for Risk-Benefit Analysis of the Dutch Pediatric Formulary 

This repository is part of a project of the Dutch Pediatric Formulary to evaluate
the benefit of using ASReview LAB for the development of Risk Benefit Analysis (RBA).
More specific, within this project the following is investigated:

1. What relevant information can be retrieved by using ASReview LAB for screening on top of
the results found within an existing Risk Benefit Analysis?
2. How many relevant papers can be identified using a broader search strategy and using ASReview LAB for screening, instead of the existing review strategy?

Within this repository one can find the scripts that are used to pre-process the
Risk Analysis Documents to ASReview LAB compatible input files.

## Installation

### For Windows

1. ##### Install Miniconda (or Anaconda) for Python 3.8 or higher  

  *   Please install Miniconda following the instructions given in the link: [https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)

  * Open the search field of Windows and enter CMD. Click on Command Prompt.

  *   Check if Miniconda is installed. On the command prompt enter the following. If you see a list of Python packages, the Miniconda is installed. Otherwise, you might see an error (e.g. command not found).

      ``` bash
      $ conda list
      ```

      Please note that you may need to restart Command Prompt to see that Miniconda works.

2. ##### Setup the project

   *   Open the search field of Windows and enter CMD. Click on Command Prompt, if it is not already open.

   *   Navigate to the folder you want to download the project. For example, I collect my projects under the folder name 'workspaces':

       ``` bash
       $ cd  <path-to-some-folder>
       # In my case the path would look like this:
       $ cd C:\Users\username\Documents\workspaces
       ```

   *   If git is not already installed, you can install it by following the instructions in this link: [https://phoenixnap.com/kb/how-to-install-git-windows](https://phoenixnap.com/kb/how-to-install-git-windows).

   *   Clone the git repository of the project by entering the below commands into the Command Prompt. You might be asked to enter your Github credentials.

       ``` bash
       $ git clone <GIT_REPO_URL>
       # In our case our repository url is:
       $ git clone https://github.com/asreview/paper-kinderformularium.git    
       ```

   *   Navigate to the project folder  

       ``` bash
       $ cd  <path-to-project-folder>  
       # In our case our project folder name will be "paper-kinderformularium"
       $ cd C:\Users\username\Documents\workspaces\paper-kinderformularium  
       ```

   *   Create a virtual environment and activate it to prevent system interference.

       ``` bash
       $ conda create -n .env
       $ conda activate .env
       ```

   *   Install the required Python libraries:

        ``` bash
        $ conda config --add channels conda-forge
        $ conda install --file installation/requirements-windows.txt
        $ conda install -c conda-forge jupyter notebook
        ```

### For Linux and Mac OS

1. ##### Install at least Python 3.8 if you do not have already.  

    *   Open Your Terminal. In Mac OS open your Launchpad (in your Dock) and search Terminal. In Linux open search area and and find Terminal.

    *   Enter the following command in Terminal to check the Python version. If this returns a number higher than 3.8, Python is installed and you can continue with step two. If not, continue to the next step.  

        ``` bash
        $ python --version
        ```

    *   Install Python. You can download an installer for Linux/Mac OS from [https://www.python.org/downloads/](https://www.python.org/downloads/).

2. ##### Setup the project

    *   On your terminal, navigate to the folder you want to download the project.

        ``` bash
        $ cd  <path-to-some-folder>
        ```

    *   Clone the git repository.

        ``` bash
        $ git clone <GIT_REPO_URL>
        ```

    *   Navigate to the project folder

        ``` bash
        $ cd  <path-to-project-folder>
        ```

    *   Create a virtual environment to prevent system interference and activate it.

        ``` bash
        $ pip install virtualenv
        $ python -m venv .env
        $ source /bin/activate
        ```

    *   Install the required Python libraries

        ``` bash
        $ pip install -r installation/requirements.txt
        ```

## The scripts

The scripts, which are inside the `src` folder, are configured to process certain docx formatted files to extract detailed reference information as a CSV file.

#### Running Jupyter notebook  

After activating the virtual environment (On Windows with `conda` command) type in Command Prompt/Terminal:

  * For windows:

        ``` bash
        $ jupyter-lab
        ```

  * For Linux and MacOS:

        ``` bash
        $ jupyter notebook  
        ```

You can run the cells in Jupyter notebook by pressing `CTRL^Enter`.

#### Input formats
- Input must be a .docx formatted file
- Input can contain the reference information in 3 different formats:
  - V1: Files contain the references as an Endnote
  - V2: Files contain the references as a table
  - V3: Files contain the references both as an Endnote and a table
- V1 files should have the headlines "references" or "referenties"
- V2 files should have the references given inside a table
- V2 file tables should has the column name "effect"
- V2 file effect column should have the information formatted as:
  - The paper title should be given bold
  - The author names should follow the title in a normal text format
  - The file should have "summary" or "samenvatting" underlined, right after the author names

#### Output format
- CSV formatted file
- Contains the columns 'pubmed_id', 'title', 'abstract', 'doi', 'final_included'

## License
The content in this repository is published under the [MIT license](https://github.com/asreview/paper-kinderformularium/blob/main/LICENSE).

## Contact
For any questions or remarks, please send an email to asreview@uu.nl.
