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

*  Please install Miniconda following the instructions given in the link: [https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)
  *  IMPORTANT: Check the box to add Miniconda to your PATH environment variable during the installation. Or you can add it to your environment manually. You can follow the instructions in this link:  [How to add Miniconda to Path](https://www.quora.com/How-can-I-add-conda-command-into-the-PATH-environment-variable-so-that-it-recognizes-the-instruction-given-that-the-executable-program-is-already-installed)
        *  You can also find the environment variables by searching in the search box, if it is not in the same place with the instructions.
  *  Please note that, at the end, you might need to restart Command Prompt to see the changes.


*  Open the search field of Windows and enter CMD. Click on Command Prompt.

*  Check if Miniconda is installed. On the command prompt enter the following. If you see a list of Python packages, the Miniconda is installed. Otherwise, you might see an error (e.g. command not found).

  ``` bash
    conda list
  ```

*  Please note that you may need to restart Command Prompt to see that Miniconda works.

2. ##### Setup the project

*  Open the search field of Windows and enter CMD. Click on Command Prompt, if it is not already open.

*  Navigate to the folder you want to download the project. For example, I collect my projects under the folder name 'workspaces':
   *  TIP: You can easily copy the path of your folder by going to the folder in windows explorer and clicking on the filepath and pressing Ctrl + c. Then you can paste it in              Command Prompt with Ctrl + v.
   ![image](https://user-images.githubusercontent.com/64579032/147078737-2098241e-9588-4a2c-b5ff-9a0a75e6164c.png)

  ``` bash
    cd  <path-to-some-folder>
    # In my case the path would look like this:
    cd C:\Users\username\Documents\workspaces
  ```

*  Check if git is installed. If the following command is giving the version number git is already installed please skip the next step.  
Or it is a warning like 'git' is not recognized as an internal or external command'. git is not installed please continue with the next step.
  ``` bash
    git --version
  ```

*  If git is not already installed, you can install it by following the instructions in this link: [https://phoenixnap.com/kb/how-to-install-git-windows](https://phoenixnap.com/kb/how-to-install-git-windows).
  *  Please note that, at the end, you might need to restart Command Prompt to see the changes.
  *  Do not forget to navigate to the folder you want to download the project.


*  Clone the git repository of the project by entering the below commands into the Command Prompt. You might be asked to enter your Github credentials.

  ``` bash
    git clone <GIT_REPO_URL>
    # In our case our repository url is:
    git clone https://github.com/asreview/paper-kinderformularium.git    
  ```

*  Navigate to the project folder  

  ``` bash
    cd  <path-to-project-folder>  
    # In our case our project folder name will be "paper-kinderformularium"
    cd C:\Users\username\Documents\workspaces\paper-kinderformularium  
  ```

*  Create a virtual environment and activate it to prevent system interference.

  ``` bash
    conda create -n .env
    conda activate .env
  ```

*  Install the required Python libraries:

  ``` bash
    conda config --add channels conda-forge
    conda install --file installation/requirements-windows.txt
    conda install -c conda-forge jupyter notebook
  ```

  *  IMPORTANT: Please note that, during installation you might need to press yes (y) to confirm the installation.

### For Linux and Mac OS

1. ##### Install at least Python 3.8 if you do not have already.  

*  Open Your Terminal. In Mac OS open your Launchpad (in your Dock) and search Terminal. In Linux open search area and find Terminal.

*   Enter the following command in Terminal to check the Python version. If this returns a number higher than 3.8, Python is installed and you can continue with step two. If not, continue to the next step.  

  ``` bash
     python --version
  ```

*  Install Python. You can download an installer for Linux/Mac OS from [https://www.python.org/downloads/](https://www.python.org/downloads/).

2. ##### Setup the project

*  On your terminal, navigate to the folder you want to download the project.

  ``` bash
    cd  <path-to-some-folder>
  ```

*  Clone the git repository.

  ``` bash
    git clone <GIT_REPO_URL>
  ```

*  Navigate to the project folder

  ``` bash
    cd  <path-to-project-folder>
  ```

*  Create a virtual environment to prevent system interference and activate it.

  ``` bash
    pip install virtualenv
    python -m venv .env
    source /bin/activate
  ```

*  Install the required Python libraries

  ``` bash
    pip install -r installation/requirements.txt
  ```

## The scripts

The scripts, which are inside the `src` folder, are configured to process certain docx formatted files to extract detailed reference information as a CSV file.

#### Where to put the data to be processed

To be able to process docx files you need to create a folder named 'data' under the project folder. Then you should put the docx files to be processed under data folder. Scripts will create a `csv` folder and under the `data` folder. You can find the corresponding CSV outputs under the `csv` folder.

#### How to run the scripts

Once you installed You can do these steps to run the Scripts.

* Open Command Prompt or Terminal

*  Navigate to the project folder from your terminal.
  ``` bash
    cd  <path-to-project-folder>
  ```

*  Activate your environment.
  Type in Command Prompt/Terminal:

  ``` bash
    # For Linux and MacOS
    source /bin/activate
    # For Windows
    conda activate .env
  ```
*  Running Jupyter notebook:

  ``` bash
    # For all Operating Systems:
    jupyter notebook  
  ```

  *  TIP: If you are new to the Jupyter Notebook you can learn the basics from this link:
      [https://www.dataquest.io/blog/jupyter-notebook-tutorial/](https://www.dataquest.io/blog/jupyter-notebook-tutorial/)    
      You can run the cells in Jupyter notebook by pressing `CTRL^Enter`, `SHIFT^Enter` or you can use the `Play` button at the top.  


#### Input formats

- Input must be a .docx formatted file
- Input can contain the reference information in 3 different formats:
  - V1: Files contain the references as an Endnote
  - V2: Files contain the references as a table
  - V3: Files contain the references both as an Endnote and a table
- V1 files should have the headlines "references" or "referenties"
- V2 files should have the references given inside a table
- V2 file tables should have the column name "effect"
- V2 file effect column should have the information formatted as:
  - The paper title should be given bold
  - The author names should follow the title in a normal text format
  - The file should have "summary" or "samenvatting" underlined, right after the author names.
  As you can see it in the following image:

 ![Example_table_cell](https://user-images.githubusercontent.com/5430558/148247563-ea808cad-2fcc-48f9-a570-68ae8af37b72.png)

- If you would like a combined output from multiple RBA files, you can achieve this by naming the files with the following notation:  
  `[identical_number][letter][.][space][identical_document_name]`    
For example:
1a.some_document_name.docx and 1b.some_document_name.docx outputs 1.some_document_name.csv

#### Output format
- CSV formatted file
- Contains the columns 'pubmed_id', 'title', 'abstract', 'doi', 'final_included'

## License
The content in this repository is published under the [MIT license](https://github.com/asreview/paper-kinderformularium/blob/main/LICENSE).

## Contact
For any questions or remarks, please send an email to asreview@uu.nl.
