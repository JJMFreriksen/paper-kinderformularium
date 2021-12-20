# Pre-processing systematic review data for the Dutch Pediatric Formulary

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
<!-- 1. ##### Prerequisites   -->
1. Please install at least Python 3.8 if you do not have already.  
    a. Open the search field of Windows and enter CMD.  
    b. Click on Command Prompt.   
    c. Write python --version in the Command Prompt and press ENTER.  
    d. If this returns a number higher than 3.8, Python is installed and you can continue with step two. If not, continue to the next step.  
    e. Install Python.
    You can download an installer for Windows from <https://www.python.org/downloads/>  
    f. Install Python for Windows. IMPORTANT: Check the box to add Python to your PATH environment variable during the installation.
    Or Add Python to your environment manually. You can follow the instructions in this link: https://datatofish.com/add-python-to-windows-path/  
    ** Notes  
      * You can also find the environment variables by searching in the search box. If it is not in the same place with the instructions.
      * Please note that at the end you might need to restart Command Prompt to see the changes.

  g. Please install Miniconda with the instructions given in this link: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html  

  h. Check if Miniconda is installed.  
  * Please note that at the end you might need to restart Command Prompt.

    $ conda list

2. #### Setup the project
   a. Navigate to the folder you want to work on. For example, I collect my projects under the folder name 'workspaces'.

       $ cd C:\Users\username\Documents\workspaces

   b. Clone the git repository. If git is not installed, you can install it by following the instructions in this link: https://phoenixnap.com/kb/how-to-install-git-windows. You might be asked to enter your Github credentials.

       $ git clone <GIT_REPO_URL>
       # In our case our repository url is:
       $ git clone https://github.com/asreview/paper-kinderformularium.git    

  c. Navigate to the project folder  

       $ cd  <path-to-project-folder>  
       # In our case our project folder is:  
       $ cd C:\Users\username\Documents\workspaces\paper-kinderformularium  

   d. Create a virtual environment to prevent system interference and activate it.

         $ conda create -n .env
         $ conda activate .env

   d. Install the required Python libraries

        $ conda install --file installation/requirements-w.txt

### For Linux and Mac OS

1. Please install at least Python 3.8 if you do not have already.  
    a. Open Your Terminal.In Mac OS open your Launchpad (in your Dock) and search Terminal. In Linux open search area and and find Terminal.
    b. Click on Terminal.   
    c. Write python --version in Terminal and press ENTER.
    d. If this returns a number higher than 3.8, Python is installed and you can continue with step two. If not, continue to the next step.  
    e. Install Python.
    You can download an installer for Windows from <https://www.python.org/downloads/>  

2. #### Setup the project
    a. Navigate to the folder you want to work on. For example, I collect my projects under the folder name 'workspaces'.

        $ cd  <path-to-some-folder>

    b. Clone the git repository.

        $ git clone <GIT_REPO_URL>

    c. Navigate to the project folder

        $ cd  <path-to-project-folder>

    d. Create a virtual environment to prevent system interference and activate it.

        $ pip install virtualenv
        $ python -m venv .env
        $ source /bin/activate

    c. Install the required Python libraries

        $ pip install -r requirements.txt

## The scripts
The scripts, which are inside the `src` folder, are configured to process certain docx formatted files to extract detailed reference information as a CSV file.

#### Running Jupyter notebook  

After activating the virtual environment type in Terminal/Command Prompt  

  * For windows
        $ jupyter-lab
  * For Linux and MacOS
        $ jupyter notebook  

You can run the cells in Jupyter notebook with ctrl^enter  

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
