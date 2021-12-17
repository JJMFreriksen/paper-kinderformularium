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

### Prerequisites  

*   Install at least Python 3.8 if you do not have already. There are several ways to install Python depending on the operating system. Alternatively you can download an installer from <https://www.python.org/downloads/>.

### Setup

*   Clone the git repository.

        $ git clone <GIT_REPO_URL>

*   Navigate to the project folder

        $ cd  <path-to-project-folder>

*   Install the required Python libraries

        $ pip install -r requirements.txt

## The scripts
The scripts, which are inside the `src` folder, are configured to process certain docx formatted files to extract detailed reference information as a CSV file.

The Jupyter notebooks inside `notebooks` folder function as an additional usage guides that shows the code in action.

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
- Contains the columns 'record_id', 'title', 'abstract', 'doi', 'final_included'

## License
The content in this repository is published under the [MIT license](https://github.com/asreview/paper-kinderformularium/blob/main/LICENSE).

## Contact
For any questions or remarks, please send an email to asreview@uu.nl.
