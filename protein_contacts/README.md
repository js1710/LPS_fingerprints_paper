# Protein Contacts
A description of the options available can be found using the following
```
python protein_contacts.py -h
```
This script is designed to iterate over a file tree such as
```
[Protein]
   [ReLPS repeat 1]
   [ReLPS repeat 2]
   [RaLPS repeat 1]
   ...
```
In the script the names of the protein folders are contained in the list `proteins`, whereas the names of the subdirectories 
are stored in the variable `lps_all`. The output is saved in a pickle script.
