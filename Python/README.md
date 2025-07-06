**NOTE:** **AlphaFold model files store per-residue confidence scores (pLDDT) in the B-factor field and can be used directly in ChimeraX. Therefore, the file conversion step is completely optional but follows as below**

Copy the python script in .text and edit input and output file names

Save the file as .py 

open commandprompt

#install the package below

pip install biopython

#Set working directory

C:\Users\...\Main_folder_containing_the_downloaded_AlphFold_subfolder> cd "C:\Users\...\Main_folder_containing_the_downloaded_AlphFold_subfolder\The_AlphFold_folder"


#Upload the cif to pdb conversion script file (Save the conversion script file in the same AlphaFold folder)

C:\Users\...\Main_folder_containing_the_downloaded_AlphFold_subfolder\The_AlphFold_folder> python name_of_the_conversion_script_file.py
