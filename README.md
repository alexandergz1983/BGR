This Python script is a Genomic Rearrangement Classifier (GRC) that analyzes conserved block location data files (LCBs: Mauve coordinates) to determine different types of genomic rearrangements, such as translocations and inversions.  
First, the script reads a CSV input (.lcbs) file containing information about the LCBs. Then, it prompts the user to enter reference genome coordinates, such as the OriC and TerC coordinates (origin and replication terminus), as well as the size of the reference genome.  
It then classifies LCBs into different types of genomic rearrangements, such as translocations and inversions, using criteria based on the distance between the start and end points of LCBs.  
Subsequently, it calculates the proximity of genomic rearrangements with respect to the OriC and TerC landmarks.  
Finally, it determines whether the LCBs exhibit symmetric or asymmetric rearrangements based on the length difference before and after the landmarks.  
The script produces an output CSV file with the classification of LCBs and genomic rearrangement information.

use:
python BGR.py --infile_lcbs africa.lcbs.
the script will request:
TerC(ej: 685460), OriC(ej:1556620), Size genome(1643831) 
