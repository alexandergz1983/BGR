#! /usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: J-Alexander Garcia-Zea

from os import sep
import pandas as pd
import numpy as np
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument('--infile_lcbs')
args = parser.parse_args()

infile_lcbs = args.infile_lcbs

infile_lcbs = pd.read_csv(infile_lcbs, sep="\t")
df_lcbs = pd.DataFrame(infile_lcbs)

TerC_ref = int(input("Enter coordinate terC reference genome: "))
OriC_ref = int(input("Enter coordinate OriC reference genome: "))
size_genome_ref = int(input("Enter Size genome reference genome: "))

# EXAMPLE: J99 genome reference

# OriC:            1556620
# TerC:            685460
# Size genome:     1643831 

####################################################################################################
#######################  Genomic Rearrangement Classifier GRC  #####################################
####################################################################################################

for i in range(2, df_lcbs.shape[1], 2):
    paired = df_lcbs.iloc[:, np.r_[0:2, i:i+2]]
    df_paired_lr = pd.DataFrame(paired)

    E = 25000

    for i in range(0, len(df_paired_lr)):

        if abs(df_paired_lr.iloc[i, 0] - df_paired_lr.iloc[i, 2]) <= E:
            df_paired_lr.loc[i, "IVT-LR"] = "syntenic"

        elif abs(df_paired_lr.iloc[i, 2]) - df_paired_lr.iloc[i, 0] > E and df_paired_lr.iloc[i, 2] < 0:
            df_paired_lr.loc[i, "IVT-LR"] = "inversion/translocation"

        elif df_paired_lr.iloc[i, 0] - abs(df_paired_lr.iloc[i, 2]) > E and df_paired_lr.iloc[i, 2] < 0:
            df_paired_lr.loc[i, "IVT-LR"] = "inversion/translocation"

        elif df_paired_lr.iloc[i, 2] - df_paired_lr.iloc[i, 0] > E:
            df_paired_lr.loc[i, "IVT-LR"] = "translocation"

        elif df_paired_lr.iloc[i, 0] - df_paired_lr.iloc[i, 2] > E:
            df_paired_lr.loc[i, "IVT-LR"] = "translocation"

            if df_paired_lr.iloc[i, 0] - abs(df_paired_lr.iloc[i, 2]) <= E and df_paired_lr.iloc[i, 2] < 0:
                df_paired_lr.loc[i, "IVT-LR"] = "inversion"

    #outfile = df_paired_lr.to_csv('Classify_LCBs_in_IVT_Moran-Lawrence-Roth.csv', sep=";", mode='a' if i > 0 else 'w')
    
    ######################################################################################################
    ### calculate proximity the Inversion, Inversion/translocation, translocation  from oriC and terC ####
    ######################################################################################################

    #df_paired_lr = pd.read_csv('Classify_LCBs_in_IVT_Moran-Lawrence-Roth.csv', sep=";")

    df_ivt = pd.DataFrame(df_paired_lr)
    medium = int(size_genome_ref/2)
    M = int(medium)

    for i in range(0, len(df_ivt)):
        if df_ivt.iloc[i, 4] == "translocation" and df_ivt.iloc[i, 2] and df_ivt.iloc[i, 3] >= 1 and df_ivt.iloc[i, 2] and df_ivt.iloc[i, 3] <= M:
            df_ivt.loc[i, "Origin Replication"] = "proximate-to-TerC"

        elif df_ivt.iloc[i, 4] == "translocation" and df_ivt.iloc[i, 2] and df_ivt.iloc[i, 3] > M and df_ivt.iloc[i, 2] and df_ivt.iloc[i, 3] <= size_genome_ref:
            df_ivt.loc[i, "Origin Replication"] = "proximate-to-OriC"

        elif df_ivt.iloc[i, 4] == "translocation" and df_ivt.iloc[i, 2] and df_ivt.iloc[i, 3] > M and df_ivt.iloc[i, 2] and df_ivt.iloc[i, 3] > 1:
            df_ivt.loc[i, "Origin Replication"] = "proximate-to-OriC"

        elif df_ivt.iloc[i, 4] == "inversion" or "inversion/translocation" and df_ivt.iloc[i, 2]*-1 and df_ivt.iloc[i, 3]*-1 >= 1 and df_ivt.iloc[i, 2]*-1 and df_ivt.iloc[i, 3]*-1 <= M:
            df_ivt.loc[i, "Origin Replication"] = "proximate-to-TerC"

        elif df_ivt.iloc[i, 4] == "inversion" or "inversion/translocation" and df_ivt.iloc[i, 2]*-1 and df_ivt.iloc[i, 3]*-1 > M and df_ivt.iloc[i, 2]*-1 and df_ivt.iloc[i, 3]*-1 <= size_genome_ref:
            df_ivt.loc[i, "Origin Replication"] = "proximate-to-OriC"

    
    ##############################################################################################################################
    ### calculate value for symmetric and asymmetric inversions, translocations, inversions/translocations from oriC and terC ####
    ##############################################################################################################################

    # These lines calculate the lengths of the LCBs for sequences 1 and 2 using the abs() function to obtain the absolute value of the difference between the start and end points of each LCB. 
    # The results are stored in the new columns "size_LCBs_1" and "size_LCBs_2".

    # size_LCBs_1= ∣ df_size_lbs[1] − df_size_lbs[0] ∣
    # size_LCBs_2= ∣ df_size_lbs[3] − df_size_lbs[2] ∣

    # df_size_lbs[1] and df_size_lbs[0] represent the start and end points of the LCB for sequence 1.
    # df_size_lbs[3] and df_size_lbs[2] represent the start and end points of the LCB for sequence 2.
    # ∣⋅∣ denotes the absolute value of the difference between the start and end points of each LCB.
            
    df_size_lbs = df_ivt
    df_size_lbs['size_LCBs_1'] = abs(df_size_lbs.iloc[:, 1] - df_size_lbs.iloc[:, 0])
    df_size_lbs['size_LCBs_2'] = abs(df_size_lbs.iloc[:, 3] - df_size_lbs.iloc[:, 2])
    
    
    # These lines calculate the total length of the LCBs before oriC and after terC for each row of the DataFrame. The apply() function is used in conjunction with a lambda function to apply 
    # conditional logic based on the starting values of the LCBs relative to the OriC_ref and TerC_ref values. 
    # The results are stored in the new columns "total_length_before_oriC" and "total_length_after_terC".
      
    df_size_lbs['total_length_before_oriC'] = df_size_lbs.apply(lambda row: row['size_LCBs_1'] if row[0] < OriC_ref else 0, axis=1) + df_size_lbs.apply(lambda row: row['size_LCBs_2'] if row[0] < OriC_ref else 0, axis=1)
    df_size_lbs['total_length_after_terC'] = df_size_lbs.apply(lambda row: row['size_LCBs_1'] if row[0] > TerC_ref else 0, axis=1) + df_size_lbs.apply(lambda row: row['size_LCBs_2'] if row[0] > TerC_ref else 0, axis=1)

    # In these lines, a threshold of 10% is defined to determine whether the lengths of the LCBs before oriC and after terC are symmetrical or asymmetrical. 
    # The determine_balance() function calculates the difference between the lengths and compares this difference with the threshold. 
    # If the difference is less than the threshold multiplied by the sum of the total lengths, the LCBs are considered symmetric; otherwise, they are considered asymmetric. 
    # The results are stored in the new column "Balance Genomic Rearrangements" of the DataFrame df_size_lbs.

    # The determine_balance(row) function evaluates whether the difference between the total lengths before and after certain reference points, expressed as difference, is less than a certain percentage (threshold) of the sum of these total lengths, expressed as total_length_sum. 
    # If this difference is less than the threshold multiplied by the sum of the total lengths, the function returns "Symmetric"; otherwise, it returns "Asymmetric".
    
    # Let D be the absolute difference between the total lengths before and after certain reference points, and let S be the sum of these total lengths. Let thresholdthreshold be the given threshold, which here is 0.1 (10%).
    # Then, the function determine_balance(row) can be expressed by the following formula:
    
    # If D < threshold x S, then Symmetric, otherwise Asymmetric

    # This means that if the absolute difference between the total lengths before and after certain reference points is less than 10% of the sum of these lengths, 
    # the entity is considered symmetrical; otherwise, it is considered asymmetrical.
    
    threshold = 0.1  # threshold 10%
    def determine_balance(row):
        difference = abs(row['total_length_before_oriC'] - row['total_length_after_terC'])
        total_length_sum = row['total_length_before_oriC'] + row['total_length_after_terC']
        if difference < threshold * total_length_sum:
            return 'Symmetric'
        else:
            return 'Asymmetric'
    
    df_size_lbs['Balance Genomic Rearrangements'] = df_size_lbs.apply(determine_balance, axis=1)

    #print(df_size_lbs)

    df_ivt.to_csv('Classify_IT_TI_BGR.csv', sep=";", mode='a' if i > 0 else 'w', index=False)
