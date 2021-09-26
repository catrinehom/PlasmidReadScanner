#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 12:26:29 2021

@author: catrinehom
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
"""
########################################
#### CODE GRAVEYARD
########################################

### Read file
read_file = "/Users/catrinehom/Desktop/plasmid_read_scanner/sequencing_summary_read_lenght_only_02.txt"
reads = pd.read_csv(read_file, sep='\t', index_col=False)
counts = reads["sequence_length_template"].value_counts(sort=False)
#counts = counts.reset_index()

### plot
sns.lineplot(data = counts, y = 'sequence_length_template', x='index')
plt.xlim(0, 30000)

### data wrangling
### I am trying out with pandas dfs data structure

counts = counts.sort_values(by=['index'])

counts_dict = dict()
inc = 10
ind = 0

for row in counts.iterrows():
    if inc < 10:
        counts_dict[ind] += row[1]["sequence_length_template"]
        inc += 1
    else:
        ind = row[1]["index"]
        counts_dict[ind] = row[1]["sequence_length_template"]
        inc = 0

counts_dict = dict()
inc = 0
for i in range(0,max(counts.index)):
    if round(i/10)*10 in counts_dict:
        try:
            counts_dict[round(i/10)*10] += counts["sequence_length_template"][i]
        except:
            counts_dict[round(i/10)*10] += 0
    else:
        try:
            counts_dict[round(i/10)*10] = counts["sequence_length_template"][i]
        except:
            counts_dict[round(i/10)*10] = 0
"""            
########################################
#### Approach 1 
########################################
# Første meget naive tilgang.
# Inddel pos i 50 intervaller
# Regn kun fra øverste punkt
# Hvis punktet er mere end 100 reads mere end punktet før, giv output

read_file = "/Users/catrinehom/Desktop/plasmid_read_scanner/sequencing_summary_pass_read_length_and_barcode.txt"
f = open(read_file, "r")
f.readline()

counts = dict()

# Convert file to dict
all_barcodes = dict()
# Split into intervals of 10 in pos, and sum the no. of reads.
for line in f:
    line = line.split()
    barcode = line[1]
    pos = round(int(line[0])/20)*20
    if barcode in all_barcodes:
        if pos in all_barcodes[barcode]:
            all_barcodes[barcode][pos] += 1
        else:
            all_barcodes[barcode][pos] = 1
    else:
        all_barcodes[barcode] = dict()

for barcode in all_barcodes: 
    counts = all_barcodes[barcode]        
            
    prev = 10000
    max_value = max(counts.values())
    flag = False
    pot_plasmids = list()
    
    for k in sorted(counts.keys()):
        # See if top value is found
        if counts[k] == max_value:
            flag = True
        # If we are passed top value, find 100+ reads
        if flag == True:
            if counts[k] > prev + 50:
                #print((k, counts[k]))
                pot_plasmids.append([k, counts[k]])
            prev = counts[k]
       
    # Convert to dataframes for plotting     
    df = pd.DataFrame.from_dict(counts, orient='index', columns=['count'])
    df['pos'] = df.index
    pot_plasmids_df = pd.DataFrame(pot_plasmids, columns =['pos', 'count'])
    
    # Plot
    sns.lineplot(data = df, y = 'count', x='pos', palette="red")
    sns.scatterplot(data=pot_plasmids_df, x="pos", y="count")
    plt.xlim(0, 60000)
    
    plt.savefig("/Users/catrinehom/Desktop/plasmid_read_scanner/plots/1_større_end_20_interval_"+'{}.png'.format(barcode), format='png', dpi=900)
    plt.close()


########################################
#### Approach 1 
########################################
# Første meget naive tilgang.
# Inddel pos i 50 intervaller
# Regn kun fra øverste punkt
# Hvis punktet er mere end 100 reads mere end punktet før, giv output

read_file = "/Users/catrinehom/Desktop/plasmid_read_scanner/sequencing_summary_pass_read_length_only_barcode02.txt"
f = open(read_file, "r")
f.readline()

counts = dict()

# Split into intervals of 50 in pos, and sum the no. of reads.
for line in f:
    pos = round(int(line)/50)*50
    if pos in counts:
        counts[pos] += 1
    else:
        counts[pos] = 1
        
        
prev = 10000
max_value = max(counts.values())
flag = False
pot_plasmids = list()

for k in sorted(counts.keys()):
    # See if top value is found
    if counts[k] == max_value:
        flag = True
    # If we are passed top value, find 100+ reads
    if flag == True:
        if counts[k] > prev + 50:
            print((k, counts[k]))
            pot_plasmids.append([k, counts[k]])
        prev = counts[k]
   
# Convert to dataframes for plotting     
df = pd.DataFrame.from_dict(counts, orient='index', columns=['count'])
df['pos'] = df.index
pot_plasmids_df = pd.DataFrame(pot_plasmids, columns =['pos', 'count'])

# Plot
sns.lineplot(data = df, y = 'count', x='pos')
sns.scatterplot(data=pot_plasmids_df, x="pos", y="count")
plt.xlim(0, 60000)

########################################
#### Approach 2 
########################################

import statistics

# Regn kun fra øverste punkt
# Tag gennemsnittet af fx tyve punkter
# Outliers er skal være større end gennemsnittet af de 50 foregående + 2x standardafvigelse
      
read_file = "/Users/catrinehom/Desktop/plasmid_read_scanner/sequencing_summary_pass_read_length_and_barcode.txt"
f = open(read_file, "r")
f.readline()

# Convert file to dict
all_barcodes = dict()
# Split into intervals of 10 in pos, and sum the no. of reads.
for line in f:
    line = line.split()
    barcode = line[1]
    pos = round(int(line[0])/30)*30
    if barcode in all_barcodes:
        if pos in all_barcodes[barcode]:
            all_barcodes[barcode][pos] += 1
        else:
            all_barcodes[barcode][pos] = 1
    else:
        all_barcodes[barcode] = dict()

for barcode in all_barcodes: 
    counts = all_barcodes[barcode]
    print(barcode)
    max_value = max(counts.values())
    flag = True
    ten_list = []
    pot_plasmids = list()
    
    for k in sorted(counts.keys()):
        # See if top value is found
        if counts[k] == max_value:
            flag = True
    
        # If we are passed top value
        if flag == True:
            # Remove noise in the end by requiring a depth
            if len(ten_list) == 30:
                # Get average of 10 pos and std dev
                ten_average = sum(ten_list)/30
                teen_std_dev = statistics.stdev(ten_list)
                if counts[k] > ten_average + 2*teen_std_dev + 5:
                    #print((k, counts[k]))
                    pot_plasmids.append([k, counts[k]])
                    #print(ten_list, ten_average, teen_std_dev)
        # Add the new pos and remove the oldest one
        ten_list.append(counts[k])
        if len(ten_list) > 30:
            ten_list.pop(0)
            
    # Convert to dataframes for plotting     
    df = pd.DataFrame.from_dict(counts, orient='index', columns=['count'])
    df['pos'] = df.index
    pot_plasmids_df = pd.DataFrame(pot_plasmids, columns =['pos', 'count'])
    
    # Plot
    sns.lineplot(data = df, y = 'count', x='pos', color='green')
    sns.scatterplot(data=pot_plasmids_df, x="pos", y="count")
    plt.xlim(0, 60000)    
    #plt.show()
     
    
    plt.savefig("/Users/catrinehom/Desktop/plasmid_read_scanner/plots/"+"2_gennemsnit_20_interval_"+'{}.png'.format(barcode), format='png', dpi=900)
    plt.close()


########################################
#### Approach 3
########################################

import statistics
import scipy

# Regn kun fra øverste punkt
# Tag gennemsnittet af fx tyve punkter
# Outliers er skal være større end den linære regression + 2x standardafvigelse
      
read_file = "/Users/catrinehom/Desktop/plasmid_read_scanner/sequencing_summary_pass_read_length_and_barcode.txt"
f = open(read_file, "r")
f.readline()

# Convert file to dict
all_barcodes = dict()
# Split into intervals of 10 in pos, and sum the no. of reads.
for line in f:
    line = line.split()
    barcode = line[1]
    pos = round(int(line[0])/50)*50
    if barcode in all_barcodes:
        if pos in all_barcodes[barcode]:
            all_barcodes[barcode][pos] += 1
        else:
            all_barcodes[barcode][pos] = 1
    else:
        all_barcodes[barcode] = dict()

for barcode in all_barcodes: 
    counts = all_barcodes[barcode]
    print(barcode)
    max_value = max(counts.values())
    flag = False
    y_list = []
    x_list = []
    pot_plasmids = list()

    for k in sorted(counts.keys()):
        # See if top value is found
        if counts[k] == max_value:
            print(counts[k], k,  max_value)
            flag = True
    
        # If we are passed top value
        if flag == True:
            # Remove noise in the end by requiring a depth
            if len(x_list) == 30:
                # Get linear regression of 10 pos and std dev
                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x_list,y_list)
                std_dev = statistics.stdev(y_list)
                if counts[k] > (k*slope) + intercept + 2*std_dev + 5:
                    pot_plasmids.append([k, counts[k]])
                    #print(x_list,y_list, counts[k], (k*slope + intercept), 4*std_err + 10)
        # Add the new pos and remove the oldest one
        y_list.append(counts[k])
        x_list.append(k)
        if len(x_list) > 30:
            y_list.pop(0)
            x_list.pop(0)
            
            
    # Convert to dataframes for plotting     
    df = pd.DataFrame.from_dict(counts, orient='index', columns=['count'])
    df['pos'] = df.index
    pot_plasmids_df = pd.DataFrame(pot_plasmids, columns =['pos', 'count'])
    
    # Plot
    sns.lineplot(data = df, y = 'count', x='pos', color='orange')
    sns.scatterplot(data=pot_plasmids_df, x="pos", y="count")
    plt.xlim(0, 60000)    
    #plt.show()
    
    
    plt.savefig("/Users/catrinehom/Desktop/plasmid_read_scanner/plots/"+"3_lineær_regression_50_interval_"+'{}.png'.format(barcode), format='png', dpi=900)
    plt.close()
