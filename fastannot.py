#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 11:48:21 2021

@author: mazeller
"""

import getopt, sys, os
import pandas as pd

def usage():
    print("python3 fastannot.py -f contigs.fasta -a kaiju.names.out")
    
def main():
    #Get args, handle errors
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:a:", ["help", "fasta=","annot="])
    except getopt.GetoptError as err:
        print(err)  
        usage()
        sys.exit(2)
        
    #Predfine
    fasta = None
    annot = None
    
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("-f", "--fasta"):
            fasta = a
        elif o in ("-a", "--annot"):
            annot = a
        else:
            assert False, "unhandled option"
            
    if (fasta == None or annot == None):
        print("Both arguments are required")
        sys.exit(2)
    i=0    
    #Load annotations as dataframe
    df = pd.read_csv(annot, delimiter = "\t", index_col = 1, header = None)
    #print(df.index)
    #Process contents
    with open(fasta) as f:
        for line in f:
            segment = line.strip();
            
            #Tack onto def line
            hitClass = ""
            if (segment[0] == ">"):
                
                #Check if annotations
                if (segment[1:] in df.index):
                    hitClass = ":" + df.loc[segment[1:], 7 ]
                    hitClass = hitClass.replace(" ","_")
                segment += hitClass
            
            #Dump to screen
            print("{}".format(segment))

if __name__ == "__main__":
    main()