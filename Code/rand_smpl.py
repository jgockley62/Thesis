#! python

##Written by Jake Gockley, Yale University, 2017
##This program downsamples Tag Counts 
#Run Example: python rand_smpl.py <Max_Tag_Number_to_use>

#Must have R scripts in run directory: Modular_Quantile_Normalizer.R
#Required support files in run directory:
#/Users/Sazerac/Dropbox/MPRA_Data/2016_12_04_Tag_DownSampling/MPRA_minP_*bc.txt

import sys;
import os;
import re;
import random;
import numpy as np;
import subprocess

##Go to file directory


#Get File lists
Files = os.popen('ls MPRA_minP_*bc.txt').read().split('\n')
del Files[-1]

for FILE in Files:

	#Will permute both replicates
	TotalLines = {}
	
	#Figure out Max Tags to permute per-permutation for given file	
	TagNumb = int(FILE.replace("MPRA_minP_","").replace("bc.txt",""))
	
	#Load Main File
	header = 0
	temp = os.getcwd()
	#print str(temp+"/"+FILE)
	Refs = str(temp+"/"+FILE)
	F0 = open(Refs)
	for line in F0:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
	
		Entry= '\t'.join(map(str,lst))
		#Print Header
		if lst[0] == "Tag":
			header = Entry
		 
		else:
			#Append to Fragment list
			if lst[1] in TotalLines:
				TotalLines[lst[1]].append(lst)
			
				#delete the empty initialized array value
				if not TotalLines[lst[1]][0]:
					del TotalLines[lst[1]][0]
				else:
					pass				
		
			#Create nested Fragment list
			else:
				TotalLines[lst[1]] = [[]]
				TotalLines[lst[1]].append(lst)
				
	F0.close();

	#Seed Perm Range and range of Barcode Refs to sample
	Samps = range(1,(int(TagNumb)/2)+1)
	BC_Ranges = range(0,TagNumb) 
	
	#Preform Permutations
	for N in Samps:
		#Open Output and print header
		OUTA = open("Permute_A.txt", "w")
		OUTB = open("Permute_B.txt", "w")
	
		print >>OUTA, header
		print >>OUTB, header

		#Permute N barcodes from each fragment Array
		for key in TotalLines:
	
			#Randomize the index numbers (2xN) for 2 replicates
			double = 2 * int(N)
			Inxs = np.random.choice(BC_Ranges, double, replace=False)
			#Print out the barcodes corresponding to the randomized index numbers
			count = 0
			for val in Inxs:
				if count < int(N):
					#print '\t'.join(TotalLines[key][val])
					print >>OUTA, '\t'.join(TotalLines[key][val])
					count += 1
				else:
					#print '\t'.join(TotalLines[key][val])
					print >>OUTB, '\t'.join(TotalLines[key][val])

		OUTA.close();
		OUTB.close();			

		# Build subprocess command And run R script (output is saved as #_Barcodes.txt)
		r_file_name = "Modular_Quantile_Normalizer.R"
		ppath = '/Library/Frameworks/R.framework/Resources/bin/R CMD BATCH '
		args = ''.join(["'","--args"," ",str(N),"'"])
		
		proc = subprocess.Popen("%s %s %s" % (ppath, args, r_file_name), stdout=subprocess.PIPE, shell=True)
		output = proc.stdout.read()
	
		#Clean up Permutation Files
		os.system("rm Permute_*.txt")
	
	CMD = "cat $(ls -tr *_Barcodes.txt) > "+str(TagNumb)+"_Barcode_Perms.txt"
	os.system(CMD)
	os.system("rm *_Barcodes.txt")	

	
	
	
