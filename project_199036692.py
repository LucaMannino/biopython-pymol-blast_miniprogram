from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=True)
import subprocess
import math
print("\t Welcome to PDB analysis program\n")
print("Please input which pdb file you would like to analise")
file_name=input("")
print("Please input name of the protein\nThis will be equal to the name of the file without .pdb if the file name hasn't been changed")
sequence_id = input("")
#empty lists for the coordinate of the 2 atoms are initiated any int value would be fine these are used for the atom atom distance
coord_1=[1,1,1]
coord_2=[1,1,1]
#reference:
#Biopython Tutorial and Cookbook was used as a reference to understand byopithon modules.
#this information was used particularly for atom_atom_dist(),first_choice() and fasta_sequence() function
#Biopython Tutorial and Cookbook was written by Jeff Chang, Brad Chapman, Iddo Friedberg, Thomas Hamelryck,Michiel de Hoon, Peter Cock, Tiago Antao, Eric Talevich, Bartek Wilczyński.
#Biopython Tutorial and Cookbook: Last Update – 25 May 2020 (Biopython 1.77)
#Information on how to use byopython.PDB was obtained from https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
#the following is a link for the Biopython Tutorial and Cookbook http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec194

#The relevant information of the pdb file is added to  variable structure using PBD parser.get_structure() method. structure variable is an object containing information on different "layers" of the protein from model to chain, from chain to residue, from residue to atom
structure = parser.get_structure(sequence_id, file_name)

#this function is used to return the atom atom distance
def atom_atom_dist(atom1, atom2, sig_figs):
    
#iterates through the different component of structure object in order to obtain the atom coordinates
#this type of navigation was taken from Biopython Tutorial and Cookbook section 11.5
    for model in structure:
        for chain in model:
            for residue in chain:
               for atom in residue:
                 if atom.get_serial_number() == int(atom1):
                  global coord_1
                  coord_1=atom.get_coord()  
                
    for model in structure:
        for chain in model:
            for residue in chain:
               for atom in residue:
                 if atom.get_serial_number() == int(atom2):
                  global coord_2
                  coord_2=atom.get_coord()  
                  
                  
#atom coordinates are used to calculate the distance
#fromula to calculate the distance the result is rounded up with an adjustable significant figure value
    distance =math.sqrt((coord_1[0]-coord_2[0])**2+(coord_1[1]-coord_2[1])**2+(coord_1[2]-coord_2[2])**2)
    sig_figs=int(sig_figs)
    return round(distance, sig_figs)




#This function is used to produce a file with useful information on the protein sequence
def first_choice(sequence_id, file_name,export_file):
#first the file is created with the name chosen by the user
   file= open(export_file,"w+")
   
# the following portion of the code is used to obtain the amino acid sequence
#this method and for loop to obtain the protein sequence was obtained from biopython.org it can be found at: https://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html
   ppb=PPBuilder()
   file.write("Protein Name: " + sequence_id + "\n")
   for pp in ppb.build_peptides(structure):
        lines = str(pp.get_sequence())
#count is initialised at 0, it is needed to count characters in each line in order to slipt the sequence in multiple lines 
   count = 0
   total_aa=str(len(lines))
   file.write("Total amino acid sequence length: " + total_aa + "\n")
#The following portion makes sure that the amino acid sequence is never longer than 35 characters each line or it would be difficult to read
   for character in lines:
        file.write(character)
        count = count +1
        if count%35==0:
           file.write("\n")

#the next for loop will enable to write into the file information about the polypeptide amino acids, this include atom numbers information necessary for atom atom distance calculation
#this model of iteration throgh the biopython information was obtained from Biopython Tutorial and Cookbook section 11.5 
   for model in structure:
        for chain in model:
            file.write(chain.get_id())
            file.write("\n")
            for residue in chain:
                  file.write(str(residue))
                  file.write("\n\n")
                  
                  for atom in residue:
                           file.write("Element:")
                           name = str(atom.get_id())
                           file.write(name[0])
                           file.write("\t\t")
                           
                           file.write("Atom Number:")
                           file.write(str(atom.get_serial_number())) #serial number is the atom number
                           file.write("\n\n")
                           
   file.close()





#this function returns the one letter amino acid sequence necessary for fasta format
#code obtained from the byopython.org website its available at: https://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html
def fasta_sequence(sequence_id, file_name,):
   
   
   ppb=PPBuilder()
#the for loop iteretes through ppb and uses get sequence method to return the 1 letter amino acid sequence from the 3 letter sequence data  
   for pp in ppb.build_peptides(structure):
        sequence = str(pp.get_sequence())
   return(sequence)

#this function produces a pymol script cosisting of color choice for secondary structures
def pymol_script():
	print("Choose a file Name to save the pymol script as")
	name=input()

	file= open(name,"w+")
	helix = ""
	beta = ""
	undef = ""
	print("You will now be able to choose which color you want the secondary strucrues to be")
	while True: 
    		print("Choose which color you want the helix to be\n input 1 for red 2 for green and 3 for yellow\ninput anything else to leave it unchanged" )
    		hel_col=input()
    		if hel_col == "1":
        		helix = "color red, ss h\n"
    		elif (hel_col == "2"):
       	 		helix = "color green, ss h\n"
    		elif hel_col == "3":
        		helix = "color yellow, ss h\n"
    		    
    		print("Choose which color you want the beta strand to be\n input 1 for red 2 for green and 3 for yellown\ninput anything else to leave it unchanged" )
    		beta_col=input()  
    		if beta_col == "1":
        		beta = "color red, ss s\n"
    		elif (beta_col == "2"):
        		beta ="color green, ss s\n"
    		elif beta_col == "3":
        		beta ="color yellow, ss s\n"
        
    		print("In thise which color you want the loops and undefined structures to be: \n input 1 for red 2 for green and 3 for yellow input 4 to restart the color selection\ninput anything else to leave it unchanged" )
    		undef_col=input()
    		if undef_col == "1":
        		undef ="color red, ss l+\n"
    		elif (undef_col == "2"):
        		undef ="color green, ss l+\n"
    		elif undef_col == "3":
        		undef ="color yellow, ss l+\n"
#this is added in case the user changes his mind on the previous selection
    		if undef_col != "4": break;



	

	file.write(helix+beta+undef)
	file.close()
	return(str(name))
#pymol information for the putty and balls and stick configuration was obtained from the pymolwiki.org
#this function allows the user to choose how to visualise the protein
#3 options are given 
def pymol_how_to_visualise(name):
#the file modified is the same file created in the pymol script function
	file =open(name, "a")
	print("how would you like to visualise the protein:\ninput1 for putty visualisation,\n 2 for balls and sticks configuration,\n 3 for classic cartoon visualisation,\n input anything else to leave unchanged")
	choice=input("")
	if choice == "1":
		file.write("show cartoon\n")
		file.write("cartoon putty\n")
		file.write("hide lines\n")
		file.write("hide sticks\n")
		file.write("hide nonbonded\n")
		file.write("hide ribbon\n")
		file.write("unset cartoon_smooth_loops\n")
		file.write("unset cartoon_flat_sheets\n")
		file.close()
	if choice == "2":
		file.write("show sticks\n")
		file.write("show spheres\n")
		file.write("hide lines\n")
		file.write("hide cartoon\n")
		file.write("hide nonbonded\n")
		file.write("hide ribbon\n")
		file.write("set stick_radius, 0.3, (all)\n")
		file.write("set sphere_scale, 0.3, (all)\n")
		file.close()
	if choice == "3":
		file.write("show cartoon\n")
		file.write("hide sticks\n")
		file.write("hide spheres\n")
		file.write("hide lines\n")

		file.write("hide nonbonded\n")
		file.write("hide ribbon\n")
	file.close()
	return()

#this while loop is used to iterate through the 4 different program functionalities
#making it a while true loop ensures that the user can choose to use muliple options one after the other

while True:
	print("\tInput 0,1,2,3 or 4 depending which task you would like to complete,\n")
	print("0. To exit the program")
	print("1. Produce an outline of a pdb file")
	print("2. Atom-Atom distance")
	print("3. Produce a pymol script")
	print("4. Local BLASTp")
	print("5. Choose a new protein to analyse")
	choice = input("")
#the first option enable the user to save an outline of the .pdb file with usefull information in a separate file which name is chosen by the user
	if choice =="1":
		print("\t Outline of pdb file")
		print("The protein outline data obtained will be exported as a file in the working directory")
		export_file=input("Input file name\n")
		first_choice(sequence_id, file_name, export_file)
	elif choice =="2":
#while True loop is used to ensure the user is inputting an appropiate value
#if the user wants to exit this portion of the program and select another task they can input 0
		while True:
			print("\tAtom-Atom Distance")
			atom1=input("first ATOM number, input 0 to go back to task selection\n")
			if atom1=="0": break;
			atom2=input("second ATOM number,input 0 to go back to task selection\n")
			if atom2=="0": break;
			sig_figs=input("how many significant figures would you like the result to be visualised?\ninput 0 to go back to task selection\n")
			if sig_figs=="0": break;			
    
			if atom1.isdigit() ==False or atom2.isdigit() ==False or sig_figs.isdigit() ==False:
        			print("wrong input")
			elif int(atom1) <=1 or int(atom2) <=1 or int(sig_figs) <=1 :
       				print("invalid input")
			else:
				print("Distance is:\n")       				
				print(atom_atom_dist(atom1, atom2, sig_figs))
				break;
	elif choice =="3":
		print("\tPymol Script")
#the function pymol script returns the file name as this is needed to allow the following function to modify the script file
		name = pymol_script()
		pymol_how_to_visualise(name)
		locate_file="locate -br '^%s'"%(name)
		
		print("You can use the pymol script by typing @filepath/filename in the pymol command line")
		print("The information needed for the path and file name should appear on the first line below:")
#locate file is used as a shell command to return information including the file path and name
		subprocess.call(locate_file, shell=True)
	elif choice =="4":
#this last option is to run a local blastp
#the user would have previously needed to download a database
        	print("\tLocal BLASTp\n")
        	print("For the local BLASTp first a query fasta format file needs to be created.")
        	fasta = input("the protein sequence will be saved as an external text file \n input the name you would like so save this query file as ")
        	file= open(fasta,"w+")
        	out_putfile =input("please imput the name you want your result to be save as, .txt file format is advised\n")
#choice for the blastp results file name
        	database_name =input("please input the database name\n")
        	file.write(">Sequence" + sequence_id + "\n") 
        	file.write(fasta_sequence(sequence_id,file_name))
        	file.close()

# command_0 initialises the protein database this is a necessary preliminary step
        	command_0 = "makeblastdb -in %s -dbtype prot" % (database_name)
# command  this command performs the blastp 
        	command = "blastp -query %s -db %s -out %s" %(fasta, database_name,out_putfile)

        	subprocess.call(command_0, shell=True)
        	subprocess.call(command, shell=True)
#this choice anables the user to choose another protein file to analise
	elif choice =="5":
		print("Choose new .pdb file")
		file_name=input("")
		print("please input name of the protein\nThis will be equal to the name of the file without .pdb if the file name hasn't been changed")
		sequence_id = input("")
#empty lists for the coordinate of the 2 atoms are initiated any int value would be fine these are used for the atom atom distance
		coord_1=[1,1,1]
		coord_2=[1,1,1]
#this line initialise the variable structure object the PBD parser.get_structure() method is used
#this method was obtained from Byopython 
		structure = parser.get_structure(sequence_id, file_name)
#if the user imputs 0 the loop will break and the program will stop running
	elif choice =="0":
		break;
#if the input was not a value between 0 and 4 the loop will start from the beginning
	else:
		print("input appropiate value")

































    
