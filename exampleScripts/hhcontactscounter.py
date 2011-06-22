#!/usr/bin/env python


import sys, os, getopt
from p3d import protein as protein
from p3d.library import convert321, calcHydropathy,calcSideChainSurface,calcMass
from collections import defaultdict as ddict
import optparse

#Bondi van-der-Waals radii
#vdw_radii ={'H':1.20,'C':1.70,'N':1.55,'O':1.52,'S':1.80}
vdw_radii ={'C':1.70,'N':1.55,'O':1.52,'S':1.80}
C_atoms=[]
N_atoms=[]
O_atoms=[]
S_atoms=[]


# based on Kyte-Doolittle hydrophobicity scale
def isHydrophobic(residue):
	if calcHydropathy(residue) > 0:
		return True
	else:
		return False




def estimateMinContacts(Protein=None,ResidueAtom=None,targetChain=None,DistanceMatrix=None,maxDistance=None, intraChainExcludeRange=None, Model=1, countHydrophobic=False,protComplex=False,atomwiseContacts=False,VDWoverlap=None):
	# atoms of focal residue
	residueAtoms = Protein.query('protein and resid {0} and chain {1} and model {2}'.format(ResidueAtom.resid, ResidueAtom.chain, Model))
	
	# filter residues within the excluded intra-chain distance range (upstream and downstream)
	exclude_string = ""
	if intraChainExcludeRange != None:
		resid = ResidueAtom.resid
		prot_len = len(Protein.atoms)
		exclude_list=[]
		for i in range(resid-intraChainExcludeRange, resid+intraChainExcludeRange+1):
			if i >= 0 and i < prot_len:
				exclude_list.append(i)
		exclude_string=str(exclude_list).strip("[]")
	
	# get surrounding atoms for this residue (within maxDistance) in target chain
	if protComplex:
		surrounding1 = Protein.query('protein and within {0} of'.format(maxDistance)\
			,residueAtoms,'and chain {1} not resid {1} and model {2}'.format(ResidueAtom.chain,exclude_string, Model))
		surrounding2 = Protein.query('protein and within {0} of'.format(maxDistance)\
			,residueAtoms,'and chain {1} not resid {1} and model {2}'.format(targetChain,exclude_string, Model))
		surrounding=surrounding1+surrounding2
	else:
		surrounding = Protein.query('protein and within {0} of'.format(maxDistance)\
		,residueAtoms,'and chain {0} and not resid {1} and model {2}'.format(targetChain, exclude_string, Model))
	
	#add a new contact to DistanceMatrix as ordered pair of focal residue position and position of surrounding residue 
	for residAtom in residueAtoms:
		if not countHydrophobic or (countHydrophobic and isHydrophobic(convert321(residAtom.aa))):
			for atom in surrounding:
				if (atom.chain == residAtom.chain and atom.resid > residAtom.resid) or atom.chain != residAtom.chain: # to avoid counting the same contact twice (in reversed order)
					
					AA1 = convert321(residAtom.aa)
					AA2 = convert321(atom.aa)
					
					d = atom.distanceTo(residAtom)
					
					if atomwiseContacts:
						atom1 = residAtom.atype.strip()
						atom2 = atom.atype.strip()
						
						ordered_tuple = (residAtom.resid, atom.resid,atom1,atom2,AA1,AA2) #if targetChain == 'B' else (atom.resid,residAtom.resid)
						
						if VDWoverlap != None:
							# the vdW radii of both atoms
							r1=0.0
							r2=0.0
							
							
							if atom1[0] in vdw_radii:
								r1 = vdw_radii[atom1[0]]
								#print "atom1 found"
							else:
								print "atom1 type not recognised:",atom1
								continue
							
							if atom2[0] in vdw_radii:
								r2 = vdw_radii[atom2[0]]
								#print "atom2 found"
							else:
								print "atom2 type not recognised:",atom2
								continue
							
							vdw_overlap = d - (r1+r2)
							#print vdw_overlap
							if vdw_overlap >= VDWoverlap:
								d = vdw_overlap
							else:
								continue
						
						DistanceMatrix[ordered_tuple] = d
					else:
						ordered_tuple = (residAtom.resid,atom.resid,AA1,AA2) #if targetChain == 'B' else (atom.resid,residAtom.resid)
						#print ordered_tuple
						#print d, DistanceMatrix[ordered_tuple]
						if DistanceMatrix[ordered_tuple] > d or DistanceMatrix[ordered_tuple] == 0:
							DistanceMatrix[ordered_tuple] = d
							#print AA1,hydropath1,AA2,hydropath2
	return DistanceMatrix




# determines which types of atoms are present in the surrounding of specified residue within maxDistance (angstrom)
def estimateAtomicSurrounding(Protein=None,ResidueAtom=None,targetChain=None,DistanceMatrix=None,maxDistance=None, intraChainExcludeRange=None, Model=None, countHydrophobic=False,protComplex=False,VDWoverlap=None):
	# atoms of focal residue
	residueAtoms = Protein.query('protein and resid {0} and chain {1} and model {2}'.format(ResidueAtom.resid, ResidueAtom.chain, model))
	
	# filter residues within the excluded intra-chain distance range (upstream and downstream)
	exclude_string = ""
	if intraChainExcludeRange != None:
		resid = ResidueAtom.resid
		prot_len = len(Protein.atoms)
		exclude_list=[]
		for i in range(resid-intraChainExcludeRange, resid+intraChainExcludeRange+1):
			if i >= 0 and i < prot_len:
				exclude_list.append(i)
		exclude_string=str(exclude_list).strip("[]")	
	
	# get surrounding atoms for this residue (within maxDistance) in target chain
	
	if protComplex:
		surrounding1 = Protein.query('protein and within {0} of'.format(maxDistance)\
			,residueAtoms,'and chain {1} not resid {1} and model {2}'.format(ResidueAtom.chain,exclude_string, model))
		surrounding2 = Protein.query('protein and within {0} of'.format(maxDistance)\
			,residueAtoms,'and chain {1} not resid {1} and model {2}'.format(targetChain,exclude_string, model))
		surrounding=surrounding1+surrounding2
	else:
		surrounding = Protein.query('protein and within {0} of'.format(maxDistance)\
		,residueAtoms,'and chain {0} and not resid {1} and model {2}'.format(targetChain, exclude_string, model))
	nAtoms=len(surrounding)
	
	print ResidueAtom.resid, ResidueAtom.aa,ResidueAtom.chain, nAtoms, calcHydropathy(convert321(ResidueAtom.aa)),
	#for residAtom in residueAtoms:
	#	tmp = residAtom.atype
	#	print tmp,
	#print "|",
	C_counter = 0
	CA_counter =0
	non_C_counter = 0
	h_counter=0
	p_counter=0
	self_counter=0
	for atom in surrounding:
		tmp = atom.atype
		tmp_aa = atom.aa
		
		#if atom.resid == ResidueAtom.resid:
		#	self_counter+=1
		
		hP = calcHydropathy(convert321(tmp_aa))
		if hP > 0 and not "CA" in tmp:
			h_counter += 1
		else:
			p_counter += 1
		
		if "CA" in tmp:
			CA_counter += 1
		
		if "C" in tmp:
			C_counter += 1
		else:
			non_C_counter += 1
		#print tmp,
	if non_C_counter == 0:
		ratio=C_counter
	else:
		ratio = (C_counter/float(non_C_counter))
	print C_counter,non_C_counter,h_counter,p_counter,CA_counter


def hh_contact_count(pdbfile, chain1, chain2,maxDistance, excludeRange,model=1,print_contact_list=False,print_summary=False,count_total_contacts=False,atomic=False,protComplex=False,atomwiseContacts=False,VDWoverlap=None,twoCarbonCount=False):
	#here, the protein structure is converted into a python object (see p3d docs)
	pdb=None
	#try:
	if chain1 != chain2:
		pdb = protein.Protein(pdbfile,chains=[chain1,chain2])
	else:
		pdb = protein.Protein(pdbfile,chains=[chain1])
	#except ValueError as ve:	
	#	print "a value error occurred"
	#	print ve.args
	#	print ve.message
	#	
	#	sys.exit(1)
	#except IndexError:
	#	print "an index  error occurred"
	#	
	#	
	#	sys.exit(1)
	#for line in pdb.info():
	#	print(line)

	
	# this data structure will hold a dicionary where the keys are ordered pairs
	# of the residue positions that are in contact, and the value is the distance
	# between those residues in angstrom. (--> closest distance between any atoms of the residues)
	DistanceMatrix = ddict(float)
	
	''' query Set alphas to have one atom per residue from both chains '''
	chain1_alphas = pdb.query('chain {0} and alpha and model {1}'.format(chain1, model))
	chain2_alphas = pdb.query('chain {0} and alpha and model {1}'.format(chain2, model))
	
	#print "#+#+#+#+#+#+#", len(chain1_alphas),len(chain2_alphas)
	
	if protComplex:
		complex_alphas = chain1_alphas + chain2_alphas
	else:
		complex_alphas = chain1_alphas
	
	chain1_hydrophobic_counter = 0
	chain2_hydrophobic_counter = 0
	chain1_polar_counter = 0
	chain2_polar_counter = 0
	contacts_counter = 0
	zero_C_counter = 0
	one_C_counter = 0
	two_C_counter = 0
	
	for alpha in complex_alphas: #these are chain1 alphas if protComplex=False
		if atomic:
			estimateAtomicSurrounding(Protein=pdb,\
							ResidueAtom=alpha,\
				 			targetChain=chain2,\
				 			DistanceMatrix=DistanceMatrix,\
				 			maxDistance=maxDistance,\
				 			intraChainExcludeRange=excludeRange,\
				 			Model=model, \
				 			countHydrophobic=False,\
				 			atomwiseContacts=atomwiseContacts)
		# counting all contacts?
		if count_total_contacts:
			if isHydrophobic(convert321(alpha.aa)):
				chain1_hydrophobic_counter += 1
			else: 
				chain1_polar_counter += 1
			# add contacts for residue of this alpha carbon atom to dictionary
			DistanceMatrix = estimateMinContacts(Protein=pdb,ResidueAtom=alpha,\
				 targetChain=chain2,DistanceMatrix=DistanceMatrix,\
				 maxDistance=maxDistance,intraChainExcludeRange=excludeRange,Model=model, countHydrophobic=False,atomwiseContacts=atomwiseContacts,VDWoverlap=VDWoverlap)
		else:
			# or counting only contacts between hydrophobic residues
			if isHydrophobic(convert321(alpha.aa)):
				chain1_hydrophobic_counter += 1
				#append contacts to DistanceMatrix
				DistanceMatrix = estimateMinContacts(Protein=pdb,ResidueAtom=alpha,\
					 targetChain=chain2,DistanceMatrix=DistanceMatrix,\
					 maxDistance=maxDistance,intraChainExcludeRange=excludeRange,Model=model,countHydrophobic=True,atomwiseContacts=atomwiseContacts,VDWoverlap=VDWoverlap)
			else: 
				chain1_polar_counter += 1
	
	
	if print_summary:
		for alpha in chain2_alphas:
			if isHydrophobic(convert321(alpha.aa)):
				chain2_hydrophobic_counter += 1
				#DistanceMatrixHH = estimateMinContacts(Protein=pdb,ResidueAtom=alpha,\
				#									 targetChain=chain1,DistanceMatrix=DistanceMatrixHH,\
				#									 maxDistance=maxDistance,intraChainExcludeRange=excludeRange,Model=model)
			else:
				chain2_polar_counter += 1
	
	#column headers for output table of contacts
	if print_contact_list:
		if atomwiseContacts:
			print "\t".join(["pos_chain1","pos_chain2","contact_dist(A)","atom1","atom2","AA1","AA2","#C-atoms/contact"])
		else:
			print "\t".join(["pos_chain1","pos_chain2","contact_dist(A)","AA1","hydropathy1","surface1(nm^2)","mass1(Da)","AA2","hydropathy2","surface2(nm^2)","mass2(Da)"])
	
	
	for coords,value in DistanceMatrix.items():
		if atomwiseContacts:
			AA1 = coords[4]
			AA2 = coords[5]
			atom1 = coords[2]
			atom2 = coords[3]
			C_atoms=0
			if "C" in atom1:
				C_atoms += 1
			if "C" in atom2:
				C_atoms += 1
			if C_atoms == 0:
				zero_C_counter += 1
			elif C_atoms == 1:
				one_C_counter += 1
			elif C_atoms == 2:
				two_C_counter += 1
			if print_contact_list:
				print "\t".join([str(coords[0]),str(coords[1]),str(value),str(atom1),str(atom2),AA1,AA2,str(C_atoms) ])
		else:
			if print_contact_list:
				AA1 = coords[2]
				AA2 = coords[3]
				hydropath1 = calcHydropathy(AA1)
				hydropath2 = calcHydropathy(AA2)
				surf1 = calcSideChainSurface(AA1)
				surf2 = calcSideChainSurface(AA2)
				mass1 = calcMass(AA1)
				mass2 = calcMass(AA2)
				
				print "\t".join([str(coords[0]),str(coords[1]),str(value),AA1,str(hydropath1),str(surf1),str(mass1),AA2,str(hydropath2),str(surf2),str(mass2) ])
		contacts_counter += 1
	
	
	
	if print_summary:
		print "===================================================================="
		print "PDB file:", pdbfile
		print "chain {0} vs chain {1}".format(chain1, chain2)
		if count_total_contacts:
			print "total contacts counted:", contacts_counter
		else:
			print "hydrophobic contacts counted:", contacts_counter
		if atomwiseContacts:
			print "C/C atom contacts:",two_C_counter
		print "--------------------------------------------------------------------"
		print "chain1:",chain1
		print "length1:", len(chain1_alphas)
		print "polar1:",str(chain1_polar_counter)
		print "hydrophobic1:",str(chain1_hydrophobic_counter)
		ch1_seq=""
		ch1_dict={}
		for alpha in chain1_alphas:
			 ch1_dict[alpha.resid] = convert321(alpha.aa)
		
		for i in range(1,len(chain1_alphas)+1):
			ch1_seq += ch1_dict[i]
		
		print "sequence1:",ch1_seq
		
		
		
		if chain1 != chain2:
			print "--------------------------------------------------------------------"
			print "chain2:",chain2
			print "length2:", len(chain2_alphas)
			print "polar2:",str(chain2_polar_counter)
			print "hydrophobic2:",str(chain2_hydrophobic_counter)
			ch2_seq=""
			ch2_dict={}
			for alpha in chain2_alphas:
				 ch2_dict[alpha.resid] = convert321(alpha.aa)
		
			for i in range(1,len(chain2_alphas)+1):
				ch2_seq += ch2_dict[i]
		
			print "sequence2:",ch2_seq
			
		print "===================================================================="
	
	if twoCarbonCount:
		return two_C_counter
	else:
		return contacts_counter



####################################
# MAIN PROGRAM #
####################################

if (__name__ == '__main__'):
	####################################
	# HANDLE ARGUMENTS, OPTIONS & HELP #
	####################################
	description="HYDROPHOBIC CONTACTS COUNTER (by Tobias Sikosek) uses p3d python module: \nFufezan & Specht (2009) p3d - python module for structural bioinformatics. BMC Bioinformatics, 10:258.  It Counts the number of contacts between hydrophobic residues, either within the same PDB chain or between different chains."
	version="2011-05-26; written by Tobias Sikosek"
	usage="\npython %prog <PDB file> <chain1> <chain2> <maxDistance to take into account> <range of up/downstream residues excluded> [ <nmr model number>] [OPTIONS]\n\nNotes:\nChain1 and chain2 have to be identical for counting contacts within the same structure.\nWithout options, only the number of (hydrophobic or total) contacts between chains 1 and 2 is returned.\n.Distance is calculated per residue pair.\n"
	op=optparse.OptionParser(description=description, usage=usage, version=version)
	op.add_option('-c','--contactList', action="store_true", dest="print_contact_list",help="display list of contacting residue pairs")
	op.add_option('-t','--totalContacts', action="store_true", dest="count_total_contacts",help="instead of only hydrophobic residues, consider all residues")
	op.add_option('-s','--summary', action="store_true", dest="print_summary",help="display summary info on input protein")
	op.add_option('-d','--maxDistance', action="store",type="float", dest="maxDistance",help="maximum distance (between closest atom pairs from different residues) that is considered a contact (in angstrom). default=%default")
	op.add_option('-i','--excludeRange', action="store",type="int", dest="excludeRange",help="number of residues up- and downstream of chain that are ignored as potential contacts. default=%default")
	op.add_option('-m','--model', action="store",type="int", dest="model",help="if PDB file contains ensemble of structures, specify model. default=%default")
	op.add_option('-r','--atomicSurrounding', action="store_true",dest="atomic",help="EXPERIMENTAL! examine surrounding atoms instead of residues")
	op.add_option('-x','--complex', action="store_true",dest="protComplex",help="EXPERIMENTAL! expecting two different chains of a complex for inter- AND intra-chain contacts.")
	op.add_option('-a','--atomContacts', action="store_true",dest="atomwiseContacts",help="EXPERIMENTAL! count contacts between atoms.")
	op.add_option('-w','--vdwOverlap', action="store",type="float",dest="VDWoverlap",help="EXPERIMENTAL! count contacts between atoms as van-der-Waals overlap. Chimera uses -0.4 angstrom as default for finding contacts.")
	op.add_option('-2','--twoCarbon', action="store_true",dest="twoCarbonCount",help="EXPERIMENTAL! count contacts between atoms, if both atoms are C atoms.")
	
	#DEFAULTS
	op.set_defaults(print_contact_list=False,count_total_contacts=False,print_summary=False,maxDistance=4,excludeRange=3,model=1,atomic=False,protComplex=False,atomwiseContacts=False,VDWoverlap=None,twoCarbonCount=False)
	opt, args = op.parse_args()
	
	# SET VARIABLES
	print_contact_list=opt.print_contact_list
	count_total_contacts=opt.count_total_contacts
	print_summary=opt.print_summary
	maxDistance=opt.maxDistance
	excludeRange=opt.excludeRange
	model=opt.model
	atomic=opt.atomic
	protComplex=opt.protComplex
	atomwiseContacts=opt.atomwiseContacts
	VDWoverlap = opt.VDWoverlap
	twoCarbonCount = opt.twoCarbonCount
	#print "VDWoverlap:",VDWoverlap
	
	
	
	# not enough command line arguments given
	if len(args) < 3:
		print "use -h for help"
		sys.exit(0)
	
	# command line arguments
	pdbfile=args[0]
	chain1=args[1]
	chain2=args[2]
	
	# for protein complex analysis we need two different chains (at least different by name, the chains can still be identical)
	if protComplex and chain1 == chain2:
		print "please specify two different chains of complex!"
		sys.exit(0)
	
	#####################################
	# START CALCULATIONS 
	#####################################
	
	result = hh_contact_count(pdbfile, chain1, chain2, maxDistance, excludeRange, model, print_contact_list, print_summary, count_total_contacts, atomic,protComplex,atomwiseContacts,VDWoverlap,twoCarbonCount)
	
	if not print_summary and not print_contact_list:
		print result
