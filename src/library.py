# p3d General biochemical definitions 
# encoding: utf-8
# Copyright (c) 2009 Ch. Fufezan.
# module is distributed under the terms of the GNU General Public License

import math

"""
Wrapper for

convert321(AA) - convert AA three letter code to one letter code
convert123(AA) - convert AA one letter code to three letter code
calcHydropathy(AA) - returns hydropathy value for AA
calcSideChainSurface(AA) - return sidechain surface in nm**2
calcMass(AA) - returns mass of AA


General References:

Library[key]

If not indicated otherwise, values hold dictionaries that have one letter AA as keys.

Keys:
	AA 				: Three leter Aminoacid codes in ['GLY', ...]
	AA_3_to_1		: three to one letter code
	AA_1_to_3		: one to three letter code
	Hydropathy		: Hydropathy index number (Hydrophilic negative)	{'A':1.8, ...}
	SidechainSurface: AA side chain surface volume in nm2				{'A':1.15, ...}
	AtomCount		: Number of Atoms									{'A':5, ...}
	Atoms			: Atom names in AA									{'A':['CA','C'...], ...}
	UniqueAtoms		: Unique atoms in given AA							{'CD':['GLN', 'LYS', 'PRO', 'GLU', 'ARG'], ...}
	Mass			: Mass in Da of given AA							{'A':89, ...}
	Non-Prot-Res	: Non-protein components, identified by Hic-up db	{'HEM':'Heme' , ....}
					  Obtained by:
					b = {}
					for entry in a:
					 qs = 'http://xray.bmc.uu.se/hicup/'+entry.strip()+'/index.html'
					 try:
					  url = urlopen(qs)
					  for line in url:
					   if line.startswith('<TITLE>'):
					    b[entry] = line.strip()
					 except:
					  print entry
					for k,v in b.items():
					 print "'"+k.strip()+"':'"+v[37:-15].strip()+"',"
					
					  And images by 
					for entry in Library['Non-Prot-Res'].keys():
					 qs = 'http://xray.bmc.uu.se/hicup/'+entry+'/'+entry.lower()+'_rcsb.gif'
					 try:
					  urllib.urlretrieve(qs, entry+'_rcsb.gif')
					 except:
					  print entry

	AA-CLASS		: AA classifications 								{'HYDROPHOBIC': ...}
	
					
"""

"""
This file is part of p3d.

    p3d is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    p3d is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""




AA3 = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'CYS', 'MET', 'PHE', 'TYR', 'TRP', 'PRO', 'SER', 'THR', 'ASN', 'GLN', 'ASP', 'GLU', 'HIS', 'LYS', 'ARG']
AA1 = ['G','A','V','L','I','C','M','F','Y','W','P','S','T','N','Q','D','E','H','K','R']

RAD_2_DEG = 180/math.pi


class InputError(Exception): pass
class NotAminoAcidError(InputError): pass
class Not3LetterAA(InputError): pass
class Not1LetterAA(InputError): pass

def convert321(xxx):
	'''
	converts aminoacid 3 letter code to 1
	'''
	if xxx in AA3:
		return Library['AA_3_to_1'][xxx]
	else:
		raise Not3LetterAA()

def convert123(x):
	'''
	converts aminoacid 1 letter code to 3
	'''
	if x in AA1:
		return Library['AA_1_to_3'][x]
	else:
		raise Not1LetterAA()

def calcHydropathy(xxx):
	if xxx in AA3:
		return Library['Hydropathy'][convert321(xxx)]
	elif xxx in AA1:
		return Library['Hydropathy'][xxx]
	else:
		raise NotAminoAcidError()	

def calcSideChainSurface(xxx):
	if xxx in AA3:
		return Library['SidechainSurface'][convert321(xxx)]
	elif xxx in AA1:
		return Library['SidechainSurface'][xxx]
	else:
		raise NotAminoAcidError()

def calcMass(xxx):
	if xxx in AA3:
		return Library['Mass'][convert321(xxx)]
	elif xxx in AA1:
		return Library['Mass'][xxx]
	else:
		raise NotAminoAcidError()

Library = {
'AA' : ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'CYS', 'MET', 'PHE', 'TYR', 'TRP', 'PRO', 'SER', 'THR', 'ASN', 'GLN', 'ASP', 'GLU', 'HIS', 'LYS', 'ARG'],
'AA_3_to_1' :{
					'ALA':'A',
					'ARG':'R',
					'ASN':'N',
					'ASP':'D',
					'CYS':'C',
					'GLN':'Q',
					'GLU':'E',
					'GLY':'G',
					'HIS':'H',
					'ILE':'I',
					'LEU':'L',
					'LYS':'K',
					'MET':'M',
					'PHE':'F',
					'PRO':'P',
					'SER':'S',
					'THR':'T',
					'TRP':'W',
					'TYR':'Y',
					'VAL':'V',
					'UNK':'X'
},
'AA_1_to_3' :{
					'G':'GLY',
					'A':'ALA',
					'V':'VAL',
					'L':'LEU',
					'I':'ILE',
					'C':'CYS',
					'M':'MET',
					'F':'PHE',
					'Y':'TYR',
					'W':'TRP',
					'P':'PRO',
					'S':'SER',
					'T':'THR',
					'N':'ASN',
					'Q':'GLN',
					'D':'ASP',
					'E':'GLU',
					'H':'HIS',
					'K':'LYS',
					'R':'ARG'
},
'pKa' :{
					'C':9.1,
					'Y':10.2,
					'D':4.75,
					'E':4.75,
					'H':6.98,
					'K':10.8,
					'R':12.5
},
'Hydropathy' :{
					# Hydropathy Index J.Kyte and R.F.Doolittle (1982) JMB 157, 105-132
					# Version 1.0 by Ch.Fufezan 2007
					# Mean hydropathy index is -0.49
					'A':1.8,
					'R':-4.5,
					'N':-3.5,
					'D':-3.5,
					'C':2.5,
					'Q':-3.5,
					'E':-3.5,
					'G':-0.4,
					'H':-3.2,
					'I':4.5,
					'L':3.8,
					'K':-3.9,
					'M':1.9,
					'F':2.8,
					'P':-1.6,
					'S':-0.8,
					'T':-0.7,
					'W':-0.9,
					'Y':-1.3,
					'V':4.2
},
'SidechainSurface' : {
					# Acce.surface in nm2 J.Kyte and R.F.Doolittle (1982) JMB 157, 105-132
					# Library by Ch.Fufezan 2007
					# Version 1.0
					# Mean surface size is 1.7025
					'A':1.15,
					'R':2.25,
					'N':1.6,
					'D':1.5,
					'C':1.35,
					'Q':1.8,
					'E':1.9,
					'G':0.75,
					'H':1.95,
					'I':1.75,
					'L':1.7,
					'K':2.0,
					'M':1.85,
					'F':2.1,
					'P':1.45,
					'S':1.15,
					'T':1.4,
					'W':2.55,
					'Y':2.3,
					'V':1.55
},
'AtomCount' :{
					'ALA':5,
					'ARG':11,
					'ASN':8,
					'ASP':8,
					'CYS':6,
					'GLN':9,
					'GLU':9,
					'GLY':4,
					'HIS':10,
					'ILE':8,
					'LEU':8,
					'LYS':9,
					'MET':8,
					'PHE':11,
					'PRO':7,
					'SER':6,
					'THR':7,
					'TRP':14,
					'TYR':12,
					'VAL':7
},
'Atoms' : {
					'ALA': ['CB', 'CA', 'C', 'O', 'N'],
					'ARG': ['C', 'CB', 'CA', 'CG', 'NE', 'O', 'N', 'CZ', 'NH1', 'NH2', 'CD'],
					'ASN': ['C', 'CB', 'CA', 'CG', 'O', 'N', 'OD1', 'ND2'],
					'ASP': ['C', 'CB', 'CA', 'CG', 'O', 'N', 'OD1', 'OD2'],
					'CYS': ['C', 'CB', 'CA', 'O', 'N', 'SG'],
					'GLN': ['C', 'CB', 'CA', 'CG', 'O', 'N', 'CD', 'NE2', 'OE1'],
					'GLU': ['C', 'CB', 'CA', 'CG', 'O', 'N', 'OE2', 'CD', 'OE1'],
					'GLY': ['CA', 'C', 'O', 'N'],
					'HIS': ['C', 'CE1', 'CB', 'CA', 'CG', 'O', 'N', 'CD2', 'ND1', 'NE2'],
					'ILE': ['C', 'CB', 'CA', 'O', 'N', 'CD1', 'CG1', 'CG2'],
					'LEU': ['C', 'CB', 'CA', 'CG', 'O', 'N', 'CD1', 'CD2'],
					'LYS': ['C', 'CB', 'CA', 'CG', 'O', 'N', 'NZ', 'CE', 'CD'],
					'MET': ['C', 'CB', 'CA', 'CG', 'O', 'N', 'CE', 'SD'],
					'PHE': ['C', 'CE1', 'CB', 'CA', 'CG', 'O', 'N', 'CZ', 'CD1', 'CD2', 'CE2'],
					'PRO': ['C', 'CB', 'CA', 'CG', 'O', 'N', 'CD'],
					'SER': ['C', 'OG', 'CB', 'CA', 'O', 'N'],
					'THR': ['C', 'CB', 'CA', 'OG1', 'O', 'N', 'CG2'],
					'TRP': ['C', 'CZ2', 'CB', 'CA', 'CG', 'CH2', 'O', 'N', 'CE2', 'CE3', 'CD1', 'CD2', 'CZ3', 'NE1'],
					'TYR': ['C', 'CE1', 'OH', 'CB', 'CA', 'CG', 'O', 'N', 'CZ', 'CD1', 'CD2', 'CE2'],
					'VAL': ['C', 'CB', 'CA', 'O', 'N', 'CG1', 'CG2']
},
'UniqueAtoms' : {
					'CD' : ['GLN', 'LYS', 'PRO', 'GLU', 'ARG'],
					'CD1': ['ILE', 'PHE', 'LEU', 'TRP', 'TYR'],
					'CD2': ['PHE', 'HIS', 'LEU', 'TRP', 'TYR'],
					'CE' : ['LYS', 'MET'],
					'CE1': ['PHE', 'HIS', 'TYR'],
					'CE2': ['PHE', 'TRP', 'TYR'],
					'CE3': ['TRP'],
					#'CG' : ['ASP', 'GLN', 'LYS', 'PRO', 'PHE', 'HIS', 'GLU', 'LEU', 'ARG', 'TRP', 'ASN', 'TYR', 'MET'],
					'CG1': ['ILE', 'VAL'],
					'CG2': ['ILE', 'THR', 'VAL'],
					'CH2': ['TRP'],
					'CZ' : ['PHE', 'ARG', 'TYR'],
					'CZ2': ['TRP'],
					'CZ3': ['TRP'],
					'ND1': ['HIS'],
					'ND2': ['ASN'],
					'NE' : ['ARG'],
					'NE1': ['TRP'],
					'NE2': ['GLN', 'HIS'],
					'NH1': ['ARG'],
					'NH2': ['ARG'],
					'NZ' : ['LYS'],
					'OD1': ['ASP', 'ASN'],
					'OD2': ['ASP'],
					'OE1': ['GLN', 'GLU'],
					'OE2': ['GLU'],
					'OG' : ['SER'],
					'OG1': ['THR'],
					'OH' : ['TYR'],
					'SD' : ['MET'],
					'SG' : ['CYS']
},
'Mass'			: {
'G':	57.02,
'A':	71.04,
'S':	87.03,
'P':	97.05,
'V':	99.07,
'T':	101.05,
'C':	103.01,
'L':	113.08,
'I':	113.08,
'N':	114.04,
'D':	115.03,
'Q':	128.06,
'K':	128.09,
'E':	129.04,
'M':	131.04,
'H':	137.06,
'F':	147.07,
'R':	156.1,
'Y':	163.06,
'W':	186.08

},
'Non-Prot-Res' : {
				'017':'(3r,3as,6ar)-hexahydrofuro(2,3-b)furan-3-yl(1s,2r)-3- (((4-aminophenyl)sulfonyl)(isobutyl)amino)-1-benzyl-2- hydroxypropylcarbamate; tmc114; uic-94017',
				'12P':'dodecaethylene glycol; polyethylene glycol peg400',
				'1EM':'(1s)-2-hydroxy-1-((nonanoyloxy)methyl)ethyl myristate',
				'2OB':'cholesteryl oleate; (3beta,9beta,14beta,17alpha)-cholest-5-en-3-yl (9z)- octadec-9-enoate',
				'2SA':'2-(9-(3,4-dihydroxy-5-phosphonooxymethyl-tetrahydro- furan-2-yl)-9h-purin-6-ylamino)-succinic acid; adenylosuccinic acid tetrahydro-furan-2-yl)-9h-purin-6-ylamino)- succinic acid',
				'3BI':'(2s)-2-(((r)-((3r)-3-carboxy-3-(((4-(((2,4- diaminopteridin-6-yl)methyl)(methyl)amino)phenyl) carbonyl)amino)propyl)(hydroxy) phosphoryl)methyl)pentanedioic acid',
				'3ON':'(3r)-3-hydroxy-8p-apocarotenol; (1r)-4-((1e,3e,5e,7z,9e,11z,13e,15e)-17-hydroxy-3,7, 12,16-tetramethylheptadeca-1,3,5,7,9,11,13,15-octaen- 1-yl)-3,5,5-trimethylcyclohex-3-en-1-ol',
				'4AD':'4-amino-1,4-dioxobutan-2-aminium adenosine-5p- monophosphate; asnamp',
				'4TA':'p1-(5p-adenosyl)p4-(5p-(2p-deoxy-thymidyl)) tetraphosphate',
				'6PL':'(4s,7r)-4-hydroxy-n,n,n-trimethyl-9-oxo-7- ((palmitoyloxy)methyl)-3,5,8-trioxa-4- phosphahexacosan-1-aminium 4-oxide; 1-palmitoyl-2-stearoyl-sn-glycero-3-phosphocholine',
				'6UL':'tetracosyl palmitate',
				'869':'(1-tert-butyl-5-hydroxy-1h-pyrazol-4-yl)(6- (dihydroxy(methyl)-lambda~4~-sulfanyl)-4p-methoxy-2- methyl-1,1p-biphenyl-3-yl)methanone; (1-tert-butyl-5-hydroxy-1h-pyrazol-4-yl)-(6- methanesul',
				'8DG':'8-oxo-2p-deoxyguanosine-5p-triphosphate',
				'A1T':'5-pentyl-n-((4p-(piperidin-1-ylcarbonyl)biphenyl-4- yl)methyl)-n-(1-(pyridin-2-ylmethyl)piperidin-4- yl)pyridine-2-carboxamide',
				'ACO':'acetyl coenzyme a',
				'ACP':'phosphomethylphosphonic acid adenylate ester; adenosine-5p-(beta, gamma-methylene)triphosphate',
				'ACR':'acarbose; 1,4-deoxy-4-((5-hydroxymethyl-2,3,4- trihydroxycyclohex-5,6-enyl)amino)fructose',
				'ADV':'alpha-beta methylene adp-ribose; ampcpr; ((5-(6-amino-purin-9-yl)-3,4-dihydroxy- tetrahydro-furan-2-ylmethoxy)-hydroxy- phosphorylmethyl)-phosphonic acid mono-(3,4,5- trihydroxy-tetrahydro-',
				'AKT':'10-decarboxymethylaclacinomycin t (dcmat); 10-(4-dimethylamino-5-hydroxy-6-methyl-tetrahydro- pyran-2-yloxy)-8-ethyl-1,8,11-trihydroxy-7,8,9,10- tetrahydro-naphthacene-5,12-dione',
				'ANP':'phosphoaminophosphonic acid-adenylate ester',
				'APC':'diphosphomethylphosphonic acid adenosyl ester; alpha,beta-methyleneadenosine-5p-triphosphate',
				'APX':'2p-monophosphoadenosine-5p-diphosphoribose',
				'ATG':'phosphothiophosphoric acid-adenylate ester',
				'ATI':'n-(3-amino-2-hydroxy-5-methylhexanoyl) valylvalylaspartic acid; amastatin',
				'ATP':'adenosine-5p-triphosphate',
				'AZX':'4-((3-chloro-4-(((2r)-3,3,3-trifluoro-2-hydroxy-2- methylpropanoyl)amino)phenyl)sulfonyl)-n,n- dimethylbenzamide',
				'B12':'cobalamin',
				'B13':'5-hydroxybenzimidazolylcob(iii)amide',
				'B4P':'bis(adenosine)-5p-tetraphosphate',
				'BA3':'bis(adenosine)-5p-triphosphate',
				'BCA':'4-hydroxybenzoyl coenzyme a',
				'BCD':'beta-cyclodextrin; cyclo-hepta-amylose',
				'BCL':'bacteriochlorophyll a',
				'BID':'bistramide a; (2s,3r)-3-hydroxy-n-(3-((2r,3s,6s,8s)-8-((3s,4e,6s)-6- hydroxy-3,5-dimethylhept-4-en-1-yl)-3-methyl-1,7- dioxaspiro(5.5)undec-2-yl)propyl)-2-methyl-4-((((2s, 3s,6r)-3-methyl-6',
				'BLA':'biliverdine ix alpha',
				'BPB':'bacteriopheophytin b',
				'BPH':'bacteriopheophytin a',
				'BT5':'biotinyl-5-amp',
				'C2E':'9,9p-((2r,3r,3as,5s,7ar,9r,10r,10as,12s,14ar)-3,5,10, 12-tetrahydroxy-5,12-dioxidooctahydro-2h,7h-difuro(3, 2-d:3p,2p-j)(1,3,7,9,2, 8)tetraoxadiphosphacyclododecine-2,9-diyl)bis(2-amino- 1,',
				'C2F':'5-methyl-5,6,7,8-tetrahydrofolic acid',
				'CAO':'oxidized coenzyme a',
				'CDL':'cardiolipin; diphosphatidyl glycerol; bis-(1,2-diacyl-sn-glycero-3- phospho)-1p,3p-sn-glycerol',
				'CDN':'cardiolipin',
				'CHL':'chlorophyll b',
				'CLA':'chlorophyll a',
				'CNC':'co-cyanocobalamin',
				'COA':'coenzyme a',
				'COH':'protoporphyrin ix containing co',
				'CPS':'3-((3-cholamidopropyl)dimethylammonio)-1- propanesulfonate; chaps',
				'CSF':'cytidine-5p-monophosphate-3-fluoro-n-acetyl- neuraminic acid; cmp-3fneuac',
				'CTO':'triacetylchitotriose',
				'CXR':'cyclic adenosine diphosphate-ribose; cyclic adp-ribose',
				'CXT':'carboxyatractyloside',
				'CZN':'(2s,8r)-8-benzyl-2-hydroperoxy-6-(4-hydroxyphenyl)-2- (2-naphthylmethyl)-7,8-dihydroimidazo(1,2-a)pyrazin- 3(2h)-one; n-coeleneterazine',
				'DAU':'2pdeoxy-thymidine-5p-diphospho-alpha-d-glucose',
				'DBV':'15,16-dihydrobiliverdin',
				'DGD':'digalactosyl diacyl glycerol (dgdg)',
				'DGT':'2p-deoxyguanosine-5p-triphosphate',
				'DHE':'heme d',
				'DMU':'decyl-beta-d-maltopyranoside; decylmaltoside',
				'DR9':'1-cis-9-octadecanoyl-2-cis-9-hexadecanoyl phosphatidyl glycerol; (2r)-3-(((((2s)-2,3-dihydroxypropyl)oxy)(hydroxy) phosphoryl)oxy)-2-((9e)-hexadec-9-enoyloxy)propyl (9e)-octadec-9-enoate',
				'DSU':'((2r,3s,4s,5s)-3,4-dihydroxy-5-(hydroxymethyl)-5-((2r, 3s,4s,5s,6r)-3,4,5-trihydroxy-6-methoxy-tetrahydro-2h- pyran-2-yloxy)-tetrahydrofuran-2-yl)methyl nonanoate',
				'DTA':'(2s,3s,4r,5r,2ps,3ps,4pr,5pr)-2,2p- (dithiobis(methylene))bis(5-(6-amino-9h-purin-9-yl) tetrahydrofuran-3,4-diol); di-(5p-thioadenosine)',
				'EPH':'l-alpha-phosphatidyl-beta-oleoyl-gamma-palmitoyl- phosphatidylethanolamine',
				'FA5':'adenosine-5p-(phenylalaninyl-phosphate)',
				'FAD':'flavin-adenine dinucleotide',
				'FMN':'flavin mononucleotide; riboflavin monophosphate',
				'FON':'folinic acid; 5-formyl-5,6,7,8-tetrahydrofolate',
				'FRG':'2-(3-methyl-4-(n-methyl-guanidino)-butyrylamino)-3-(4- phenylethynyl-phenyl)-propionic acid methyl ester; (r)-n-(2-(1-(aminoiminomethyl)-3-piperidinyl)-1- oxoethyl)-4-(phenylethynyl)-l-phen',
				'GDC':'guanosine-5p-diphosphate-beta-l-galactose',
				'GDD':'guanosine-5p-diphosphate-alpha-d-mannose',
				'GDX':'guanosine 5p-(trihydrogen diphosphate), pp-d- mannopyranosyl ester',
				'GEL':'1-o-octyl-2-heptylphosphonyl-sn-glycero-3- phosphoethanolamine',
				'GNP':'phosphoaminophosphonic acid-guanylate ester',
				'GTP':'guanosine-5p-triphosphate',
				'GUD':'glucose-uridine-c1,5p-diphosphate',
				'GVR':'(2r)-n-hydroxy-3-naphthalen-2-yl-2- ((naphthalen-2-ylsulfonyl)amino)propanamide',
				'H4M':'5,10-dimethylene tetrahydromethanopterin',
				'HAS':'heme-as',
				'HDC':'3r-hydroxydecanoyl-coenzyme a; 3r-hydroxydecanoyl-coa',
				'HEA':'heme-a',
				'HEC':'heme c',
				'HEM':'protoporphyrin ix containing fe; heme; (dihydrogen 3,7,12,17-tetramethyl-8,13- divinyl-2,18-porphinedipropionato(2-))iron',
				'HPI':'n-(1-carboxy-3-phenylpropyl)phenylalanyl-alpha- asparagine',
				'IB2':'p1-p2-methylene-p3-thio-diadenosine triphosphate; ado-p-ch2-p-ps-ado',
				'IHP':'inositol hexakisphosphate; myo-inositol hexakisphosphate',
				'KAN':'kanamycin a',
				'LAN':'lanosterol',
				'LBV':'3-((2z)-2-((3-(2-carboxyethyl)-5-((z)-((3e,4s)-3- ethylidene-4-methyl-5-oxopyrrolidin-2-ylidene)methyl)- 4-methyl-1h-pyrrol-2-yl)methylene)-4-methyl-5-((z)-(3- methyl-5-oxo-4-vinyl-1,5-dihy',
				'LHG':'1,2-dipalmitoyl-phosphatidyl-glycerole',
				'LI1':'1-(2,6,10.14-tetramethyl-hexadecan-16-yl)-2-(2,10,14- trimethylhexadecan-16-yl)glycerol; lipid fragment',
				'LMT':'dodecyl-beta-d-maltoside; lauryl-beta-d-maltoside',
				'LMU':'dodecyl-alpha-d-maltoside',
				'LP3':'(7s)-4,7-dihydroxy-n,n,n-trimethyl-10-oxo-3,5,9- trioxa-4-phosphaheptacosan-1-aminium 4-oxide',
				'LUX':'(3r,3pr,6ps,9r,9pr,13r,13ps)-4p,5p-didehydro-5p,6p,7p, 8p,9,9p,10,10p,11,11p,12,12p,13,13p,14,14p,15,15p- octadecahydro-beta,beta-carotene-3,3p-diol',
				'MCA':'methylmalonyl-coenzyme a',
				'MCN':'pterin cytosine dinucleotide',
				'MD1':'phosphoric acid 4-(2-amino-4-oxo-3,4,5,6,-tetrahydro- pteridin-6-yl)-2-hydroxy-3,4-dimercapto-but-3-en-yl ester guanylate ester',
				'MGD':'2-amino-5,6-dimercapto-7-methyl-3,7,8a,9-tetrahydro-8- oxa-1,3,9,10-tetraaza-anthracen-4-one guanosine dinucleotide; molybdopterin guanosine dinucleotide',
				'MGE':'(1s)-2-(alpha-l-allopyranosyloxy)-1-((tridecanoyloxy) methyl)ethyl palmitate; monogalactosyl-diacylglycerol',
				'MGP':'7-methyl-guanosine-5p-triphosphate',
				'MGT':'7n-methyl-8-hydroguanosine-5p-triphosphate',
				'MLR':'maltotriose; amylotriose',
				'MRR':'(r)-2-methylmyristoyl-coenzyme a; (5-(6-aminopurin-9-yl)-2-((((3-(2-(2-(r)-2- methyltetradecanoyl)-sulfanylethylcarbamoyl ethylcarbamoyl)-3-hydroxy-2,2-dimethyl-propoxy)- hydroxy-phosphoryl',
				'MYA':'tetradecanoyl-coa; myristoyl-coa',
				'NAD':'nicotinamide-adenine-dinucleotide',
				'NAI':'1,4-dihydronicotinamide adenine dinucleotide; nadh',
				'NAP':'nadp nicotinamide-adenine-dinucleotide phosphate; 2p-monophosphoadenosine 5p-diphosphoribose',
				'NDP':'nadph dihydro-nicotinamide-adenine-dinucleotide phosphate',
				'NEX':'(1r,3r)-6-((3e,5e,7e,9e,11e,13e,15e,17e)-18-((1s,4r, 6r)-4-hydroxy-2,2,6-trimethyl-7-oxabicyclo(4.1.0)hept- 1-yl)-3,7,12,16-tetramethyloctadeca-1,3,5,7,9,11,13, 15,17-nonaenylidene)-1,5,5-t',
				'OTP':'(2e,6e,10e,14e,18e,22e,26e)-3,7,11,15,19,23,27,31- octamethyldotriaconta-2,6,10,14,18,22,26,30-octaenyl trihydrogen diphosphate; octaprenyl pyrophosphate',
				'PC1':'1,2-diacyl-sn-glycero-3-phosphocholine; 3-sn-phosphatidylcholine',
				'PCW':'1,2-dioleoyl-sn-glycero-3-phosphocholine; (z,z)-4-hydroxy-n,n,n-trimethyl-10-oxo-7-((1-oxo-9- octadecenyl)oxy)-3,5,9-trioxa-4-phosphaheptacos-18-en- 1-aminium-4-oxide',
				'PE3':'3,6,9,12,15,18,21,24,27,30,33,36,39- tridecaoxahentetracontane-1,41-diol; polyethylene glycol',
				'PEB':'phycoerythrobilin',
				'PEF':'di-palmitoyl-3-sn-phosphatidylethanolamine; 3-(aminoethylphosphoryl)-(1,2-di-palmitoyl)-sn- glycerol',
				'PEH':'di-stearoyl-3-sn-phosphatidylethanolamine',
				'PEK':'(1s)-2-(((2-aminoethoxy)(hydroxy)phosphoryl)oxy)-1- ((stearoyloxy)methyl)ethyl (5e,8e,11e,14e)-icosa-5,8, 11,14-tetraenoate; phosphatidylethanolamine, 2-arachidonoyl-1-stearoyl- sn-glycerol',
				'PEU':'2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56, 59,62,65,68,71,74,77,80-heptacosaoxadooctacontan-82-ol; peg 8000',
				'PEV':'(1s)-2-(((2-aminoethoxy)(hydroxy)phosphoryl)oxy)-1- ((palmitoyloxy)methyl)ethyl stearate; phosphatidylethanolamine; 1-palmitoyl-2-oleoyl-sn- glycero-3-phosphoethanolamine',
				'PGM':'1-myristoyl-2-hydroxy-sn-glycero-3-(phospho-rac-(1- glycerol)); lysophosphatidylglycerol',
				'PGV':'(1r)-2-(((((2s)-2,3-dihydroxypropyl)oxy)(hydroxy) phosphoryl)oxy)-1-((palmitoyloxy)methyl)ethyl (11e)- octadec-11-enoate; phosphatidylglycerol, 2-vaccenoyl-1-palmitoyl-sn- glycerol-3-phosph',
				'PHO':'pheophytin a',
				'PIB':'2-(butanoyloxy)-1-(((hydroxy((2,3,4,6-tetrahydroxy-5- (phosphonooxy)cyclohexyl)oxy)phosphoryl) oxy)methyl)ethyl butanoate; d-myo-phosphatidylinositol 3-phosphated (+)-sn-1,2-di- o-butanoylg',
				'PID':'peridinin',
				'PMT':'phosphoric acid mono-(3-(3-((5-(4-amino-2-oxo-2h- pyrimidin-1-yl)-3,4- dihydroxy-tetrahydro-furan-2- ylmethoxy)-hydroxy-phosphoryloxy)-3-oxo- propylcarbamoyl)-3-hydroxy-2,2- dimethyl-propyl',
				'PQ9':'5-((2e,6e,10e,14e,18e,22e)-3,7,11,15,19,23,27- heptamethyloctacosa-2,6,10,14,18,22,26-heptaenyl)-2,3- dimethylbenzo-1,4-quinone',
				'PSC':'(7r,17e,20e)-4-hydroxy-n,n,n-trimethyl-9-oxo-7- ((palmitoyloxy)methyl)-3,5,8-trioxa-4-phosphahexacosa- 17,20-dien-1-aminium 4-oxide; phosphatidylcholine, 2-linoleoyl-1-palmitoyl-sn- gycerol',
				'PTY':'phosphatidylethanolamine',
				'PUL':'(1s,2s,3e,5e,7e,10s,11s,12s)-12-((2r,4e,6e,8z,10r,12e, 14e,16z,18s,19z)-10,18-dihydroxy-12,16,19-trimethyl- 11,22-dioxooxacyclodocosa-4,6,8,12,14,16,19-heptaen-2- yl)-2,11-dihydroxy-1,10-di',
				'QKH':'methyl-cyclo-hepta-amylose; methyl-beta-cyclodextrin; (1s,3r,5r,6s,8r,10r,11s,13r, 15r,16s,18r,20r,21r,23r,25r,26r,28r,30r,31s,33r,35r, 36r,37r,38r,39r,40r,41r,42s,43r,44s,45r,46s,47r,48s,',
				'QSI':'5p-o-(n-(l-glutaminyl)-sulfamoyl)adenosine',
				'RFP':'rifampicin',
				'RG1':'rhodopin glucoside',
				'RPL':'(c8-s)-hydantocidin 5p-phosphate; (8,9-dihydroxy-3-(4-carboxy-hydroxy-hydroxymethyl- amino-butyl)-2,4-dioxo-6-oxa-1,3-diaza-spiro(4.4)non- 7-ylmethyl) phosphate',
				'SAP':'adenosine-5p-diphosphate monothiophosphate',
				'SLT':'5-(acetylamino)-3,5-dideoxynon-2-ulopyranonosyl-(2- >3)-beta-d-lyxo-hexopyranosyl-(1->4)hexopyranose; lactose sialic acid; alpha(2,3) sialyl lactose',
				'SLU':'5p-o-(n-(dehydroluciferyl)-sulfamoyl) adenosine',
				'SMA':'stigmatellin a',
				'SPO':'spheroidene',
				'SRM':'siroheme',
				'T5X':'2-c-(3-((4-amino-2-methylpyrimidin-5-yl)methyl)-5-(2- (((r)-hydroxy(phosphonooxy)phosphoryl)oxy)ethyl)-4- methyl-1,3-thiazol-3-ium-2-yl)-5-o-phosphono-d-xylitol; d-xylulose-5-phosphate thia',
				'TDK':'3-((4-amino-2-methylpyrimidin-5-yl)methyl)-2-((1s)-1- hydroxy-1-((r)-hydroxy(methoxy)phosphoryl)ethyl)-5-(2- (((s)-hydroxy(phosphonooxy)phosphoryl)oxy)ethyl)-4- methyl-1,3-thiazol-3-ium; 2-',
				'TGL':'tristearoylglycerol; triacylglycerol',
				'TP6':'3-(1,3,7-trihydro-9-d-ribityl-2,6,8-purinetrione-7- yl) 1-phosphate',
				'TPG':'2,2,7-trimethyl-guanosine-5p-triphosphate-5p-guanosine',
				'U10':'ubiquinone-10',
				'U20':'uridine-5p-diphosphate-3-o-(r-3-hydroxymyristoyl)-n- acetyl-d-glucosamine; (2r,3r,4r,5s,6r)-3-(acetylamino)-2-(((r)-(((s)-(((2r, 3s,4r,5r)-5-(2,4-dioxo-3,4-dihydropyrimidin-1(2h)-yl)- 3,4-d',
				'U2F':'uridine-5p-diphosphate-2-deoxy-2-fluoro-alpha-d- glucose alpha-d-glucose',
				'UD1':'uridine-diphosphate-n-acetylglucosamine',
				'UD2':'uridine-diphosphate-n-acetylgalactosamine',
				'UMA':'uridine-5p-diphosphate-n-acetylmuramoyl-l-alanine',
				'UMQ':'undecyl-maltoside; undecyl-beta-d-maltopyranoside',
				'UPG':'uridine-5p-diphosphate-glucose; uridine-5p-monophosphate glucopyranosyl-monophosphate ester',
				'UPP':'phenyl-uridine-5p-diphosphate',
				'VIA':'5-(2-ethoxy-5-((4-methylpiperazin-1-yl) sulfonyl)phenyl)-1-methyl-3-propyl-1h,6h,7h- pyrazolo(4,3-d)pyrimidin-7-one; sildenafil, viagra',
				'VIV':'(2r)-2,5,7,8-tetramethyl-2-((4r,8r)-4,8,12- trimethyltridecyl)chroman-6-ol',
				'VS1':'3-((n-(morpholin-n-yl)-carbonyl)-phenylalaninyl- amino)-5- phenyl-pentane-1-sulfonylbenzene',
				'XAT':'(3s,5r,6s,3ps,5pr,6ps)-5,6,5p,6p-diepoxy-5,6,5p,6p- tetrahydro-beta,beta-carotene-3,3p-diol; violaxanthin',
				'XPE':'3,6,9,12,15,18,21,24,27-nonaoxanonacosane-1,29-diol; decaethylene glycol',
				'XTG':'4-nitrophenyl 6-thio-6-s-alpha-d-xylopyranosyl-beta-d- glucopyranoside; 4-nitrophenyl-(6-s-alpha-d-xylopyranosyl)-beta-d-',
				'YOK':'((2,2p-(4-carboxyethyl-1,2- phenylenebis(nitrilomethylidyne))bis(phenolato))(2-)- n,np,o,op)-iron; salophen-10-propionate iron chelate',
				'ZBA':'12,13-epoxytrichothec-9-ene-3,4,8,15-tetrol-4,15- diacetate-8-isovalerate; none',
				'ZIO':'(3s,5r,6s,7r,8r,11r,12s,13r,14s,15s)-6- hydroxy-5,7,8,11,13,15-hexamethyl-4,10-dioxo-14- ((3,4,6-trideoxy-3-(dimethylamino)-beta-d-xylo- hexopyranosyl)oxy)-1,9-dioxaspiro(2.13) hexadec-12-y',
				},
'AA_CLASS' : {
				'HYDROPHOBIC'      : ['G','A','V','L','I','P'],
				'AROMATIC'         : ['F','Y','W'],
				'HYDROPHILIC_0CRG' : ['S','T','C','M','N','Q'],
				'HYDROPHILIC_-CRG' : ['D','E'],
				'HYDROPHILIC_+CRG' : ['K','R','H']
			}

}

if __name__ == "__main__":
    print("yes")