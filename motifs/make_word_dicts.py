import sys
import os
import miscUtils
import openbabel

################################################################################
def get_bonds(obConversion, f):
    bonds = []
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, f)
    
    for i in range(mol.NumBonds()):
        bond = mol.GetBond(i)
        res = bond.GetBeginAtom().GetResidue()
        res2 = bond.GetEndAtom().GetResidue()
        bondAtoms = [res.GetAtomID(bond.GetBeginAtom()).strip(), res2.GetAtomID(bond.GetEndAtom()).strip()]
        bondAtoms.sort()
        bonds += [bondAtoms]
    
    return bonds

################################################################################
if(len(sys.argv) < 2):
    print 'Usage: python make_word_dicts.py <directory>'
    sys.exit(1)

fl = miscUtils.getFileList(sys.argv[1], '*.pdb')

obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats('pdb', 'pdb')

wordsDict = {}
bondWordDict = {}

for i, f in enumerate(fl):
    print 'Working on ' + f + ' ' + str(i+1) + ' of ' + str(len(fl)) + '.'
    
    (dir, fn) = os.path.split(f)
    word = os.path.split(dir)[-1]
    resn = os.path.splitext(fn)[0]
    bonds = get_bonds(obConversion, f)
    if(word not in wordsDict):
        wordsDict[word] = {}
    for bond in bonds:
        key = (resn, bond[0], bond[1])
        if(resn not in wordsDict[word]):
            wordsDict[word][resn] = '+'.join([bond[0], bond[1]])
        else:
            wordsDict[word][resn] += '+' + '+'.join([bond[0], bond[1]])
        bondWordDict[key] = word
    
open('wordsDict.py', 'w').write('WORDS_DICT = ' + str(wordsDict) + '\n\nBOND_WORD_DICT = ' + str(bondWordDict) + '\n\n')
