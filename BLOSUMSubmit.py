#!/usr/bin/python
import cgi
import cgitb
import sys
import superCodonTools
import tempfile
import super_nt_maker
import cPickle
import numpy
import os

################################################################################
def check_and_parse_input_args(form):
    aminoAcids = form['AminoAcids'].value.split(',')
    identityPercent = form['IdentityPercent'].value.strip('%')
    aas = []
    errs = []
    
    if( len(aminoAcids) <= 0 ):
        errs += ['Need to specify at least one amino acid.']
    
    for aa in aminoAcids:
        aa = aa.upper()
        if aa not in superCodonTools.DEFAULT_AA_ORDER:
            if(aa not in superCodonTools.THREE_LETTER_AA_DICTIONARY):
                errs += ['Unrecognized amino acid: ' + aa]
            else:
                aas += [superCodonTools.THREE_LETTER_AA_DICTIONARY[aa]]
        else:
            aas += [aa]
            
    try:
        identityPercent = int(identityPercent)
    except:
        errs += ['Identity Percent needs to be an integer.']    

    if(errs != []):
        superCodonTools.print_errors(errs)
        sys.exit(1)
    
    return (identityPercent, aas)

################################################################################
def get_scaled_blosum_distributions(identityPercent):
    identityPercent = float(identityPercent) / 100.0
    of = open(superCodonTools.BLOSUM_DESIRED_DISTRIBUTIONS_PICKLE, 'rb')
    desiredDists = cPickle.load(of)

    # Load the BLOSUM distributions and scale them so that the 
    # identity amino acid has a value of identityPercent.
    for aa in superCodonTools.DEFAULT_AA_ORDER[0:-1]:
        desiredDist = desiredDists[aa][1]
        desiredDist2 = [desiredDists[aa][0], []]
        for aa2 in superCodonTools.DEFAULT_AA_ORDER[0:-1]:
            # If this is the identity amino acid, i.e. we are an ALA
            # and will mutate to an ALA, then set to identityPercent
            if(aa2 == aa):
                desiredDist2[1] += [identityPercent]
            else:
                # Scale the desired frequency of the non-identity amino acids
                # relative to the total non-identity frequency.
                desiredDist2[1] += [desiredDist[aa2]/(1-desiredDist[aa]) * (1-identityPercent)]
        # Normalize to a total freq of 1.  It should be close already, but may be 
        # floating point errors.
        desiredDist2[1] = list(numpy.array(desiredDist2[1]) / sum(desiredDist2[1]))
        desiredDists[aa] = desiredDist2

    return desiredDists    

################################################################################
def write_distributions(ofn, aas, identityPercent=0.5):
    of = open(ofn, 'w')

    desiredDistributions = get_scaled_blosum_distributions(identityPercent)
    
    for aa in superCodonTools.DEFAULT_AA_ORDER[0:-1]:
        of.write(',' + aa)
    of.write('\n')
    for aa in aas:
        of.write(desiredDistributions[aa][0])
        for i in range(len(superCodonTools.DEFAULT_AA_ORDER[0:-1])):
            of.write(',' + str(desiredDistributions[aa][1][i]))
        of.write('\n')
    
    of.close()

################################################################################
if __name__ == '__main__':
    cgitb.enable()
    print 'Content-type: text/html\n\n'
    
    form = cgi.FieldStorage()
    (identityPercent, aas) = check_and_parse_input_args(form)
    
    outputDir = tempfile.mkdtemp(prefix = superCodonTools.OUT_DIR_PREFIX, dir=superCodonTools.RESULTS_DIR)
    ofn = os.path.join(outputDir, superCodonTools.DIST_FILE)
    write_distributions(ofn, aas, identityPercent)
    
    super_nt_maker.run(ofn, outputDir)
    ofn = os.path.join('../', 'results', os.path.split(outputDir)[-1], 'index.html')
    print superCodonTools.REDIRECT_CODE % (ofn, ofn, ofn)
    
