#!/usr/bin/python
import cgi
import cgitb
import superCodonTools
import sys
import super_nt_maker
import tempfile
import os

ERR_TOLERANCE = 0.005

################################################################################
def check_and_parse_input_args(form):
    distributions = form['DesiredAADists'].value.split('\n')
    aas = []
    errs = []
    dists = []
    
    # Allow for a starting comma
    aas = distributions[0].strip().split(',')
    if(aas[0] == []):
        aas = aas[1:]
    # Check that these are just standard AA's    
    for aa in aas:
        aa = aa.upper()
        if(aa not in superCodonTools.DEFAULT_AA_ORDER):
            if(aa not in superCodonTools.THREE_LETTER_AA_DICTIONARY):
                errs += ['Unrecognized amino acid: ' + aa]
            else:
                aa = superCodonTools.THREE_LETTER_AA_DICTIONARY[aa]
        
        aas += [aa]
    
    for i, dist in enumerate(distributions[1:]):
        currVals = []
        
        # Allow for the first element to be the name of this distribution.
        vals = dist.split(',')
        if(len(vals) > len(aas)):
            name = vals[0]
            vals = vals[1:]
        else:
            name = 'Distribution ' + str(i)
        
        for v in vals:
            try:
                currVals += [float(v)]
            except:
                errs += ['Error reading value ' + v]
        if(len(currVals) != len(aas)):
            err = 'Number of AA percentages does not match number AA\'s.\n'
            err += 'Number of AA\'s %d, # of AA percentages %d.' % (len(aas), len(currVals))
            errs += [err]
        if(abs(1-sum(currVals)) > ERR_TOLERANCE):
            errs += ['Percentages should sum to 1.  Instead, they sum to: %f' % (sum(currVals))]
        
        dists += [(name, currVals)]
    
    if(errs != []):
        superCodonTools.print_errors(errs)
        sys.exit(1)
    
    return(aas, dists)

################################################################################
def write_distributions(ofn, aas, distributions):
    of = open(ofn, 'w')
    for aa in aas:
        of.write(',' + aa)
    of.write('\n')
    
    for i, dist in enumerate(distributions):
        if(len(dist) == 0):
            continue
        of.write(dist[0])
        for v in dist[1]:
            of.write(',' + str(v))
        of.write('\n')
    
    of.close()

################################################################################
if __name__ == '__main__':
    cgitb.enable()
    print 'Content-type: text/html\n\n'
    
    form = cgi.FieldStorage()
    (aas, distributions) = check_and_parse_input_args(form)
    
    outputDir = tempfile.mkdtemp(prefix = superCodonTools.OUT_DIR_PREFIX, dir=superCodonTools.RESULTS_DIR)
    ofn = os.path.join(outputDir, superCodonTools.DIST_FILE)
    write_distributions(ofn, aas, distributions)
    
    super_nt_maker.run(ofn, outputDir, objectiveFunction='5')
    ofn = os.path.join('../', 'results', os.path.split(outputDir)[-1], 'index.html')
    print superCodonTools.REDIRECT_CODE % (ofn, ofn, ofn)
    
