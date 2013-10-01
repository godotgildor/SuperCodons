#!/usr/bin/python
import urllib
import cgi
import cgitb
import superCodonTools
import re
import super_nt_maker
import os
import tempfile
import cPickle

AB_PICKLE_FILE = './ab_dists.pickle'
ABYSIS_CREDIT = 'Antibody distributions from <a href="http://www.bioinf.org.uk/abysis/">abYsis</a>'

################################################################################
def check_and_parse_input_args(form, distributions):
    abNumbering = form['ab_numbering'].value
    abPositions = form['ab_positions'].value.split(',')
    species = form['species'].value
    errs = []
    positions = []
    
    key = (species, abNumbering)
    
    for position in abPositions:
        position = position.strip().upper()
        if(position not in distributions[key]):
            errs += [position + ' not a valid ' + abNumbering + ' position.']
        positions += [position]
        
    if(errs != []):
        superCodonTools.print_errors(errs)
        import sys
        sys.exit(1)
        
    return(abNumbering, positions, species)

################################################################################
def write_distributions(ofn, distributions):
    of = open(ofn, 'w')
    
    for aa in superCodonTools.DEFAULT_AA_ORDER[0:-1]:
        of.write(',' + aa)
    of.write('\n')
    for dist in distributions:
        of.write(dist[0])
        for aa in superCodonTools.DEFAULT_AA_ORDER[0:-1]:
            of.write(',' + str(dist[1][aa]))
        of.write('\n')
    
    of.close()

################################################################################
if __name__ == '__main__':
    cgitb.enable()
    print 'Content-type: text/html\n\n'
    
    form = cgi.FieldStorage()
    
    distributions = cPickle.load(open(AB_PICKLE_FILE, 'rb'))
    (abNumbering, positions, species) = check_and_parse_input_args(form, distributions)
    #abNumbering = 'Kabat'
    #positions = ['L78','L76','L77','L80']
    #species = '32'
    
    dists = []
    key = (species, abNumbering)
    for position in positions:
        dists += [(position, distributions[key][position])]
        
    outputDir = tempfile.mkdtemp(prefix = superCodonTools.OUT_DIR_PREFIX, dir=superCodonTools.RESULTS_DIR)
    ofn = os.path.join(outputDir, superCodonTools.DIST_FILE)
    write_distributions(ofn, dists)

    super_nt_maker.run(ofn, outputDir, objectiveFunction='4', numThreads=7, aaLimitsStr='Stop:0.1', extraComment=ABYSIS_CREDIT)
    ofn = os.path.join('../', 'results', os.path.split(outputDir)[-1], 'index.html')
    print superCodonTools.REDIRECT_CODE % (ofn, ofn, ofn)
