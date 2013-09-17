#!/usr/bin/python
import cgi
import cgitb
import superCodonTools
import tempfile
import os
import super_nt_maker

################################################################################
def check_and_parse_input_args(form):
    percent = None
    multipleAlignmentData = []
    errs = []
    try:
        percent = int(form['Threshhold'].value)
    except:
        errs += ['Threshhold needs to be an integer value.']
    if(form['MultipleAlignment'].value != ''):
        multipleAlignmentData = form['MultipleAlignment'].value.split()
    if(any(len(item)!=len(multipleAlignmentData[0]) for item in multipleAlignmentData)):
        errs += ['All lines in the multiple alignment need to be of the same size.']
        
    if(errs != []):
        superCodonTools.print_errors(errs)
        import sys
        sys.exit(1)
        
    return(percent, multipleAlignmentData)

################################################################################
def create_distributions(percent, mad):
    distributions = []
    for i in range(len(mad[0])):
        currDistribution = {}
        for seq in mad:
            aa = seq[i].upper()
            if(aa in superCodonTools.DEFAULT_AA_ORDER):
                if(aa not in currDistribution):
                    currDistribution[aa] = 0
                currDistribution[aa] += 1
        total = float(sum(currDistribution.values()))
        currDistribution = dict([(aa, currDistribution[aa]/total) for aa in currDistribution])
        if(max(currDistribution.values()) <= percent):
            distributions += [(i+1,currDistribution)]
            
    return distributions
    
################################################################################
def write_distributions(ofn, distributions):
    of = open(ofn, 'w')
    
    for aa in superCodonTools.DEFAULT_AA_ORDER:
        of.write(',' + aa)
    of.write('\n')
    
    for dist in distributions:
        of.write('Position ' + str(dist[0]))
        for aa in superCodonTools.DEFAULT_AA_ORDER:
            if(aa in dist[1]):
                of.write(',' + str(dist[1][aa]))
            else:
                of.write(',0.0')
        of.write('\n')

    of.close()

################################################################################
if __name__ == '__main__':
    cgitb.enable()
    print 'Content-type: text/html\n\n'

    form = cgi.FieldStorage()
    (percent, mad) = check_and_parse_input_args(form)
    distributions = create_distributions(percent, mad)
    
    outputDir = tempfile.mkdtemp(prefix = superCodonTools.OUT_DIR_PREFIX, dir=superCodonTools.RESULTS_DIR)
    ofn = os.path.join(outputDir, superCodonTools.DIST_FILE)
    write_distributions(ofn, distributions)
    
    super_nt_maker.run(ofn, outputDir)
    ofn = os.path.join('../', 'results', os.path.split(outputDir)[-1], 'index.html')
    print superCodonTools.REDIRECT_CODE % (ofn, ofn, ofn)
    