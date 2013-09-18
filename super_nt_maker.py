#!/usr/bin/python
import sys
import os
import numpy
import miscUtils
import scipy.optimize
import scipy.stats
import multiprocessing
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import superCodonTools

OPTIONS = {'desired_distributions' :   (True, 'A file which gives the desired destribution of amino acids.'),
           'objective_function':       (False, 'The objective function to minimize.', '4'),
           'num_threads':              (False, 'The number of threads to launch.', 4),
           'output_dir':               (True,  'The name of a dir to place the output html and images.'),
           'aaLimits':                 (False,  'Any limits on the percentage of certain AA\'s?','Stop:0.1'),
           'num_super_nts':            (False,  'The number of Super Nucleotides to use.', 4)
           }
                      
NUCLEOTIDES = ['A', 'C', 'T', 'G']

MAX_NUM_THREADS = 128

CODON_TABLE = { 'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Stop', 'TAG': 'Stop',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'TGT': 'C', 'TGC': 'C', 'TGA': 'Stop', 'TGG': 'W',
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
                
def calc_kl_divergence(desiredDist, actualDist):
    p = desiredDist[0:21]
    q = actualDist[0:21]
    p[p <= 0.0] = numpy.finfo('f').tiny
    q[q <= 0.0] = numpy.finfo('f').tiny
    kl_div = 0.5*numpy.sum(p*(numpy.log2(p) - numpy.log2(q)) + q*(numpy.log2(q) - numpy.log2(p)))
    
    return kl_div


def calc_wei_wang_obj_fun_4(desiredDist, actualDist):
    p = desiredDist[0:21]
    q = actualDist[0:21]
    
    val = numpy.sum(numpy.subtract(1, numpy.cos((p-q)*numpy.pi)))
    
    return val


def calc_wei_wang_obj_fun_5(desiredDist, actualDist):
    p = desiredDist[0:21]
    q = actualDist[0:21]
    p = p + numpy.finfo('f').tiny
    q = q + numpy.finfo('f').tiny
    val = numpy.sum(q*(numpy.log(q) - numpy.log(p)) + 0.5*(p-q)*(p-q))
    
    return val

def calc_bth_obj_fun_5(desiredDist, actualDist):
    p = desiredDist[0:21]
    q = actualDist[0:21]
    p = p + numpy.finfo('f').tiny
    q = q + numpy.finfo('f').tiny
    maxIndex = (desiredDist==max(desiredDist))
    maxDiff = abs(desiredDist[maxIndex] - actualDist[maxIndex])
    val = numpy.sum(p*(numpy.log(p) - numpy.log(q))) + maxDiff
    
    return val

def calc_relative_entropy(desiredDist, actualDist):
    scaleFactor = 1.0/(12.0*numpy.log(10))
    eps = 10**-6
    relEntropy = scaleFactor*numpy.sum((desiredDist - actualDist) * numpy.log((desiredDist + eps)/(actualDist+eps)))
    return relEntropy

def get_aa_distribution(ntVats):
    '''This function will go through all the nucleotides in each vat
    and calculate the probability of each amino acid.'''
    aaDist = {}
    for i0, p0 in enumerate(ntVats[0]):
        for i1, p1 in enumerate(ntVats[1]):
            for i2, p2 in enumerate(ntVats[2]):
                codon = NUCLEOTIDES[i0] + NUCLEOTIDES[i1] + NUCLEOTIDES[i2]
                aa = CODON_TABLE[codon]
                if(aa not in aaDist):
                    aaDist[aa] = 0.0
                    
                percentage = p0*p1*p2
                aaDist[aa] += percentage
                
    return numpy.array([aaDist[aa] for aa in superCodonTools.DEFAULT_AA_ORDER])

def get_aa_distributions(ntVats):
    '''For the given set of super-nucleotides, calculate the possible amino 
    acid distributions.'''
    aaDists = {}
    for i, vat1 in enumerate(ntVats):
        for j, vat2 in enumerate(ntVats):
            for k, vat3 in enumerate(ntVats):
                aaDists[(i, j, k)] = get_aa_distribution([vat1, vat2, vat3])
    
    return aaDists

def prune_dists_for_cutoff(aaDists, aaLimits):
    '''Remove any distributions that exceed our given maximum stop percentage.'''
    aaDistsPruned = {}
    for superNtNumber in aaDists:
        aaDist = aaDists[superNtNumber]
        # Loop over all of the limits that were set and make sure
        # this distribution satisfies them all.
        passesLimits = True
        for limit in aaLimits:
            if limit in superCodonTools.AA_ORDER_DICT:
                if(aaDist[superCodonTools.AA_ORDER_DICT[limit]] > aaLimits[limit]):
                    passesLimits = False
                    break
        if(passesLimits):
            aaDistsPruned[superNtNumber] = aaDist
        
    return aaDistsPruned
    
def bound_check(superNt):
    '''The bounds setting on the minimizer seem to be a little soft.
       Make sure all probs are between 0 and 1.'''
    modified = False
    for i, nt in enumerate(superNt):
        if(nt < 0.0):
            superNt[i] = 0.0
            modified = True
        if(nt > 1.0):
            superNt[i] = 1.0
            modified = True
    
    # If we modified a value, we'll renormalize.  This is
    # probably over-kill, but what the hey.
    if(modified):
        total = sum(superNt)
        for i, nt in enumerate(superNt):
            superNt[i] = nt/total
    
    return superNt
    
    
def unpack_vat_info(currVats):
    '''This function will a) unpack data into actual vats.
                          b) make sure that all vats have % > 0 and < 1'''
    superNts = []
    # Unpack the values into super-nucleotides.
    for i in range(0,len(currVats), 3):
        superNt = [currVats[i], currVats[i+1], currVats[i+2], 1 - sum(currVats[i:i+3])]
        superNt = bound_check(superNt)
        superNts += [superNt]
        
    return superNts

def find_best_distribution_match(desiredDistribution, aaDists, objFun):
    '''This function will look for the best distribution in aaDists to match to
    the desiredDistribution.'''    
    bestMatch = {}
    bestMatch['score'] = sys.float_info.max
    bestMatch['desiredDistribution'] = desiredDistribution
    bestMatch['calcDistribution'] = None
    bestMatch['superNtNumber'] = None
    
    for superNtNumber in aaDists:
        aaDist = aaDists[superNtNumber]
        score = objFun(desiredDistribution, aaDist)
        if(score < bestMatch['score']):
            bestMatch['score'] = score
            bestMatch['calcDistribution'] = aaDist
            bestMatch['superNtNumber'] = superNtNumber
    
    return bestMatch

def find_best_dist_fits(currVats, desiredDistributions, aaLimits, objFun, simple):
    '''This is the full objective function to be minimized. It basically
    will calculate the objective function on each distribution, and then
    sum those results for an uber-objective function value.'''
    superNts = unpack_vat_info(currVats)
          
    # Get the distributions that these super-nucleotides can produce.
    aaDists = get_aa_distributions(superNts)
    
    # Now remove any distributions that exceed our maximum stop codon percentage.
    if(aaLimits != {}):
        aaDists = prune_dists_for_cutoff(aaDists, aaLimits)
    
    
    # For each of our desired distributions, find the distribution from
    # these super-nucleotides that best matches per our objective function,
    # and keep track of a total score.
    retVal = {}
    retVal['total_score'] = 0.0
    retVal['matches'] = []
    retVal['superNts'] = superNts
    for title in desiredDistributions:
        desiredDistribution = desiredDistributions[title]
        bestMatch = find_best_distribution_match(desiredDistribution, aaDists, objFun)
        bestMatch['title'] = title
        bestMatch['stats'] = {'corr': scipy.stats.pearsonr(desiredDistribution, bestMatch['calcDistribution'])[0], 
                              'relEntropy': calc_relative_entropy(desiredDistribution, bestMatch['calcDistribution'])}
        
        retVal['matches'] += [bestMatch]
        retVal['total_score'] += bestMatch['score']
    
    if(simple):
        return retVal['total_score']
    
    return retVal

def load_desired_distributions(distFile):
    '''This function will read in a set of desired distributions.  The first
    line should give a list of the amino acid order, and subsequent lines
    should give the distributions.  The desired Stop percentage is given
    a default of 0.0 if not specified.'''
    aas = None
    dists = {}
    for l in open(distFile):
        w = l.strip().split(',')
        if(aas == None):
            aas = w[1:]
        else:
            title = w[0]
            currDist = [float(f) for f in w[1:]]
            currDist = dict(zip(aas, currDist))
            if 'Stop' not in currDist:
                currDist['Stop'] = 0.0
            dists[title] = numpy.array([currDist[aa] for aa in superCodonTools.DEFAULT_AA_ORDER])
    
    return dists

def get_random_vats(numSuperNucleotides):
    '''This function will simply create four sets of super-nucleotides, each
    composed of a random percentage of A, C, T, and G.'''
    vats = []
    for i in range(numSuperNucleotides):
        currVat = numpy.random.random(4)
        currVat = currVat / sum(currVat)
        vats += list(currVat[0:3])

    return vats

def get_best_result(results):
    bestResult = {}
    bestResult['score'] = sys.float_info.max
    
    for result in results:
        if(result.fun < bestResult['score']):
            bestResult['score'] = result.fun
            bestResult['val'] = list(result.x)
    
    # Round to the nearest percent
    retVal = [round(a,2) for a in bestResult['val']]
    
    return retVal

def get_obj_fun(objFunName):
    indiv_obj_fun = None
    
    if(objFunName == '4'):
        indiv_obj_fun = lambda x, y: calc_wei_wang_obj_fun_4(x,y)
    elif(objFunName == '5'):
        indiv_obj_fun = lambda x, y: calc_wei_wang_obj_fun_5(x,y)
    elif(objFunName == 'bth'):
        indiv_obj_fun = lambda x, y: calc_bth_obj_fun_5(x,y)
    else:
        indiv_obj_fun = lambda x, y: calc_kl_divergence(x,y)
    
    return indiv_obj_fun

def optimize_supernts(desiredDistributions, aaLimits, whichObjFun, initial_guess):
    myBounds = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
                (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
                (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
                (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),]
                
    myConstraints = ({'type': 'ineq',
                      'fun' : lambda x: numpy.array([1 - (x[0] + x[1] + x[2])])},
                     {'type': 'ineq',
                      'fun' : lambda x: numpy.array([1 - (x[3] + x[4] + x[5])])},
                     {'type': 'ineq',
                      'fun' : lambda x: numpy.array([1 - (x[6] + x[7] + x[8])])},
                     {'type': 'ineq',
                      'fun' : lambda x: numpy.array([1 - (x[9] + x[10] + x[11])])},
                     )

    indiv_obj_fun = get_obj_fun(whichObjFun)
                     
    objFun = lambda x: find_best_dist_fits(x, desiredDistributions, 
                                           aaLimits, indiv_obj_fun, True)
    result = scipy.optimize.minimize(objFun, 
                                     initial_guess, 
                                     method='SLSQP',
                                     bounds = myBounds, 
                                     constraints = myConstraints)
    return result

def make_plots(bestDistributionFits, outputDir):
    for dists in bestDistributionFits['matches']:
        calcDist = 100.0*dists['calcDistribution']
        # Note that the mismatch in the deisred vs. calc below is on purpose
        # so that we make sure that both dists use the same list of aa in same order.
        desiredDist = 100.0*dists['desiredDistribution']
        
        barWidth = 0.35
        xLocation = numpy.arange(len(calcDist))
        fig = pyplot.figure()
        ax = fig.add_subplot(1,1,1)
        rects1 = ax.bar(xLocation, calcDist, barWidth, color='r')
        rects2 = ax.bar(xLocation+barWidth, desiredDist, barWidth, color='y')
        ax.legend((rects1[0], rects2[0]), ('Actual Distribution', 'Desired Distribution'))
        ax.set_ylabel('Percentage')
        ax.set_xticks(xLocation+barWidth)
        ax.set_xticklabels(superCodonTools.DEFAULT_AA_ORDER)
        fig.suptitle(dists['title'])
        ax.set_title(dists['superNtNumber'])
        pyplot.savefig(os.path.join(outputDir, dists['title'] + '.png'), dpi=100)
        
################################################################################
def get_initial_vats(num_attempts, numSuperNucleotides):
    for i in range(num_attempts):
        initial_vats += [get_random_vats(numSuperNucleotides)]
    
    return initial_vats
    

################################################################################
def print_results(bestResult, bestDistributionFits, outputDir):
    ofn = os.path.join(outputDir, 'index.html')
    of = open(ofn, 'w')
    of.write(superCodonTools.HTML_PREFIX + '\n')
    of.write('      <div class="row">\n')
    of.write('        <h2>Super Nucleotides</h2><br/>\n')
    of.write('        <table class="table table-bordered">\n')
    of.write('          <thead>\n')
    of.write('            <tr>\n')
    of.write('              <td>Super NT #</td>\n')
    for nt in NUCLEOTIDES:
        of.write('              <td>' + nt + '</td>\n')
    of.write('            </tr>\n')
    of.write('          </thead>\n')
    of.write('          <tbody>\n')
    of.write('            <tr class="success">\n')
    of.write('              <td>0</td>\n')
    total = 0
    for v in bestResult[0:3]:
        v = int(100*v)
        total += v
        of.write('              <td>' + str(v) + '%</td>\n')
    of.write('              <td>' + str(100-total) + '%</td>\n')
    of.write('            </tr>\n')
    of.write('            <tr class="error">\n')
    of.write('              <td>1</td>\n')
    total = 0
    for v in bestResult[3:6]:
        v = int(100*v)
        total += v
        of.write('              <td>' + str(v) + '%</td>\n')
    of.write('              <td>' + str(100-total) + '%</td>\n')
    of.write('            </tr>\n')
    of.write('            <tr class="warning">\n')
    of.write('              <td>2</td>\n')
    total = 0
    for v in bestResult[6:9]:
        v = int(100*v)
        total += v
        of.write('              <td>' + str(v) + '%</td>\n')
    of.write('              <td>' + str(100-total) + '%</td>\n')
    of.write('            </tr>\n')
    of.write('            <tr class="info">\n')
    of.write('              <td>3</td>\n')
    total = 0
    for v in bestResult[9:12]:
        v = int(100*v)
        total += v
        of.write('              <td>' + str(v) + '%</td>\n')
    of.write('              <td>' + str(100-total) + '%</td>\n')
    of.write('            </tr>\n')
    of.write('          </tbody>\n')
    of.write('        </table></br>\n')
    for dist in bestDistributionFits['matches']:
        of.write('      <div class="row">\n')
        of.write('        <div class="span12">\n')
        of.write('          <h2>' + dist['title'] + '</h2>\n')
        of.write('        </div>\n')
        of.write('        <div class="span9">\n')
        of.write('          <h3> Super Codon = ' + str(dist['superNtNumber']) + '</h3>\n')
        of.write('        </div>\n')
        of.write('        <div class="span3">\n')
        of.write('          <h3> r = %0.2f, R = %0.2f</h3>\n' % (dist['stats']['corr'], dist['stats']['relEntropy']))
        of.write('        </div>\n')
        of.write('        <img src="' + dist['title'] + '.png"/>\n')
        of.write('      </div>\n')
    of.write('      </div>\n')
    of.write(superCodonTools.HTML_SUFFIX + '\n')

################################################################################
def get_aa_limits(aaLimitsStr):
    aaLimits = {}
    if(aaLimitsStr != ''):
        w = aaLimitsStr.split(',')
        
        for limit in w:
            w2 = limit.split(':')
            try:
                aaLimits[w2[0]] = float(w2[1])
            except:
                pass
    
    return aaLimits

################################################################################
def run(desiredDistributionsFN, outputDir, objectiveFunction='4', numThreads=8, aaLimitsStr='Stop:0.1', numSuperNucleotides=4):
    desiredDistributions = load_desired_distributions(desiredDistributionsFN)
    aaLimits = get_aa_limits(aaLimitsStr)
    
    pool = multiprocessing.Pool()
    results = []
    initial_vats = get_initial_vats(numThreads, numSuperNucleotides)
    #optimize_supernts(desiredDistributions, aaLimits, objectiveFunction, initial_vats[0])
    for initial_vat in initial_vats:
        results += [pool.apply_async(optimize_supernts, args=(desiredDistributions, aaLimits, objectiveFunction, initial_vat))]
    
    pool.close()
    pool.join()
    
    bestResult = get_best_result([r.get(timeout=1) for r in results])
    indiv_obj_fun = get_obj_fun(objectiveFunction)
    bestDistributionFits = find_best_dist_fits(bestResult, desiredDistributions, 
                                               aaLimits, 
                                               indiv_obj_fun, False)

    make_plots(bestDistributionFits, outputDir)
    print_results(bestResult, bestDistributionFits, outputDir)
    

################################################################################
if __name__ == '__main__':
    args = miscUtils.check_options(sys.argv, OPTIONS)
    args['num_threads'] = int(args['num_threads'])
    if(args['num_threads'] > MAX_NUM_THREADS):
        args['num_threads'] = MAX_NUM_THREADS

    run(args['desired_distributions'], args['output_dir'], args['objective_function'], args['num_threads'], args['aaLimits'], args['num_super_nts'])
    
