HTML_PREFIX = '''<!DOCTYPE html>
<html>
  <head>
    <title>Super Codons</title>
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <!-- Bootstrap -->
    <link href="../../css/bootstrap.css" rel="stylesheet" media="screen">
    <link href="../../css/font-awesome.css" rel="stylesheet">
    <link href="../../css/default.css" rel="stylesheet" media="screen">
  </head>
  <body>
    <script src="http://code.jquery.com/jquery-latest.js"></script>
    <script src="../../js/bootstrap.min.js"></script>
    <div class="container top-level">
      <div class="row">
        <div class="span9">
          <h1>Super Codons</h1>
        </div>
        <div class="span3">
          <img src="../../img/supercodon_logo_small.png"/>
        </div>
      </div>

'''

HTML_SUFFIX = '''
    </div>
  </body>
</html>'''

DEFAULT_AA_ORDER = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 
                    'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 
                    'T', 'V', 'W', 'Y', 'Stop']
                    
AA_ORDER_DICT = dict([(i,a) for a, i in enumerate(DEFAULT_AA_ORDER)])
                    
THREE_LETTER_AA_DICTIONARY = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'}

BLOSUM_DESIRED_DISTRIBUTIONS_PICKLE = 'blosum_AA_dists.pickle'
# This info is now stored in the blosum_AA_dists.pickle
#BLOSUM62_DESIRED_DISTRIBUTIONS = {'A': 'Ala,0.5,0.0150832488,0.0205570903,0.0283382296,0.0155169445,0.0551593606,0.010505385,0.0302366108,0.0317407292,0.0419076153,0.0127204545,0.0184879943,0.0204566789,0.0182503255,0.0222599007,0.0593006359,0.0352631674,0.048065611,0.0038204377,0.0123295799',
#                                  'C': 'Cys,0.0623362799,0.5,0.0156431799,0.0150276258,0.0201433337,0.030162853,0.0090059735,0.0429464578,0.0196474948,0.0614724937,0.0147839543,0.0171860721,0.0141547165,0.0121263668,0.0154301819,0.0408926241,0.0364693943,0.0533337675,0.0056736837,0.0135635468',
#                                  'D': 'Asp,0.0334990359,0.006168076,0.5,0.0759750342,0.0117496873,0.0389799159,0.0147419875,0.019083248,0.0377929082,0.0234951569,0.007175764,0.0574966102,0.0191233173,0.0254722014,0.0245205233,0.043329054,0.0293055238,0.0203366245,0.0025065675,0.0092487641',
#                                  'E': 'Glu,0.0391312574,0.0050210607,0.0643800528,0.5,0.0111555659,0.025382352,0.0179001935,0.0159673831,0.054118333,0.0262327824,0.008893831,0.0289396533,0.0186183291,0.046343149,0.0352740982,0.0386861221,0.0268388586,0.0222481112,0.0034691045,0.0113997622',
#                                  'F': 'Phe,0.0280717405,0.0088175614,0.0130442399,0.0146151625,0.5,0.0205497273,0.013905069,0.05226159,0.0162787715,0.092913311,0.0204211643,0.0128692642,0.0090127819,0.0093091596,0.0159882473,0.0204860999,0.0199448531,0.044204458,0.0145701963,0.0727366024',
#                                  'G': 'Gly,0.080014283,0.0105870442,0.0346991254,0.0266642313,0.0164774964,0.5,0.0131819857,0.0190535544,0.0349349492,0.0287055868,0.0100820472,0.0393365822,0.0187674851,0.0188250746,0.0236927784,0.0527519071,0.0300727671,0.0250645958,0.0056051004,0.0114834057',
#                                  'H': 'His,0.0327076448,0.0067845645,0.0281658378,0.040359329,0.0239302544,0.0282924033,0.5,0.017153038,0.035061207,0.0291437636,0.0112989955,0.042229143,0.0141058501,0.0309726923,0.0366372408,0.0326351181,0.0219608949,0.0191590174,0.0044788785,0.0449241269',
#                                  'I': 'Ile,0.0322184118,0.0110726706,0.0124782203,0.0123212189,0.0307815514,0.0139958322,0.0058704956,0.5,0.0158217285,0.1150305971,0.025347017,0.0100490594,0.0101748309,0.0090046819,0.0125676109,0.0174075596,0.0272391855,0.1209950663,0.0036566994,0.0139675629',
#                                  'K': 'Lys,0.0398010555,0.0059612713,0.0290815358,0.0491440575,0.0112833013,0.0301987878,0.0141210585,0.018619182,0.5,0.0292970818,0.0108077719,0.0290209188,0.0187585042,0.036829847,0.0741385158,0.0368867402,0.0279088217,0.0230249486,0.0032346759,0.0118819245',
#                                  'L': 'Leu,0.0357688352,0.0126954089,0.0123060743,0.0162145581,0.0438355695,0.016890013,0.0079895107,0.0921412723,0.0199415321,0.5,0.0399032071,0.0110828414,0.0114512864,0.0130932226,0.0195850431,0.019652358,0.0269062459,0.0767296192,0.0059260768,0.0178873253',
#                                  'M': 'Met,0.0320213022,0.0090049434,0.0110849526,0.0162134018,0.0284153969,0.0174959203,0.0091356383,0.059881435,0.0216967578,0.1176880357,0.5,0.0126402649,0.0097467638,0.0176658047,0.0191738793,0.0204378099,0.0241028354,0.0552068758,0.0047446023,0.0136433798',
#                                  'N': 'Asn,0.0319503533,0.0071864848,0.0609757911,0.0362183033,0.0122935257,0.0468634433,0.0234401742,0.0162982232,0.0399962557,0.0224401181,0.008677721,0.5,0.0140985339,0.0250846786,0.0324302199,0.051574291,0.0366504129,0.0196915591,0.0026486385,0.0114812725',
#                                  'P': 'Pro,0.0553846587,0.0092727712,0.0317721733,0.0365042952,0.0134880784,0.0350277418,0.0122663773,0.0258529756,0.0405018509,0.0363242608,0.0104828274,0.0220872904,0.5,0.0217410298,0.0245951019,0.0427860359,0.0346934356,0.031972406,0.0036339818,0.0116127081',
#                                  'Q': 'Gln,0.035726596,0.005743885,0.0305996768,0.0656985032,0.0100732211,0.0254044288,0.0194743371,0.0165431567,0.0574967437,0.0300300274,0.0137378285,0.0284147408,0.0157197914,0.5,0.0461533791,0.0351482287,0.0256157465,0.0216501821,0.0042215908,0.0125479363',
#                                  'R': 'Arg,0.0346608073,0.0058135379,0.0234301204,0.0397759244,0.0137611008,0.0254321542,0.0183231826,0.0183652633,0.0920622636,0.0357295982,0.0118601199,0.0292199506,0.0141452205,0.036711145,0.5,0.0334655973,0.0262911334,0.0233484483,0.0039241735,0.0136802586',
#                                  'S': 'Ser,0.0700276293,0.0116844713,0.0313992024,0.0330837083,0.0133723036,0.0429437478,0.0123782133,0.0192919863,0.0347378281,0.0271902306,0.0095875594,0.0352417506,0.0186619821,0.0212027647,0.0253800927,0.5,0.052645163,0.0264176483,0.0032179379,0.0115357803',
#                                  'T': 'Thr,0.0484417248,0.0121221832,0.0247045589,0.0266999982,0.0151448853,0.0284789108,0.0096897152,0.0351173211,0.030574692,0.043305149,0.0131531585,0.0291334075,0.0176031842,0.0179756362,0.0231948804,0.0612416214,0.5,0.0474278478,0.0037230569,0.0122680685',
#                                  'V': 'Val,0.0475527689,0.0127672576,0.0123466617,0.0159398195,0.0241737641,0.0170943957,0.0060880369,0.1123409228,0.018166126,0.0889390603,0.0216969233,0.011272898,0.0116832133,0.0109416261,0.0148348758,0.0221322427,0.0341567453,0.5,0.0033383678,0.0145342942',
#                                  'W': 'Trp,0.0307023597,0.0110325974,0.0123613917,0.0201894706,0.0647233382,0.0310522974,0.0115608728,0.0275788964,0.0207305927,0.0557974131,0.0151468419,0.0123167311,0.0107866609,0.0173306061,0.0202530946,0.021899098,0.0217801084,0.027117614,0.5,0.067640015',
#                                  'Y': 'Tyr,0.0294979174,0.0078518204,0.013578635,0.0197509404,0.0961906709,0.0189393682,0.0345211632,0.0313612014,0.0226700342,0.0501391113,0.0129666767,0.0158945132,0.0102617689,0.0153353734,0.0210194782,0.0233711444,0.02136589,0.0351476004,0.0201366923,0.5'}

RESULTS_DIR = '/super_codon_results/'
OUT_DIR_PREFIX = 'sc_'

DIST_FILE = 'desired_aa_distributions.csv'

REDIRECT_CODE = '''<!DOCTYPE HTML>
<html lang="en-US">
    <head>
        <meta charset="UTF-8">
        <meta http-equiv="refresh" content="1;url=%s">
        <script type="text/javascript">
            window.location.href = "%s"
        </script>
        <title>Page Redirection</title>
    </head>
    <body>
        If you are not redirected automatically, follow the <a href='%s'>link.</a>
    </body>
</html>'''

################################################################################
def print_errors(errs):
    print HTML_PREFIX
    print '<FONT SIZE=7 COLOR="red">Error!</FONT><br>'
    for e in errs:
        print e + '<br>\n'
    print HTML_SUFFIX

