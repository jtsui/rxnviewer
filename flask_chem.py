from flask import *
app = Flask(__name__)
import csv
from indigo.indigo import *
from indigo.indigo_renderer import *
from indigo.indigo_inchi import *
import json
import pprint
from operator import itemgetter
from collections import defaultdict
import sys

pr = pprint.PrettyPrinter(indent=2)
indigo = Indigo()
renderer = IndigoRenderer(indigo)
columns = ['acc', 'chem', 'smiles', 'products', 'proc', 'rxn_id', 'desc']
input_file = 'data/ero_gbbct.csv'
acc_line = defaultdict(list)
MYURL = ''

# load entire csv into memory
with open(input_file, 'rb') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t',
                            quoting=csv.QUOTE_MINIMAL,
                            fieldnames=columns)
    for line in reader:
        if ' ' in line['acc'] or ',' in line['acc']:
            continue
        acc_line[line['acc']].append(line)
accessions = sorted(acc_line.keys())


def savePic(substrate, product, file_name):
    output_path = "static/%s.png" % file_name
    try:
        substrate = [indigo.loadMolecule(x) for x in substrate]
        product = [indigo.loadMolecule(x) for x in product]
    except IndigoException:
        return False
    rxn = indigo.createReaction()
    for s in substrate:
        rxn.addReactant(s)
    for p in product:
        rxn.addProduct(p)
    indigo.setOption("render-output-format", "png")
    indigo.setOption("render-image-size", 1000, 300)
    renderer.renderToFile(rxn, output_path)
    return True


def generatePics(acc, output):
    lines = acc_line.get(acc, [])
    for line in lines:
        substrate = [line['smiles']]
        products = json.loads(line['products'])
        output += '<br>===========================================<br>'
        output += 'Chemical: %s<br>Smiles: %s<br>Process: %s' % (line['chem'], substrate[0], line['proc'])
        output += '<br>===========================================<br>'
        i = 0
        if products == []:
            output += '<br>No Products Found.<br>'
            continue
        for prods, prob, ero_id, rxn_id in sorted(products, key=itemgetter(1), reverse=True):
            output += '<br>-------------------------------'
            output += '<br>>>>>>>>>>>>>>>>>>'
            output += '<br>ERO: %s' % (ero_id)
            output += '<br>PROB: %s' % (prob)
            output += '<br>RXN_ID: %s' % (rxn_id)
            output += '<br><img src=%s>' % url_for('static', filename='ero_%s.png' % ero_id)
            output += '<br><<<<<<<<<<<<<<<<<'
            output += '<br>-------------------------------'
            forward = prods['forward']
            reverse = prods['reverse']

            output += '<br><br>=> Forward - Applying ERO w/ chemical %s as substrate' % line['chem']
            if forward == [[]]:
                output += '<br>No products found.<br>'
            else:
                for rxn in forward:
                    success = savePic(substrate, rxn, i)
                    if success:
                        output += '<br><img src=%s>' % url_for('static', filename='%s.png' % i)
                        output += '<br>%s -> %s<br>' % (line['smiles'], ' + '.join(rxn))
                        i += 1
            output += '<br><br><= Reverse - Applying ERO w/ chemical %s as product' % line['chem']
            if reverse == [[]]:
                output += '<br>No products found.<br>'
            else:
                for rxn in reverse:
                    success = savePic(rxn, substrate, i)
                    if success:
                        output += '<br><img src=%s>' % url_for('static', filename='%s.png' % i)
                        output += '<br>%s -> %s<br>' % (' + '.join(rxn), line['smiles'])
                        i += 1
    return output


@app.route('/')
def root():
    output = 'ERO Reaction Viewer for Genbank Data. Usage /view/accession_number.<br>'
    output += '<br>Accessions:<br>'
    for acc in accessions:
        output += '<a href=%s%s>%s</a><br>' % (MYURL, acc, acc)
    return output


@app.route('/view/<accession>')
def pic(accession):
    global MYURL
    lines = acc_line.get(accession, False)
    if lines is False:
        return 'No entry with accession: %s' % accession
    next_acc = accessions[(accessions.index(accession) + 1) % len(accessions)]
    output = '<a href=%s%s>View next accession %s</a>' % (MYURL, next_acc, next_acc)
    for line in lines:
        output += '<br><br>Accession: %s<br>Description: %s<br>' % (accession, line['desc'])
        output = generatePics(accession, output)
    return output

if __name__ == '__main__':
    if len(sys.argv) == 2:
        the_file, myport = sys.argv
        MYURL = 'http://pathway.berkeley.edu:27330/view/'
        app.run(host='0.0.0.0', port=int(myport))
    else:
        MYURL = 'http://localhost:5000/view/'
        # debug=True will run with reloader enabled
        app.run(debug=True)
