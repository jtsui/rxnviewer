import pprint
from pymongo import *
import progressbar
from indigo.indigo import *
from indigo.indigo_renderer import *

pr = pprint.PrettyPrinter(indent=2)
indigo = Indigo()
renderer = IndigoRenderer(indigo)


def savePic(query_reaction, file_name):
    output_path = "static/ero_%s.png" % file_name
    rxn = indigo.loadQueryReaction(query_reaction)
    indigo.setOption("render-output-format", "png")
    indigo.setOption("render-image-size", 500, 150)
    renderer.renderToFile(rxn, output_path)


def pbar(size):
    bar = progressbar.ProgressBar(maxval=size,
                                  widgets=[progressbar.Bar('=', '[', ']'),
                                           ' ', progressbar.Percentage(),
                                           ' ', progressbar.ETA(),
                                           ' ', progressbar.Counter(),
                                           '/%s' % size])
    return bar


eros = Connection('pathway.berkeley.edu', 27017).actv01['eros']
bar = pbar(eros.find().count())
ero_iter = eros.find()
bar.start()
i = 0
for ero in ero_iter:
    ero_id = ero['_id']
    savePic(ero['readable'][1:-1].strip(), ero_id)
    i += 1
    bar.update(i)
bar.finish()
