#!/usr/bin/env python2

import xml.etree.ElementTree as ET
parser = ET.XMLParser()
parser._parser.UseForeignDTD(True)
parser.entity['rarr'] = u'<'
e = ET.parse('pythia-xml/ParticleData.xml', parser=parser).getroot()

#
# Find all particles Pythia knows about within a given mass range
#

# ps = e.findall('particle')
# for p in ps:
#     if float(p.get("m0")) > 0.7 and float(p.get("m0")) < 0.8:
#         print p.get("name"), "mass", p.get("m0"), "width", p.get("mWidth")

#
# Find particular decay channels from the Pythia xml
#

particles = e.findall('particle')
for p in particles:
    decays = p.findall("channel")
    for d in decays:
        products = [int(prod) for prod in d.get("products").split(" ") if prod != ""]
        # Find all dipion decays (change as needed)
        if len(products) == 3 \
           and abs(products[0]) == 211 and abs(products[1]) == 211  and abs(products[1]) == 211:
            #
            print p.get("id"), p.get("name"), "mass", p.get("m0"), "width", p.get("mWidth"), "bRatio", d.get("bRatio")

