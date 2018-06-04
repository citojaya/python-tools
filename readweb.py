#!/usr/bin/python
import urllib2
import re

htmlfile = urllib2.urlopen("https://www.ebay.com.au/sch/i.html?_from=R40&_trksid=p2380057.m570.l1313.TR11.TRC1.A0.H0.Xiphone+6+battery.TRS0&_nkw=iphone+6+battery&_sacat=0")
htmltext = htmlfile.read()

print htmltext
print "WORKS"
