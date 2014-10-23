#!/bin/bash

#genetic.13x13 maternal.3x3 resid.13x13

invert -i genetic.13x13 -o gen.inv
invert -i maternal.3x3 -o mat.inv
invert -i permanent.2x2 -o perm.inv

ident 13 > ident.13
ident 2 > ident.2
ident 3 > ident.3

cmult -a gen.inv -b ident.13 -c genetic.13x13.inverted
cmult -a mat.inv -b ident.3 -c maternal.3x3.inverted
cmult -a perm.inv -b ident.2 -c permanent.2x2.inverted

rm gen.inv mat.inv perm.inv
rm ident.13 ident.2 ident.3
