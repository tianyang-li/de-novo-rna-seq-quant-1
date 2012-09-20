#!/bin/sh

cd lib
mkdir tmp
cd tmp
ar x ../libaccpmparam.a
ar x ../libaccpmcore.a
ar x ../libaccpmla.a
ar x ../libaccpm.a
ar x ../libaccpmoracle.a

ar cru ../liboboe.a *.o
ranlib ../liboboe.a

cd ..
rm -rf tmp
