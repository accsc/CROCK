#
# CROCK v0.1 Makefile
# 2013. Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
# Department of Pharmacology UAH.
#
#
cc=gcc


all: gausimpose.c
	$(cc) gausimpose.c -o CROCK.exe -I . -w -lm -O3 -lz

clean:
	rm -f CROCK.exe 

