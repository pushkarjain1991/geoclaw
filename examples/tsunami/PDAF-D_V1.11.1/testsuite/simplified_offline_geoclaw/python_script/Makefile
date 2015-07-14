# Makefile for Clawpack code in this directory.
# # This version only sets the local files and frequently changed
# # options, and then includes the standard makefile pointed to by CLAWMAKE.

CLAW2MAKE = $(CLAW)/geoclaw/examples/tsunami/bowl-radial/Makefile

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

include $(CLAW2MAKE)

.PHONY: cleanall clean

das: $(EXE)
	python clean_connector2.py

cleanmy:
	$(MAKE) clean
	rm -f *.pyc *.txt *.data *.xyz
	rm -rf `find . -type d -name 'ens_*'`
