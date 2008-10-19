DATE=`/bin/date +%F`
VER=$(DATE)
BASE=matscat

SRC =   README MANIFEST *.m

BASE_SRC = \
  checked_resonances.m compute_resonances.m compare_eigs.m \
  plot_potential1.m plot_potential.m plot_fields.m \
  form_operators.m problem_size.m \
  eval_potential.m 

BASES_SRC = \
  $(BASE_SRC) form_operators_sys.m plot_potential1s.m \
  extract_eltV1.m extract_eltV2.m extract_eltR.m

SQUARE_SRC = squarepot.m square_well.m $(BASE_SRC)
SPLINE_SRC = splinepot.m spline_well.m $(BASE_SRC)
SQUARES_SRC = squarepotsys.m square_well.m $(BASES_SRC)
SPLINES_SRC = splinepotsys.m spline_well.m $(BASES_SRC)

mfiles:
	cat $(SQUARE_SRC) > test/squarepot.m
	cat $(SPLINE_SRC) > test/splinepot.m
	cat $(SQUARES_SRC) > test/squarepotsys.m
	cat $(SPLINES_SRC) > test/splinepotsys.m

resonant1d: mfiles
	cp test/*.m ~/public_html/resonant1d/mfiles
	chmod a+r ~/public_html/resonant1d/mfiles/*

dist: tgz
	(mv matscat1d-$(VER).tar.gz distweb)
	(cd distweb; ./build)

tgz:  matscat1d-$(VER).tar.gz

matscat1d-$(VER).tar.gz: 
	ls $(SRC) | sed s:^:matscat1d-$(VER)/: >MANIFEST
	(cd ..; ln -s $(BASE) matscat1d-$(VER))
	(cd ..; tar -czvf matscat1d-$(VER).tar.gz `cat $(BASE)/MANIFEST`)
	(cd ..; rm matscat1d-$(VER))
	(mv ../matscat1d-$(VER).tar.gz .)
