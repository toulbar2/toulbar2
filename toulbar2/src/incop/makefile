# Makefile for INCOP library 


#CCC=/usr/local/gcc/bin/g++ 
CCC=g++


LIBS= -lm 


NEWOPT = -O3 -DNDEBUG
# -O6

#DEBUG = -g -Wall 
DEBUG= -Wall -static -static-libgcc -O3 -DNDEBUG  
# DEBUG= -g -Wall -pg
# DEBUG= -Wall


CFLAGS=$(DEBUG) $(NEWOPT) 
LDFLAGS=  $(LIBS) -static -static-libgcc



# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o 

# ------------------------------------------------------------




timer.o:		timer.cc timer.h makefile
	$(CCC) -c $(CFLAGS)  -o timer.o timer.cc

urbcsp2.o :		urbcsp2.cc makefile
	$(CCC) -c $(CFLAGS)  -o urbcsp2.o urbcsp2.cc



incoputil.o: 		incoputil.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o incoputil.o incoputil.cc


incopalgo.o:	        incopalgo.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o incopalgo.o incopalgo.cc


move.o: 		move.cc move.h  incop.h makefile
	$(CCC) -c $(CFLAGS)  -o move.o move.cc

autotuning2.o: 		autotuning2.cc  autotuning2.h incop.h  makefile
	$(CCC) -c $(CFLAGS)  -o autotuning2.o autotuning2.cc

csproblem.o:		csproblem.cc incop.h  makefile
	$(CCC) -c $(CFLAGS)  -o csproblem.o csproblem.cc

extcsp.o:		extcsp.cc extcsp.h  incop.h makefile
	$(CCC) -c $(CFLAGS)  -o extcsp.o extcsp.cc

narycsproblem.o: 	narycsproblem.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o narycsproblem.o narycsproblem.cc

narycspmain.o: 	narycspmain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o narycspmain.o narycspmain.cc

autonarycspmain.o: 	autonarycspmain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o autonarycspmain.o autonarycspmain.cc



autonarycsp:		 incopalgo.o  timer.o csproblem.o incoputil.o  autonarycspmain.o narycsproblem.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) autonarycspmain.o narycsproblem.o incopalgo.o  timer.o csproblem.o incoputil.o  autotuning2.o -o autonarycsp $(LDFLAGS)


narycsp:		 incopalgo.o  timer.o csproblem.o incoputil.o  narycspmain.o narycsproblem.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) narycspmain.o narycsproblem.o incopalgo.o  timer.o csproblem.o incoputil.o  autotuning2.o -o narycsp $(LDFLAGS)


color.o:		color.cc color.h  incop.h makefile
	$(CCC) -c $(CFLAGS)  -o color.o color.cc


colormain.o: 		colormain.cc incop.h color.h makefile
	$(CCC) -c $(CFLAGS)  -o colormain.o colormain.cc

colorcsp:		incopalgo.o  color.o timer.o csproblem.o incoputil.o colormain.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  color.o  timer.o csproblem.o colormain.o incoputil.o autotuning2.o -o colorcsp $(LDFLAGS)

impassemain.o: 		impassemain.cc incop.h color.h makefile
	$(CCC) -c $(CFLAGS)  -o impassemain.o impassemain.cc

autocolor:	incopalgo.o  color.o  timer.o csproblem.o incoputil.o autotuning2.o autocolormain.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  color.o  timer.o csproblem.o autotuning2.o autocolormain.o incoputil.o -o autocolor $(LDFLAGS)

autocolormain.o: 		autocolormain.cc incop.h color.h makefile
	$(CCC) -c $(CFLAGS)  -o autocolormain.o autocolormain.cc




impassecolor:		incopalgo.o  color.o  timer.o csproblem.o incoputil.o impassemain.o 
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  color.o timer.o csproblem.o impassemain.o incoputil.o  -o impassecolor $(LDFLAGS)

randomcsp.o: 		randomcsp.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o randomcsp.o randomcsp.cc


randommain.o: 		randommain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o randommain.o randommain.cc

randomcsp:		incopalgo.o  extcsp.o urbcsp2.o timer.o csproblem.o incoputil.o randommain.o randomcsp.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  extcsp.o urbcsp2.o timer.o csproblem.o randommain.o randomcsp.o incoputil.o  autotuning2.o -o randomcsp $(LDFLAGS)

autorandommain.o: 		autorandommain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o autorandommain.o autorandommain.cc

autorandom:		incopalgo.o  extcsp.o urbcsp2.o timer.o csproblem.o incoputil.o autorandommain.o randomcsp.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  extcsp.o urbcsp2.o timer.o csproblem.o autorandommain.o randomcsp.o incoputil.o autotuning2.o -o autorandom $(LDFLAGS)



celar.o:		celar.cc celar.h  incop.h makefile
	$(CCC) -c $(CFLAGS)  -o celar.o celar.cc


celarmain.o: 		celarmain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o celarmain.o celarmain.cc


celarincr.o: 		celarincr.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o celarincr.o celarincr.cc

celarincr:		incopalgo.o  celar.o celarincr.o  timer.o csproblem.o incoputil.o celarmain.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  celar.o celarincr.o  timer.o csproblem.o celarmain.o incoputil.o autotuning2.o -o celarincr $(LDFLAGS)

celarfullincr.o: 		celarfullincr.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o celarfullincr.o celarfullincr.cc

celarfullincr:		incopalgo.o  celar.o celarfullincr.o timer.o csproblem.o incoputil.o celarmain.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  celar.o celarfullincr.o  timer.o csproblem.o celarmain.o incoputil.o autotuning2.o -o celarfullincr $(LDFLAGS)


autocelarmain.o: 		autocelarmain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o autocelarmain.o autocelarmain.cc


autocelar:		incopalgo.o  celar.o celarincr.o timer.o csproblem.o incoputil.o autocelarmain.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  celar.o celarincr.o timer.o csproblem.o autotuning2.o autocelarmain.o incoputil.o -o autocelar $(LDFLAGS)



latinsquare.o: 		latinsquare.cc latinsquare.h  move.h incop.h makefile
	$(CCC) -c $(CFLAGS)  -o latinsquare.o latinsquare.cc

sblatinsquare.o: 		sblatinsquare.cc latinsquare.h  move.h incop.h makefile
	$(CCC) -c $(CFLAGS)  -o sblatinsquare.o sblatinsquare.cc


latinsquaremain.o: 		latinsquaremain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o latinsquaremain.o latinsquaremain.cc


latinsquare:		incopalgo.o  latinsquare.o timer.o csproblem.o incoputil.o latinsquaremain.o move.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  latinsquare.o  timer.o csproblem.o latinsquaremain.o incoputil.o move.o autotuning2.o  -o latinsquare $(LDFLAGS)

sblatinsquare:		incopalgo.o  latinsquare.o sblatinsquare.o timer.o csproblem.o incoputil.o sblatinsquaremain.o move.o  autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  latinsquare.o sblatinsquare.o timer.o csproblem.o sblatinsquaremain.o incoputil.o autotuning2.o move.o  -o sblatinsquare $(LDFLAGS)


autolatinsquare:		incopalgo.o  latinsquare.o timer.o csproblem.o incoputil.o autolatinsquaremain.o move.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  latinsquare.o  timer.o csproblem.o autolatinsquaremain.o incoputil.o move.o autotuning2.o -o autolatinsquare $(LDFLAGS)

autosblatinsquaremain.o: 		autosblatinsquaremain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o autosblatinsquaremain.o autosblatinsquaremain.cc

autosblatinsquare:		incopalgo.o  latinsquare.o sblatinsquare.o timer.o csproblem.o incoputil.o autosblatinsquaremain.o move.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  latinsquare.o sblatinsquare.o  timer.o csproblem.o autosblatinsquaremain.o incoputil.o move.o autotuning2.o -o autosblatinsquare $(LDFLAGS)

autolatinsquaremain.o: 		autolatinsquaremain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o autolatinsquaremain.o autolatinsquaremain.cc






queen.o: 		queen.cc queen.h  move.h incop.h makefile
	$(CCC) -c $(CFLAGS)  -o queen.o queen.cc


queenmain.o: 		queenmain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o queenmain.o queenmain.cc

queen:		incopalgo.o  queen.o timer.o csproblem.o incoputil.o queenmain.o move.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  queen.o timer.o csproblem.o queenmain.o incoputil.o move.o -o queen $(LDFLAGS)


maxclique.o: maxclique.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o maxclique.o maxclique.cc

maxcliquemain.o:	maxcliquemain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o maxcliquemain.o maxcliquemain.cc

autocliquemain.o:	autocliquemain.cc incop.h makefile
	$(CCC) -c $(CFLAGS)  -o autocliquemain.o autocliquemain.cc

maxclique:		incopalgo.o  move.o timer.o csproblem.o incoputil.o maxclique.o maxcliquemain.o color.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  move.o timer.o csproblem.o incoputil.o maxclique.o maxcliquemain.o color.o autotuning2.o 	 -o maxclique $(LDFLAGS)

autoclique:		incopalgo.o  move.o timer.o csproblem.o incoputil.o maxclique.o autocliquemain.o color.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  move.o timer.o csproblem.o incoputil.o maxclique.o autocliquemain.o color.o autotuning2.o -o autoclique $(LDFLAGS)


wextcsp.o:		wextcsp.cc makefile
	$(CCC) -c $(CFLAGS)  -o wextcsp.o wextcsp.cc


wextcspmain.o:	wextcspmain.cc makefile
	$(CCC) -c $(CFLAGS)  -o wextcspmain.o wextcspmain.cc

wextcsp:		incopalgo.o  timer.o csproblem.o incoputil.o wextcsp.o wextcspmain.o  extcsp.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  timer.o csproblem.o incoputil.o wextcsp.o wextcspmain.o extcsp.o -o wextcsp $(LDFLAGS)


maxcspmain.o:		maxcspmain.cc makefile
	$(CCC) -c $(CFLAGS)  -o maxcspmain.o maxcspmain.cc

maxcsp:		incopalgo.o  timer.o csproblem.o incoputil.o maxcspmain.o  extcsp.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  timer.o csproblem.o incoputil.o  maxcspmain.o extcsp.o -o maxcsp $(LDFLAGS)


carseq.o: 		carseq.cc makefile
	$(CCC) -c $(CFLAGS) -o carseq.o carseq.cc

ncarseq.o: 		ncarseq.cc makefile
	$(CCC) -c $(CFLAGS) -o ncarseq.o ncarseq.cc

ficarseq.o: 		ficarseq.cc makefile
	$(CCC) -c $(CFLAGS) -o ficarseq.o ficarseq.cc

autocarseqmain.o: 		autocarseqmain.cc makefile
	$(CCC) -c $(CFLAGS) -o autocarseqmain.o autocarseqmain.cc

carseqmain.o: 		carseqmain.cc makefile
	$(CCC) -c $(CFLAGS) -o carseqmain.o carseqmain.cc

carseq:     incopalgo.o  timer.o csproblem.o incoputil.o carseqmain.o  carseq.o  move.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  timer.o csproblem.o incoputil.o  carseqmain.o carseq.o move.o autotuning2.o -o carseq $(LDFLAGS)

autocarseq:     incopalgo.o  timer.o csproblem.o incoputil.o autocarseqmain.o  carseq.o  move.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  timer.o csproblem.o incoputil.o  autocarseqmain.o carseq.o move.o autotuning2.o -o autocarseq $(LDFLAGS)


ncarseq:     incopalgo.o  timer.o csproblem.o incoputil.o carseqmain.o  ncarseq.o  move.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  timer.o csproblem.o incoputil.o  carseqmain.o ncarseq.o move.o autotuning2.o -o ncarseq $(LDFLAGS)

autoncarseq:     incopalgo.o  timer.o csproblem.o incoputil.o autocarseqmain.o  ncarseq.o  move.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  timer.o csproblem.o incoputil.o  autocarseqmain.o ncarseq.o move.o autotuning2.o -o autoncarseq $(LDFLAGS)

ficarseq:     incopalgo.o  timer.o csproblem.o incoputil.o carseqmain.o  ficarseq.o  move.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  timer.o csproblem.o incoputil.o  carseqmain.o ficarseq.o move.o autotuning2.o -o ficarseq $(LDFLAGS)

autoficarseq:     incopalgo.o  timer.o csproblem.o incoputil.o autocarseqmain.o  ficarseq.o  move.o autotuning2.o
	$(CCC) $(CFLAGS) $(LDFLAGS) incopalgo.o  timer.o csproblem.o incoputil.o  autocarseqmain.o ficarseq.o move.o autotuning2.o -o autoficarseq $(LDFLAGS)
