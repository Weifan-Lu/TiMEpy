# Makefile for GOTIC2
#    by Koji Matsumoto
#    2004.10.25

# -----< Fortran Compiler Configuration >-----

FC=/usr/bin/f95
FFLAGS   = -O3

PACKAGE = gotic2-source
PROG    = GOTIC2_mod2020

SRCS    = agravi.f  angld.f   astro.f   atilti.f  azmth.f   \
          banner_mod2020.f  chop.f    cintgl_mod2020.f  convl1.f  convl2.f  \
          convl3_mod2020.f  convl4_mod2020.f  csearc.f  dgravi.f  dtilti.f  \
          dverti.f  earth_mod2020.f   factl.f   fconf.f   gellip.f  \
          getapj.f  getapw.f  glsph.f   gmodel.f  greenf_mod2020.f  \
          greenp_mod2020.f  lpout_mod2020.f   lpout2_mod2020.f  main_mod2020.f    minmax.f  \
          mjdymd.f  omeshf.f  predic.f  quadr.f   rd1st.f   \
          rd2nd.f   rd3rd.f   rd4th.f   rdcmp_mod2020.f   rdhead.f  \
          rdispi.f  rdprein_mod2020.f readin_mod2020.f  refphs.f  reslt2.f  \
          result.f  setap.f   spharm.f  stide_mod2020.f   strani_mod2020.f  \
          tdispi.f  tjulid.f  vset.f

OBJS    = $(SRCS:.f=.o)

FILES   = Makefile $(SRCS)

VER     = `date +%Y%m%d`

BINDIR  = ..

SHELL   = /bin/sh
RM      = rm -f

### rules ###

.SUFFIXES:
.SUFFIXES: .o .f

all: $(PROG) install

$(PROG): $(OBJS)
	$(FC) -o $@ $(OBJS)

.f.o:
	$(FC) $(FFLAGS) -c $<

install:
	cp $(PROG) $(BINDIR)

### useful commands ###

.PHONY: clean
clean:
	$(RM) $(OBJS)
	$(RM) core gmon.out *~ #*#

tar:
	@echo $(PACKAGE)-$(VER) > .package
	@$(RM) -r `cat .package`
	@mkdir `cat .package`
	@ln $(FILES) `cat .package`
	tar cvf - `cat .package` | gzip -9 > `cat .package`.tar.gz
	@$(RM) -r `cat .package` .package

zip:
	zip -9 $(PACKAGE)-$(VER).zip $(FILES)
