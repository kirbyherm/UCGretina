# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

# Need this for writing output files larger than 2GB on 32-bit Linux
# (untested on other systems)
CPPFLAGS := -D_FILE_OFFSET_BITS=64 -std=c++17 -Wno-shadow

# Use -D to define LHTARGET, AD, SCANNING, NEUTRONS macros
# for the C preprocesssor
ifdef NEUTRONS
  CPPFLAGS += -DNEUTRONS
endif
ifdef LHTARGET
ifdef AD
  name := UCGretina_LH_AD
  CPPFLAGS += -DLHTARGET
  CPPFLAGS += -DAD
else
  name := UCGretina_LH
  CPPFLAGS += -DLHTARGET
endif
else
ifdef SCANNING
  name := UCGretina_Scan
  CPPFLAGS  += -DSCANNING
else
  name := sunShare_Calib
endif
endif

ifdef HIGHMULT
  CPPFLAGS  += -DHIGHMULT
endif

ifdef POL
  CPPFLAGS  += -DPOL
endif

############################################
## uncomment for geant v. 4.7.1 and above
###########################################
CPPFLAGS += -DG4V47

############################################
## uncomment for geant v. 4.8.0 and above
###########################################
CPPFLAGS += -DG4V48

############################################
## uncomment for geant v. 4.9.0
###########################################
CPPFLAGS += -DG4V49

############################################
## uncomment for geant v. 4.9.6
###########################################
CPPFLAGS += -DG4V496

############################################
## uncomment for geant v. 4.10
###########################################
CPPFLAGS += -DG4V10
##ROOT stuff:
CPPFLAGS += `root-config --cflags --libs`
CPPFLAGS += -g -ggdb -O0 -Wno-unused-variable           

LDFLAGS := `root-config --glibs`

G4TARGET := $(name)
G4EXLIB := true

ifndef G4LISTS_BASE                                                                                                                                                                           
  EXTRALIBS += -L$(G4LIB)/plists/$(G4SYSTEM)
  G4LISTS_BASE = $(G4INSTALL)/physics_lists

endif

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
