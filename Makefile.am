#general c++ compiling rule
bindir = /var/www/localhost/cgi-bin
SUBDIRS = . nff
AM_CPPFLAGS = -O3 -march=native -fomit-frame-pointer -fexpensive-optimizations
bin_PROGRAMS = xabs.cgi xrforward.cgi t0
t0_SOURCES = t0.cpp
xrfoward_cgi_SOURCES = xrfoward.cpp
#xrfoward_cgi_LDADD = multilayer.o elements.o 
xabs_cgi_SOURCES = xabs.cpp
#xabs_cgi_LDADD = multilayer.o elements.o 
LDADD = multilayer.o elements.o 
