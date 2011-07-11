#general c++ compiling rule
bindir = /var/www/localhost/cgi-bin
SUBDIRS = . nff
AM_CPPFLAGS = -O3 -march=native -fomit-frame-pointer -fexpensive-optimizations
bin_PROGRAMS = xabs.cgi xrforward.cgi
xrforward_cgi_SOURCES = xrforward.cpp
#xrfoward_cgi_LDADD = multilayer.o elements.o 
xabs_cgi_SOURCES = xabs.cpp
#xabs_cgi_LDADD = multilayer.o elements.o 
LDADD = multilayer.o elements.o 
