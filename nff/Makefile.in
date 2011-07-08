htdir = /var/www/localhost/htdocs
cgidir = /var/www/localhost/cgi-bin/nff

clean:
all:
	echo "Nothing to do"

install:
	install -m755 -o apache -d ${htdir}; \
	install -m755 -o apache -t ${htdir} *.html;
	install -m755 -o apache -d ${cgidir}; \
	install -m755 -o apache -t ${cgidir} *.dat *.nff;


