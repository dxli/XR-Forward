htdir = /var/www/localhost/htdocs
cgidir = /var/www/localhost/cgi-bin/nff

HTTP_GROUP="$(shell [[ `id http 2>/dev/null` ]] && echo http || echo apache)"

clean:
all:
	echo "Nothing to do"

install:
	install -m755 -o "${HTTP_GROUP}" -d ${htdir}; \
	install -m755 -o "${HTTP_GROUP}" -t ${htdir} *.html;
	install -m755 -o "${HTTP_GROUP}" -d ${cgidir}; \
	install -m755 -o "${HTTP_GROUP}" -t ${cgidir} *.dat *.nff;


