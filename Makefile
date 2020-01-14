#!/usr/bin/make -f

CC       ?= gcc
CFLAGS   ?= -O2
CPPFLAGS ?=
LDFLAGS  ?=

JAVA_VERSION ?=java-11-openjdk-amd64
JAVA_HOME    ?=/usr/lib/jvm/${JAVA_VERSION}
OS           ?=linux

LIBNAME   =libthesiaslib
SOVERSION =0
SONAME    =${LIBNAME}.so.${SOVERSION}

CFLAGS   += -fPIC -I${JAVA_HOME}/include -I${JAVA_HOME}/include/${OS}/
CPPFLAGS +=
LDFLAGS  += -shared -fPIC

ifeq ($(UNAME), Darwin)
	LDFLAGS += -Wl,-install_name,${SONAME}
else
	LDFLAGS += -Wl,-soname,${SONAME}
endif

export JAVA_TOOL_OPTIONS=-Dfile.encoding=UTF8

PREFIX = /usr

all: thesias

thesias: libthesiaslib.so thesias.jar

libthesiaslib.so:
	${CC} ${CPPFLAGS} ${CFLAGS} -c src/*.c
	${CC} ${LDFLAGS} *.o -o ${SONAME}

thesias.jar:
	javac -d class java/*.java
	jar cfe thesias.jar GraficT -C class . -C misc LogoThesias.png

.PHONY: clean install uninstall

install:
	install -d $(DESTDIR)$(PREFIX)/lib/jni/
	install ${SONAME} $(DESTDIR)$(PREFIX)/lib/jni/
	install -d $(DESTDIR)$(PREFIX)/share/java/
	install thesias.jar $(DESTDIR)$(PREFIX)/share/java/
	install -d $(DESTDIR)$(PREFIX)/bin/
	install misc/THESIAS $(DESTDIR)$(PREFIX)/bin/
	install -d $(DESTDIR)$(PREFIX)/share/man/man1/
	install misc/THESIAS.1 $(DESTDIR)$(PREFIX)/share/man/man1/
	ln -s $(DESTDIR)$(PREFIX)/lib/jni/${SONAME} $(DESTDIR)$(PREFIX)/lib/jni/${LIBNAME}.so

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/lib/jni/${SONAME}
	rm -f $(DESTDIR)$(PREFIX)/lib/jni/${LIBNAME}.so
	rm -f $(DESTDIR)$(PREFIX)/share/java/thesias.jar
	rm -f $(DESTDIR)$(PREFIX)/bin/THESIAS
	rm -f $(DESTDIR)$(PREFIX)/share/man/man1/THESIAS.1

clean:
	rm -f *.o
	rm -f *.so
	rm -f *.so.*
	rm -f *.jar
	rm -rf class/

