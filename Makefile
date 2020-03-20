#!/usr/bin/make -f

CC       ?= gcc
CFLAGS   ?= -O2
CPPFLAGS ?=
LDFLAGS  ?=

CFLAGS   += -fPIC
CPPFLAGS +=
LDFLAGS  += -shared -fPIC

LIBNAME   =libthesiaslib
SOVERSION =0
SONAME    =${LIBNAME}.so.${SOVERSION}

UNAME_S ?=$(shell uname -s | tr '[:upper:]' '[:lower:]')

ifeq ($(UNAME_S), darwin)
	JAVA_HOME  ?=`/usr/libexec/java_home`
	LDFLAGS    += -Wl,-install_name,${SONAME}
else
	JAVA_HOME  ?=/usr/lib/jvm/default-java
	LDFLAGS    += -Wl,-soname,${SONAME}
endif

CFLAGS += -I${JAVA_HOME}/include -I${JAVA_HOME}/include/${UNAME_S}/

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
	install misc/THESIAS.desktop $(DESTDIR)$(PREFIX)/share/applications/
	install misc/THESIAS_icon.svg $(DESTDIR)$(PREFIX)/share/pixmaps/
	ln -s $(DESTDIR)$(PREFIX)/lib/jni/${SONAME} $(DESTDIR)$(PREFIX)/lib/jni/${LIBNAME}.so

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/lib/jni/${SONAME}
	rm -f $(DESTDIR)$(PREFIX)/lib/jni/${LIBNAME}.so
	rm -f $(DESTDIR)$(PREFIX)/share/java/thesias.jar
	rm -f $(DESTDIR)$(PREFIX)/bin/THESIAS
	rm -f $(DESTDIR)$(PREFIX)/share/man/man1/THESIAS.1
	rm -f $(DESTDIR)$(PREFIX)/share/applications/THESIAS.desktop
	rm -f $(DESTDIR)$(PREFIX)/share/pixmaps/THESIAS_icon.svg

clean:
	rm -f *.o
	rm -f *.so
	rm -f *.so.*
	rm -f *.jar
	rm -rf class/

