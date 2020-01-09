#!/usr/bin/make -f

CC       ?= gcc
CFLAGS   ?= -O2
CPPFLAGS ?=
LDFLAGS  ?=

JVM_VER   =java-11-openjdk-amd64
JVM_HOME  =/usr/lib/jvm/${JVM_VER}/include

CFLAGS   += -fPIC -I${JVM_HOME} -I${JVM_HOME}/linux/
CPPFLAGS +=
LDFLAGS  += -shared -fPIC

export JAVA_TOOL_OPTIONS=-Dfile.encoding=UTF8

LIBNAME   =libthesiaslib
SOVERSION =0
SONAME    =${LIBNAME}.so.${SOVERSION}

PREFIX = /usr/

all: thesias

thesias: libthesiaslib.so thesias.jar

libthesiaslib.so:
	${CC} ${CPPFLAGS} ${CFLAGS} -c src/*.c
	${CC} ${LDFLAGS} -Wl,-soname,${SONAME} *.o -o ${SONAME}

thesias.jar:
	javac -d class java/*.java
	jar cfe thesias.jar GraficT -C class . -C misc LogoThesias.png

.PHONY: clean install uninstall

install:
	install -d $(DESTDIR)$(PREFIX)/lib/
	install libthesiaslib.so.0 $(DESTDIR)$(PREFIX)/lib/
	install -d $(DESTDIR)$(PREFIX)/share/java/
	install thesias.jar $(DESTDIR)$(PREFIX)/share/java/
	install -d $(DESTDIR)$(PREFIX)/bin/
	install misc/THESIAS $(DESTDIR)$(PREFIX)/bin/
	install -d $(DESTDIR)$(PREFIX)/share/man/man1/
	install misc/THESIAS.1 $(DESTDIR)$(PREFIX)/share/man/man1/
	ln -s $(DESTDIR)$(PREFIX)/lib/libthesiaslib.so.0 $(DESTDIR)$(PREFIX)/lib/libthesiaslib.so

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/lib/libthesiaslib.so.0
	rm -f $(DESTDIR)$(PREFIX)/lib/libthesiaslib.so
	rm -f $(DESTDIR)$(PREFIX)/share/java/thesias.jar
	rm -f $(DESTDIR)$(PREFIX)/bin/THESIAS
	rm -f $(DESTDIR)$(PREFIX)/share/man/man1/THESIAS.1

clean:
	rm -f *.o
	rm -f *.so
	rm -f *.so.*
	rm -f *.jar
	rm -rf class/

