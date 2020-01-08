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

all: thesias

thesias: libthesiaslib.so thesias.jar

libthesiaslib.so:
	${CC} ${CPPFLAGS} ${CFLAGS} -c src/*.c
	${CC} ${LDFLAGS} *.o -o libthesiaslib.so

thesias.jar:
	javac -d class java/*.java
	cd class && jar cfe thesias.jar GraficT *.class
	mv class/thesias.jar .

.PHONY: clean install uninstall

install:
	install libthesiaslib.so /usr/lib/
	install thesias.jar /usr/share/java/
	install misc/THESIAS /usr/bin/
	install misc/THESIAS.1 /usr/share/man/man1/

uninstall:
	rm -f /usr/lib/libthesiaslib.so
	rm -f /usr/share/java/thesias.jar
	rm -f /usr/bin/THESIAS
	rm -f /usr/share/man/man1/THESIAS.1

clean:
	rm -f *.o
	rm -f *.so
	rm -f *.jar
	rm -rf class/

