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

.PHONY: clean

clean:
	rm -f *.o
	rm -f *.so
	rm -f *.jar
	rm -rf class/

