#!/bin/sh

while true; do
    /usr/bin/inotifywait *.f90 &&
	make &&
	notify-send 'Recompiled'
done
