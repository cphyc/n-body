#!/bin/sh

while true; do
    /usr/bin/inotifywait *.f90 &&
	make &&
	notify-send 'Recompiled' ||
	    notify-send 'Error in compilation' && sleep 10
done
