#!/bin/sh
grep 'Total Score' $1 | cut -d' ' -f 4,6-19 | sort -n -r | head -1
