#!/bin/bash
# 201-8-26 wrap around mkdirhier/mkdir so that it won't report error when the directory already exists
if test -x /usr/bin/mkdirhier
then
    /usr/bin/mkdirhier $1 >& /dev/null
else
    mkdir $1 >& /dev/null
fi
echo "return code is $?"
