#!/bin/sh

workdir=$1
pegasus-run --conf $workdir/pegasus*.properties $workdir
