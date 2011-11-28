#!/bin/sh

workflowDir=$1
watch -d -n 5 "pegasus-status -l $workflowDir |tail -n 30"
