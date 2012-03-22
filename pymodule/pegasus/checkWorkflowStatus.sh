#!/bin/sh

workflowDir=$1
watch -d -n 8 "pegasus-status   $workflowDir |tail -n 20"
