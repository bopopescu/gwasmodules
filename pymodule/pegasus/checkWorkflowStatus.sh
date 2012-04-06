#!/bin/sh

workflowDir=$1
watch -d -n 8 "pegasus-status --noutf8  $workflowDir |tail -n 20"
