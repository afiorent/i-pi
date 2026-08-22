#!/bin/bash
# remove everything a run produces, leaving only the inputs
cd "$(dirname "$0")" || exit 1
rm -f escape.* log.* run.out RESTART '#'*
