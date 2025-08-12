#!/bin/bash
clear && for f in $(find . -iname "*.go"); do echo "file: $f"; cat $f; done
