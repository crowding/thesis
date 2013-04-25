#!/bin/bash

#wait for changes to the file set and always rebuild
while kqwait $(git ls-files --cached); do
    make $@
done
