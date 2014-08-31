all: $(MAKEFILE_LIST) pdf html

FILES := $(filter-out $(MAKEFILE_LIST),$(shell git ls-files --cached))

.SUFFIXES:
MAKEFLAGS += --no-builtin-rules

.DELETE_ON_ERROR:

LYX := /Applications/LyX.app/Contents/MacOS/lyx

monk.makefile: monk/monk.py Monkfile
	./monk/monk.py @Monkfile --files $(FILES) > $@ || rm $@

include monk.makefile

texclean:
	grep OUTPUT *.fls | cut -d ' ' -f 2- | uniq | xargs rm -f

clean:
	git clean -dfx

.PHONY: all
