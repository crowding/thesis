FILES := $(filter-out $(MAKEFILE_LIST),$(shell git ls-tree --name-only HEAD .))

LYX := /Applications/LyX.app/Contents/MacOS/lyx

monk.makefile: monk/monk.py Monkfile
	./monk/monk.py @Monkfile --files $(FILES) > $@ || rm $@

include monk.makefile

clean:
	git clean -dfx

all: $(MAKEFILE_LIST) figures pdf
