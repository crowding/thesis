FILES := $(filter-out $(MAKEFILE_LIST),$(shell git ls-tree --name-only HEAD .))

all: $(MAKEFILE_LIST) figures pdf

monk.makefile: monk/monk.py Monkfile
	./monk/monk.py @Monkfile --files $(FILES) > $@ || rm $@

include monk.makefile

clean:
	git clean -dfx

