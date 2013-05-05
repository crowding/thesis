evalinfo = $(info $(1))$(eval $(1))	

define MAKELINK
$(eval,
$(2): $(1)
	rm -f $$@
	ln -h $$^ $$@

all: $(2)

)
endef

all:

#what we're doing here is pulling data from elsewhere into a folder small enoug to share for writing purposes.
$(call MAKELINK,~/bibdesk/bibliography.bib,bibliography.bib)
$(call MAKELINK,../pools/Constant.Rdata,Constant.Rdata)
$(call MAKELINK,../pools/Segment.Rdata,Segment.Rdata)
$(call MAKELINK,../pools/Titrate.Rdata,Titrate.Rdata)
$(call MAKELINK,../Segment.eyemovements.Rdata,Segment.eyemovements.Rdata)
$(call MAKELINK,../pools/SF.Rdata,SF.Rdata)
$(call MAKELINK,../pools/TF.Rdata,TF.Rdata)
$(call MAKELINK,../figure_1/A.pdf,task_timeline.pdf)
$(call MAKELINK,../common.manipulations.R,common.manipulations.R)
$(call MAKELINK,../programming.R,programming.R)
$(call MAKELINK,../demo_stimuli.mov,demo_stimuli.mov)
$(call MAKELINK,../demo_counter.mov,demo_counter.mov)
$(call MAKELINK,../demo_single.mov,demo_single.mov)
$(call MAKELINK,../demo_segment.mov,demo_segment.mov)
$(call MAKELINK,../demo_eyemovements.mov,demo_eyemovements.mov) 

