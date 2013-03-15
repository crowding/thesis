#!/bin/Rscript
source(ifelse(file.exists("modeling"), "modeling/library.R", "library.R"))
library(plyr)

##initial parameters

MeilstrupBoyntonModel <-
  function(
    file="data.Rd", 
    filter=c(expType="spacing", subject="pbm")
    ) {

    #Load the data
    load(file)

    #filter the current subject and experiment type
    data <- merge(data, filter)

    #get some shorter names and drop the other columns
    columns_used <- c(folded_direction_content="content",
                      folded_response_with_carrier="response",
                      folded_displacement="dx",
                      spacing="target_spacing");
    data <- subset(rename(data, columns_used), select=columns_used)

    
    
    #fit the motion model
    #####
    

    #calculate number and probability of responses from data
    
    #draw predictions from the motion model for us to plot
}

run_as_command(MeilstrupBoyntonModel)

