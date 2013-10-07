splits <- c("subject", "content",
            "displacement",
            "target_number_all",
            "target_number_shown",
            "spacing", "eccentricity", "bias")

relevant <- splits %v% c("n_cw", "n_obs")

model_split <- "subject"

filter_data <- mkchain(
    subset(exp_type %in% c("spacing", "content"))
    , match_df(., subset(count(., "subject"), freq>2000), on="subject")
  )

pars <- NA

format_data <- mkchain[., energy](
    do.rename(folding=FALSE)
    , mutate(bias=1)
#    , attach_motion_energy(energy) #leave this off for now, let ME models implement
#    models implement ir , add_energies
    , mutate(data, displacement=wrap(displacement,spacing))
    , mkrates(splits)
    , recast_data
    )
