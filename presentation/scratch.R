load("../data.RData")

staircase_example <- match_df(data, list(colnames=))

( ggplot(data, aes(x=trial.ix, y=folded_displacement))
  + 
