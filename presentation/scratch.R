load("../data.RData")
library(ggplot2)

staircase_example <- match_df(data, list(colnames=))

staircase_example <-
  match_df(data,
           data.frame(subject="pbm", exp_type="spacing",
                      folded_direction_content=0.15, target_number_all=12))

(ggplot(staircase_example, aes(trial.order, folded_displacement,
                               shape=folded_response_with_carrier,
                               color=folded_response_with_carrier,
                               fill=folded_response_with_carrier))
 + geom_point()
 + labs("PBM, $C=0.15$, spacing=$3.5^\\circ$")
 + scale_x_continuous("Trial no.")
 + scale_y_continuous("$\\Delta x")
 + with_arg(name="Response", breaks=c(FALSE, TRUE), labels=c("CW", "CCW"),
            scale_shape_manual(values=c(24, 25)), scale_color_discrete())
 )
