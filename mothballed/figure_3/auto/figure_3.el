(TeX-add-style-hook "figure_3"
 (lambda ()
    (LaTeX-add-bibliographies
     "bibliography")
    (LaTeX-add-labels
     "sec:grid"
     "sec:grid_methods"
     "sec:stim-vari-scal"
     "sec:grid-motion-energy"
     "sec:grid-results"
     "sec:spatial-frequency"
     "sec:grid-temporal-frequency"
     "sec:step-size"
     "fig:grid-tableaux"
     "fig:grid-scaling"
     "fig:grid-physical")
    (TeX-add-symbols
     "biblio")
    (TeX-run-style-hooks
     "latex2e"
     "subfiles10"
     "subfiles"
     "../manuscript.tex")))

