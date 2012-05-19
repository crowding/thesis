(TeX-add-style-hook "figure_1"
 (lambda ()
    (LaTeX-add-bibliographies
     "bibliography")
    (LaTeX-add-labels
     "sec:constant"
     "fig:task"
     "fig:sigmoids"
     "fig:stimuli"
     "fig:scaling"
     "fig:constant")
    (TeX-add-symbols
     "biblio")
    (TeX-run-style-hooks
     "latex2e"
     "subfiles10"
     "subfiles"
     "../manuscript.tex")))

