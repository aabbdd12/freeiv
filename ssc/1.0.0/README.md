# freeiv 1.0.0 -- the copy submitted to the SSC archive

The same files as freeiv 1.0.0 on GitHub (tag v1.0.0), flat, with the paper as an
ancillary file (freeiv_paper.pdf). It is here to test the SSC package before
the archive publishes it. As on SSC, the package file lists every file with an
`f` line: `net install` installs the programs, the help files and the dialog
only; the example data and the paper are copied by `net get` to the current
folder, never to the system directories.

```stata
net install freeiv, from("https://raw.githubusercontent.com/aabbdd12/freeiv/main/ssc/1.0.0") replace
net get freeiv, from("https://raw.githubusercontent.com/aabbdd12/freeiv/main/ssc/1.0.0") replace
```

The examples of the help read the data from the current folder, else from the
SSC archive, else from GitHub. To go back to the GitHub version:

```stata
net install freeiv, from("https://raw.githubusercontent.com/aabbdd12/freeiv/v1.0.0") replace
```

The paper: Araar, A. (2026). *freeiv: Instrument-free estimation of a linear
structural model with an endogenous regressor* (Version 1.0.0). Zenodo.
<https://doi.org/10.5281/zenodo.23190841> (all versions:
<https://doi.org/10.5281/zenodo.22770175>).
