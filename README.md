# rendseq

<h4 align="center">
A python package to facilitate your End-Enriched RNA Sequencing (Rend)
data analysis dreams!</h4>


<div align="center">
  <!-- badges: start -->
  <div id="banner" style="overflow: hidden;justify-content:space-around;">
  <a href="https://gwlilabmit.github.io/rendseq/">
  <img src="https://img.shields.io/badge/ReadTheDocs-blue"
    alt="Documentation Link"></img>
  </a>
  <a href="https://www.tidyverse.org/lifecycle/#experimental">
  <img src="https://img.shields.io/badge/lifecycle-experimental-orange.svg"
    alt="Lifecycle: experimental"></img>
  </a>
<a href="https://codecov.io/gh/miraep8/rendseq">
  <img src="https://codecov.io/gh/miraep8/rendseq/branch/main/graph/badge.svg?token=SIGSJGCZPI"/>
</a>
  </div>
  <hr>
  <!-- badges: end -->

  <p>
    <a href="#overview">Overview</a> •
    <a href="#installation">Installation</a> •
    <a href="#contribute">Contribute</a>
    a href="#our_documentation">Documentation</a>
  </p>
</div>
  <!-- badges: end -->

# Overview
To learn more about RendSeq - eg what it is, what types of questions it can help you solve and how you can generate your own RendSeq data set please visit the
 [Li lab website](http://gwli.scripts.mit.edu/group/), or check out
the [Li Lab's 2018 Cell Paper](http://gwli.scripts.mit.edu/group/wp-content/uploads/2019/01/Lalanne_Cell2018.pdf)
by Lalanne et al.

# Installation

To get started - install our package via pip

  `pip install rendseq`

To learn more about how to use this package -please check out our [documentation page](https://rendseq.readthedocs.io/en/latest/).  (built by Quarto)

We also have example notebooks which can be accessed in the ["example notebooks"](https://github.com/gwlilabmit/rendseq/tree/main/example%20noteboooks) subfolder".  This folder also includes some example data to see proper formatting and to play around with.  These example notebooks can also be viewed on our [documentation page]().

# Contribute

We welcome collaborators on this project!  To get started check out our instructions for how to get started as a contributor.

You can also direct all correspondence about this project to our dedicated mailing list: <mark >rendseq *at* mit *dot* edu</mark>

# Our Documentation

Documentation which explains how this package works and also all the included functions, with example code can be found at our [Github Pages site](https://gwlilabmit.github.io/rendseq/)

Our docs are made with [Quarto](https://quarto.org/) and hosted on Github Pages.

In order to update the docs yourself this is what you need to do:

- Install [Quarto](https://quarto.org/docs/get-started/).
- Install [quartodocs](https://machow.github.io/quartodoc/get-started/overview.html).
- Run `quartodoc build` in the docs folder.
- Run `quarto render` in the docs folder.


To add new functions so that they are rendered and included in the docs - add them to the `quartodoc` section in the `_quarto.yml` file. The docstrings will then be added and processed automatically.

Github pages then treats the `docs\` folder as the source for the Github Pages website.

Hopefully the above ^ will be integrated into a Github action soon, but in the meantime this should work fine.
