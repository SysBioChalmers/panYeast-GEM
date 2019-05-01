# panYeast-GEM: The pan-genome metabolic model of _Saccharomyces cerevisiae_

[![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2FpanYeast-gem.svg)](https://badge.fury.io/gh/sysbiochalmers%2FpanYeast-gem) [![Join the chat at https://gitter.im/SysBioChalmers/panYeast-GEM](https://badges.gitter.im/SysBioChalmers/panYeast-GEM.svg)](https://gitter.im/SysBioChalmers/panYeast-GEM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

* Brief Model Description:

This repository contains the pan-genome metabolic model of _Saccharomyces cerevisiae_. **This README.md is mostly copied from `yeast-GEM`.**

* Model KeyWords:
**redefine**
**GEM Category:** Species; **Utilisation:** predictive simulation, multi-omics integrative analysis, _in silico_ strain design, model template; **Field:** metabolic-network reconstruction; **Type of Model:** curated, reconstruction; **Model Source:** [Yeast 7.6](https://sourceforge.net/projects/yeast/); **Taxonomy:** _Saccharomyces cerevisiae_; **Metabolic System:** General Metabolism; **Condition:** aerobic, glucose-limited, defined media, maximization of growth.

* Last update: 2019-02-19

* Main Model Descriptors:

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
|:-------:|:--------------:|:---------:|:----------:|:-----:|
|_Saccharomyces cerevisiae_|[Yeast 7.6](https://sourceforge.net/projects/yeast/)|4014|2766|1360|

This repository is administered by xx ([@xx](https://github.com/xx)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

## Installation

### Required Software - User:
** Copied from yeast-GEM **
* Matlab user:
  * A functional Matlab installation (MATLAB 7.3 or higher).
  * The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
* Python user:
  * Python 2.7, 3.4, 3.5 or 3.6
  * [cobrapy](https://github.com/opencobra/cobrapy)

### Required Software - Contributor:

* Both of the previous Matlab requirements.
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN).
* A [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Dependencies - Recommended Software:
* For Matlab, the [libSBML MATLAB API](https://sourceforge.net/projects/sbml/files/libsbml/MATLAB%20Interface/) (version 5.17.0 is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.

### Installation Instructions
* For users: Clone it from [`master`](https://github.com/SysBioChalmers/yeast-GEM) in the Github repo, or just download [the latest release](https://github.com/SysBioChalmers/yeast-GEM/releases).
* For contributors: Fork it to your Github account, and create a new branch from [`devel`](https://github.com/SysBioChalmers/yeast-GEM/tree/devel).

## Usage
**Explain how strain specific models can be generated**

## Model Files

**Include explanation that a yeast-GEM derived pan-submodel and a non-S288c pan-submodel are kept separately (with the first one being automatically generated from yeast-GEM), while the final pan model combines these two submodels **
The model is available in `.xml`, `.txt`, `.yml`, `.mat` and `.xlsx` (the last 2 extensions only in `master`). Additionally, the following 2 files are available:
* `dependencies.txt`: Tracks versions of toolboxes & SBML used for saving the model.

## Complementary Scripts
**Should include following Matlab scripts:**
* Convert yeast-GEM to pan-submodel
* Combine pan-submodels into one pan-model
* Generate strain specific models from pan-model and table matching genes with panIDs
* Potential gap-filling to result in a functional model

## Complementary Data
** Should include the following:**
* Pan-genome data, matching panIDs with strain specific gene IDs

## Contributors

* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Sweden
* [Feiran Li](https://www.chalmers.se/en/staff/Pages/feiranl.aspx) ([@feiranl](https://github.com/feiranl)), Chalmers University of Technology, Sweden
* [Hongzhong Lu](https://www.chalmers.se/en/Staff/Pages/luho.aspx) ([@hongzhonglu](https://github.com/hongzhonglu)), Chalmers University of Technology, Sweden
* [Benjamín J. Sánchez](https://www.chalmers.se/en/staff/Pages/bensan.aspx) ([@BenjaSanchez](https://github.com/benjasanchez)), Chalmers University of Technology, Sweden
