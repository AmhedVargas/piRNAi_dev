# piRNAi app
Shiny app that aids the design of synthetic fragments required for piRNA interference. See the app in action [here](https://wormbuilder.dev/piRNAi_dev/). Also, see the [piRNA-DB repository](https://github.com/AmhedVargas/piRNAi-DB) for more information regarding the creation of target sequences.

## General description
Piwi-interacting RNAs, or simply piRNAs, are a class of small RNAs (21 nucleotides long) that are active in the *C. elegans* germline and act in the regulation of gene expression. One of the mechanisms by which they can perfom this control is by targeting cytosolic mRNA transcripts by means Watson-Crick pairing. Knowing the rules of pairing allows to design synthetic piRNAs that can potentially silence any gene endogenous or exogenous to *C.elegans*. Indeed this reasoning is the basis for a novel method of germline-specific gene regulation in *C. elegans* denominated piRNA interfence (piRNAi).

This app constitutes the first tool that aid in the design of piRNAi fragments targeting endogenous *C. elegans* and *C. briggsae* genes. 

## Structure of the repository
The shiny app core lies onto two R scripts, the "server.R" and "ui.R" codes, and a text-based database of piRNAi target sequences (contained within the "DataBase" folder). While the ui code handles the user queries, the server code computes them and load one fragment of the database via un unique identifier (see below).

**Basic processing diagram**

ui.R: Display list of available genes -> User query: Select gene and isoform -> server.R: Load piRNA sequences by gene and isoform and compute best piRNAi sequences given user input -> ui.R: Display results

Please note that the database is composed of 66542 text files (one per isoform of a given *C. elegans* or *C. briggsae* gene) which contains piRNAi sequences, their location and their level of uniqueness. Instead of using any kind of structured database format, using simple text files requires no additional instalation without sacrificing too much computational performance.

## Dependencies
**R**

While the program has been succesfully tested in R 3.6.1 (2019-07-05, and later versions) on a x86_64-pc-linux-gnu platform, the application should work fine in any other version as long as the following libraries are succesfully installed:

*   install.packages("shiny")
*   install.packages("shinythemes")
*   install.packages("DT")
*   install.packages("Biostrings")

## Deployment and implementation
You can see the app in action by going to [here](https://wormbuilder.dev/piRNAi_dev/).

To run the app locally, make sure you have installed R along with all its dependences. Follow by cloning this repository:

`git clone https://github.com/AmhedVargas/piRNAi_dev`

Run the program via command line specifying an open port, e.g., 5100, and open a web browser to access the app.

`R -e "shiny::runApp('piRNAi_dev',host="0.0.0.0",port=5100)"`

Alternatively, you can run the app in a graphical environment such as Rstudio.

## Usage
**Designer tab**
Simple design: Write the name or the WormBaseID of a *C. elegans/C. briggsae* gene. Follow by selecting an isoform to target and the specificity required, then let the program do the rest.
Advanced design: In this form you can either write the piRNAi target sequences you want to use in a given piRNA cluster (4 options each with different piRNA disposition), or search across the piRNAi target database sequences that fit your needs.

# To modify from here


**Download tab**
Download here the coordinates of piRNAi targeting sequences for *C. elegans* and *C. briggsae*.

## Some characteristics to consider on the selection of piRNAi target sequences
By design, piRNAi sequences target gene isoforms at no preferential location, i.e. 5’ or 3’ ends of the transcript.

![piRNAi distribution inside gene bodies](https://github.com/AmhedVargas/piRNAi_v2/tree/master/img/SF1.jpg?raw=true)

Similarly, the GC content of most of the piRNAi target sequences are optimal for synthesis.  

![piRNAi GC content](https://github.com/AmhedVargas/piRNAi_v2/tree/master/img/SF2.jpg?raw=true)

Finally, assuming a simple set of rules (at least 6 piRNAi targets separated by at least 5 bp), the total number of gene isoforms that could be effectively targeted by piRNAi sequences is above 90% for most of the tracks.

![Isoforms effectively targeted by piRNAi](https://github.com/AmhedVargas/piRNAi_v2/tree/master/img/SF3_combined.jpg?raw=true)

## Construction of piRNAi text-based database from piRNAi-DB tracks
In the "src" folder you can find the scripts required to produce the "DataBase" directory. However, please note that it requires as input the results produced by the [piRNA-DB repository](https://github.com/AmhedVargas/piRNAi-DB).
 
## Troubleshoot

Please feel free to [e-mail me](mailto:amhed.velazquez@kaust.edu.sa) for any question, doubt or error in the app.

