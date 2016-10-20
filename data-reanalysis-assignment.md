Data Reanalysis Assignment and Rubric
================

DUE 2pm 2016-11-10
------------------

The objective of this assignment is to use your skills in ***R*** to replicate as closely as you can a set of statistical analyses and results reported in a paper of your choosing. You should have already selected and confirmed with me a paper for which you are going to replicate some of the analyses and you should have the data in hand.

You do not need to replicate ALL of the analysis presented in the paper (the more the better), but at minimum you need to repeat at least 3 analysis, including at least one descriptive statistical analysis or visualization and one inferential statistical analysis.

For this assignment, you should prepare several files to share with me via ***GitHub*** in a new repo called "data-reanalysis-assignment": \[1\] a PDF of the paper you are reanalyzing data from, \[2\] a ".csv" file with the original data for that paper as you either downloaded them as supplementary material or received them from the paper's author, and \[3\] an ".Rmd" file where you thoroughly describe and run the code for all of the steps in your reanalysis. I should be able to take the ".Rmd" file and knit it and produce a nicely formatted ".html" report describing what you did and seeing your results.

You should also embed in your ".Rmd" file, near your own results, any images of figure from the original paper that you replicate so that I can see them together. These should be included in a folder called "img" within your repo. You can include code like the following to reference files in your "img" folder for inclusion in your document.

`<img src="img/imagename.filetype" width="###px"/>`

Where you replace *imagename.filetype* with the name of your file, e.g., "figure-1.jpeg" and *\#\#\#* with a integer number of pixels, e.g., width="200px"

If you include the headers above in this Rmd document (which lay out for ***R*** a particular set of instructions for knitting), then all of the figures associated with your chunks of code will be put into the "img" folder.

Elements of Your Report:
------------------------

You should start your reanalysis report with a short description of the study and of the specific data and reanalyses you will be doing, to orient your reader. Outline (briefly) the goal of the original paper, the data set used, and the analyses conducted, then describe which ones you will replicate. You should also demonstrate how you read in your datafile and show a few lines of raw data in your output (e.g., using `head()`).

I will be looking for you to clearly take your reader through all of the elements of data manipulation, analysis, and, if appropriate, visualization. You should provide as much coding detail, explanation, and output tables as necessary to compare your results to those published!
