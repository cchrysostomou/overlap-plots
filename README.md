# overlap-plots

Visualizing overlap and correlation between datasets

## Description

This uses a Matlab script to plot both the overlap and correlation of expression between multiple datasets. It was written to plot the overlap
of CDR3 usage in antibodies across multiple samples. In this specific case, correlation of expression refers to the number of times
a unique CDR3 sequence was observed across all datasets. But this script can be applied to any variable that has a variable of interest and 
records its occurrence across multiple samples. 

The plot can be visualized as such: Similar to a venn diagram, the program will analyze the overlap between all combinations of samples. 
So for three samples, it will analyze the features found only in sample 1, or sample 2, or sample 3; the features found between sample 1 and 2, 
sample 1 and 3, and sample 2 and 3; and finally it will analyze the features found in all samples. For each subset (7 subsets described in this example),
it will plot a circle. The radius of the circle is proportional to the total number of features found in that subset (area of circle = overlap).
Then, it will plot the correlation of each individual feature in a subset across all samples; this correlation is illustrated by plotting the frequency
of each feature in a specific sample. Essentiall, for each feature in a subset, a heatmap is plotted within the circle. Features that are plotted with 
dark colors in all three samples (subset 7 -> overlap between all three) represents features observed at high frequencies in all samples and are probably
higly correlated. 

## Example of plots 
[Cartoon of description above](example-plots/cartoon.png)

## Inputs

In order to run the program, a matlab variable must be used. The variable assumes that each row represents a unique feature and that different columns
represent different samples with the number of observations for each feature. For example:

Variable = {'name','sample1','sample2';'alpha',10,0;'beta',10,20} : alpha occurs in sample 1 10 times and 0 times in sample 2, etc. 

The matlab variable may either be a matrix or a cell


## General Usage

Before trying these functions, please ensure that the folder or all files within matlab-code have been added to the Matlab path

	%load a sample matlab variable 
	load('sample-data/strong_overlap.mat')
	%plot the overlap in this sample. columns 2,3 and 4 correspond to the number of observations for each unique row (column 1) in each sample
	PlotOverlap(strong-overlap,[2,3,4], 1, 2, 1, {'a','b','c'},'enlarge',2,'trpThr',[0.5,1])


##### Required fields

1) data - A matlab variable defining observations of each unique feature across multiple samples (see inputs above)

2) colUse - A vector defining which columns in the 'data' variable (1) contain the observations for each sample you want to compare

3) contains_header -  An integer 0 or 1 defining whether the 'data' variable (1) contains a header

4) type - An integer 0, 1, or 2 defining the method to use for plotting 

	Type = 1: plot by the ranking of each feature
	Type = 2: plot by the frqeuency of each feature 
	Type = 3: plot by the CDF of each feature (that is sum the frequencies of all features before it)
	
5) UseLogData - An integer 0 or 1 defining whether to consider the log of frequencies since most samples are highly polarized

6) Names - A 1-D cell defining the names for each sample 

##### Optional fields

The following optional fields will allow you to change how the plot will appear. They should go after all required fields and be listed as 'string' value 
{'enlarge':number} - This shows how far subsets should be plotted from one another. If circles are overlapping with one another, 
						 then this value needs to be increased. If they are plotted too far from one another and appear very small, then it needs to be decreased.
{'trpThr':[min,max]} - This will define the minimimum transparency value and maximum transparency value to use while plotting
{'colors': {[64x1],[64x1],[64x1]} } - This should be a cell of three colormaps. These colormaps will replace the default colormaps used in the program.

## Dependencies

Running of this program currently requires Matlab 2010 or higher

