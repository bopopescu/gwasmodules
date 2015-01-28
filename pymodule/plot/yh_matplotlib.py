#!/usr/bin/env python
"""
2008-10-07 simple functions related to matplotlib
"""
import os,sys
from scipy import stats

def assignMatPlotlibHueColorToLs(name_ls, debug=0):
	"""
	2008-10-07
		assign continuous HSL color spectrum to a list (in order).
		
	"""
	if debug:
		sys.stderr.write("Assigning matplotlib hue color to a list ...")
	import ImageColor
	no_of_names = len(name_ls)
	value_step = 1./(no_of_names-1)	#-1 to cover both minimum and maximum in the spectrum
	name2fc = {}
	for i in range(no_of_names):
		name = name_ls[i]
		value = i*value_step
		hue_value = max(min(int(round((1-value)*255)), 255), 0)	#max(min()) makes sure it's 0-255
		fc = ImageColor.getrgb('hsl(%s'%hue_value+',100%,50%)')#"hsl(hue, saturation%, lightness%)" where hue is the colour given as an
		# angle between 0 and 360 (red=0, green=120, blue=240),
		#saturation is a value between 0% and 100% (gray=0%, full color=100%), and lightness is a value between 0% and 100% (black=0%, normal=50%, white=100%).
		fc = [color_value/255. for color_value in fc]	#matplotlib accepts rgb in [0-1] range
		name2fc[name] = fc
	if debug:
		sys.stderr.write("Done.\n")
	return name2fc

def drawName2FCLegend(ax, name_ls, name2fc=None, shape_type=1, no_face_color=False, no_edge_color=False, title=None, font_size=4, alpha=1):
	"""
	2008-10-08
		add option title and linewidth
	2008-10-07
		draw a legend according to name_ls and colors according to name2fc.
		the shape of the legend is symmetric. It's either circle or square.
	"""
	sys.stderr.write("Drawing name2fc legend  ...")
	import matplotlib
	from matplotlib.patches import Polygon, Circle, Ellipse, Wedge
	import numpy
	
	if name2fc is None:
		name2fc = assignMatPlotlibHueColorToLs(name_ls)
	no_of_names = len(name_ls)
	value_step = min(1./(no_of_names), 0.5)	#can't exceed half of plot
	
	radius = value_step/4.
	center_x_pos = 0.25
	center_y_pos = value_step/2.
	xs = [center_x_pos-radius, center_x_pos+radius, center_x_pos+radius, center_x_pos-radius]	#X-axis value for the 4 points of the rectangle starting from lower left corner. 
	ys = numpy.array([center_y_pos-radius, center_y_pos-radius, center_y_pos+radius, center_y_pos+radius])
	for i in range(no_of_names):
		name = name_ls[i]
		fc = name2fc[name]
		facecolor = fc
		edgecolor = fc
		
		if no_face_color:
			facecolor='w'
		if no_edge_color:
			edgecolor='w'
			linewidth = 0
		else:
			linewidth = matplotlib.rcParams['patch.linewidth']
		if shape_type==1:
			patch = Circle((center_x_pos, center_y_pos), radius=radius, linewidth=linewidth, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
			center_y_pos += value_step
		elif shape_type==2:
			patch = Polygon(zip(xs, ys), linewidth=linewidth, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
			ys += value_step	#increase y-axis
		else:
			patch = Circle((center_x_pos, center_y_pos), radius=radius, linewidth=linewidth, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
			center_y_pos += value_step
		ax.add_patch(patch)
		ax.text(0.75, center_y_pos-value_step, name, horizontalalignment ='left', verticalalignment='center', size=font_size)
	if title:
		ax.set_title(title, fontsize=font_size)
	sys.stderr.write("Done.\n")

def autoLabelBarChartWithHeight(axe, rects):
	"""
	2010-4-1
		from http://matplotlib.sourceforge.net/examples/pylab_examples/barchart_demo.html
		
		attach some text labels for each rectangular in rects (return of pylab.bar())
	"""
	for rect in rects:
		height = rect.get_height()
		axe.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%.2f'%float(height), ha='center', va='bottom', fontsize='xx-small')

def setFontAndLabelSize(base_size=9):
	"""
	2012.12.3 make the legend & axes.title font size 1.4 X base_size
	2010-4-14
		function to change matplotlib font size globally
		#have to be set in the very beginning of a program, otherwise, no effect.
		
		for more additional customization check: http://matplotlib.sourceforge.net/users/customizing.html
	"""
	from matplotlib import rcParams
	rcParams['font.size'] = base_size	#default is 12
	rcParams['legend.fontsize'] = int(base_size*1.6)	#default is large 
	#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
	rcParams['axes.labelsize'] = base_size	#default is medium
	rcParams['axes.titlesize'] = int(base_size*1.6)	#default is large 
	rcParams['xtick.labelsize'] = base_size	#default is medium
	rcParams['ytick.labelsize'] = base_size	#default is medium

def setDefaultFigureSize(figsize=None):
	"""
	2012.9.6
		function to change matplotlib default figure size globally
		#have to be set in the very beginning of a program, otherwise, no effect.
		for more additional customization check: http://matplotlib.sourceforge.net/users/customizing.html
		
		below is a different way to changing figure size. not sure how inch translates into pixels.
			fig = matplotlib.pyplot.gcf()
			fig.set_size_inches(185, 60)
		
		For this figsize (X,Y) if dpi=N, then final dimension is (X*N, Y*N)
	"""
	from matplotlib import rcParams
	if figsize is None:
		figsize = (10,10)
	rcParams['figure.figsize'] = figsize

def setPlotDimension(left=0.1, right=0.9, bottom=0.1, top=0.9):
	"""
	2012.12.4
	
	for more additional customization check: http://matplotlib.sourceforge.net/users/customizing.html
		
	# The figure subplot parameters.  All dimensions are a fraction of the
	# figure width or height
	#figure.subplot.left    : 0.125  # the left side of the subplots of the figure
	#figure.subplot.right   : 0.9    # the right side of the subplots of the figure
	#figure.subplot.bottom  : 0.1    # the bottom of the subplots of the figure
	#figure.subplot.top     : 0.9    # the top of the subplots of the figure
	#figure.subplot.wspace  : 0.2    # the amount of width reserved for blank space between subplots
	#figure.subplot.hspace  : 0.2    # the amount of height reserved for white space between subplots
	
	could also be in set as 
		axe_pvalue = pylab.axes([-0.1, -0.1, 1.2, 1.2], frameon=False)	#left gap, bottom gap, width, height.
		pylab.figure(axe_pvalue.figure.number)
	"""
	from matplotlib import rcParams
	rcParams['figure.subplot.left'] = left
	rcParams['figure.subplot.right'] = right
	rcParams['figure.subplot.bottom'] = bottom
	rcParams['figure.subplot.top'] = top
	


def restoreMatplotlibRCDefaults():
	"""
	2010-4-14
		counter part of setFontAndLabelSize()
	"""
	import matplotlib
	matplotlib.rcdefaults()

def drawHist(data_ls=None, title=None, xlabel_1D=None, xticks=None, outputFname=None, min_no_of_data_points=50, needLog=False, \
			dpi=200, max_no_of_bins=20, **kwargs):
	"""
	2012.11.27 rename max_no_of_data_points to min_no_of_data_points
	2012.8.6 argument min_no_of_bins renamed to max_no_of_bins, which is actually what it does.
	2011-11-28
		add argument min_no_of_bins
	2011-8-24
		add argument kwargs, xticks
	2011-4-18
		a wrapper for histogram drawing using matplotlib
	"""
	sys.stderr.write("Drawing histogram of %s data points to %s..."%(len(data_ls), outputFname))
	import pylab
	pylab.clf()
	fig = pylab.figure(figsize=(10, 10))
	#ax = pylab.axes()
	#ax = fig.gca()
	no_of_data_points = len(data_ls)
	if no_of_data_points>=min_no_of_data_points:
		no_of_bins = max(10, min(max_no_of_bins, no_of_data_points/10))
		n, bins, patches = pylab.hist(data_ls, no_of_bins, log=needLog)
		if title:
			pylab.title(title)
		else:
			title = constructTitleFromDataSummaryStat(data_ls)
			#title = 'n=%s, mean %.4f, median %.4f'%(len(data_ls), meanValue, medianValue)
			pylab.title(title)
		if xlabel_1D is not None:
			pylab.xlabel(xlabel_1D)
		if xticks:
			x_ls = bins[:-1]	#the bins has one extra element.
			pylab.xticks(x_ls, xticks)
		pylab.savefig(outputFname, dpi=dpi)
	
	sys.stderr.write("Done.\n")

def drawBarChart(x_ls, y_ls, title=None, xlabel_1D=None, xticks=None, outputFname=None, bottom=0, needLog=False, dpi=200, **kwargs):
	"""
	2011-8-15
		a wrapper for barChart drawing using matplotlib
	"""
	sys.stderr.write("Drawing barChart of %s data points to %s..."%(len(y_ls), outputFname))
	import pylab
	pylab.clf()
	pylab.bar(x_ls, y_ls, width=0.8, bottom=bottom, log=needLog, **kwargs)
	pylab.title(title)
	if xlabel_1D is not None:
		pylab.xlabel(xlabel_1D)
	if xticks:
		pylab.xticks(x_ls, xticks)
	pylab.savefig(outputFname, dpi=dpi)
	
	sys.stderr.write("Done.\n")

def drawBoxPlot(data_2D_ls, title=None, xlabel_1D=None, xticks=None, outputFname=None, dpi=200, **kwargs):
	"""
	2011-8-15
		a wrapper for boxplot drawing using matplotlib
		
		data_2D_ls could be a list of lists or a 2D array.
			The former is more efficient because boxplot converts
			a 2-D array into a list of vectors internally anyway.
	"""
	sys.stderr.write("Drawing boxplots of %s data points to %s..."%(len(data_2D_ls), outputFname))
	import pylab
	pylab.clf()
	pylab.boxplot(data_2D_ls, notch=0, sym='+', vert=1, whis=1.5, positions=None, \
				widths=None, **kwargs)
	#patch_artist = False (default) produces boxes with the Line2D artist
	#patch_artist = True produces boxes with the Patch artist
	pylab.title(title)
	if xlabel_1D is not None:
		pylab.xlabel(xlabel_1D)
	if xticks:
		pylab.xticks(range(len(xticks)), xticks)
	pylab.savefig(outputFname, dpi=dpi)
	
	sys.stderr.write("Done.\n")

def logSum(ls):
	"""
	2011-4-27
		given a list, take the sum and then log10(sum).
		an alternative to the reduce_C_function argument in drawHexbin()
	"""
	import math
	return math.log10(sum(ls))

def drawHexbin(x_ls, y_ls, C_ls, fig_fname=None, gridsize=100, title=None, xlabel=None, ylabel=None,\
			colorBarLabel=None, reduce_C_function=None, dpi=300, mincnt=None, marginals=False, xscale='linear',\
			scale='linear'):
	"""
	2011.12.22
		xscale: [ linear | log]
			Use a linear or log10 scale on the horizontal axis.
		scale: [ linear | log]
			Use a linear or log10 scale on the vertical axis.
		mincnt: None | a positive integer
			If not None, only display cells with more than mincnt number of points in the cell
		marginals: True|False
			if marginals is True, plot the marginal density as colormapped rectagles along the bottom of the x-axis and left of the y-axis
	2011-4-27
		draw 2D histogram (reduce_C_function=logSum) or any 3D plot (3rd Dimension is determined by reduce_C_function).
		default of reduce_C_function: numpy.median.
		moved from variation/src/misc.py
	2010-7-1
		add argument reduce_C_function()
	2010-6-28
		add argument xlabel & ylabel
	2010-5-11
	"""
	import pylab, numpy
	import matplotlib.cm as cm
	if reduce_C_function is None:
		reduce_C_function = numpy.median
	pylab.clf()
	pylab.hexbin(x_ls, y_ls, C=C_ls, gridsize=gridsize, reduce_C_function=reduce_C_function, cmap=cm.jet,\
				mincnt=mincnt, marginals=marginals, xscale=xscale)#, scale=scale
	pylab.axis([min(x_ls), max(x_ls), min(y_ls), max(y_ls)])
	if title is None:
		title = "gridsize %s, %s probes."%(gridsize, len(x_ls))
	pylab.title(title)
	if xlabel:
		pylab.xlabel(xlabel)
	if ylabel:
		pylab.ylabel(ylabel)
	cb = pylab.colorbar()
	if colorBarLabel:
		cb.set_label(colorBarLabel)
	if fig_fname:
		pylab.savefig(fig_fname, dpi=dpi)



def drawScatter(x_ls=None, y_ls=None, fig_fname=None, title=None, xlabel=None, ylabel=None,\
			dpi=300, logX=False, logY=False, **keywords):
	"""
	2013.07.12 added argument logX, logY
	2011-10-17
	"""
	sys.stderr.write("Drawing scatter of %s points vs %s points to %s ..."%(len(x_ls), len(y_ls), fig_fname))
	import pylab, numpy
	import matplotlib.cm as cm
	pylab.clf()
	fig = pylab.figure(figsize=(10, 10)) # You were missing the =
	ax = fig.add_subplot(1, 1, 1)
	if logY:
		ax.set_yscale('log')
	if logX:
		ax.set_xscale('log')
	pylab.plot(x_ls, y_ls, '.', **keywords)
	if title is None:
		title = constructTitleFromTwoDataSummaryStat(x_ls, y_ls)
	pylab.title(title)
	if xlabel:
		pylab.xlabel(xlabel)
	if ylabel:
		pylab.ylabel(ylabel)
	if fig_fname:
		pylab.savefig(fig_fname, dpi=dpi)
	sys.stderr.write("Done.\n")

def constructTitleFromTwoDataSummaryStat(x_ls=None, y_ls=None):
	"""
	2015.01.28 added wilcoxon signed rank test
	2012.8.21
	
	"""
	if not y_ls and x_ls:
		title = constructTitleFromDataSummaryStat(x_ls)
	elif x_ls and y_ls:
		import numpy
		n = len(x_ls)
		xMedianValue =  numpy.median(x_ls)
		yMedianValue =  numpy.median(y_ls)
		corr = numpy.corrcoef(x_ls, y_ls)[0,1]
		wilcoxonSignedRankTestMinRankSumDif, wilcoxonSignedRankTestPvalue = stats.wilcoxon(x=x_ls, y=y_ls, zero_method='wilcox', correction=False)
		title = "n=%s p=%.4g r2=%.3f xMed=%.3g yMed=%.3g"%(n, wilcoxonSignedRankTestPvalue, 
															corr, xMedianValue, yMedianValue)
	else:
		title = ""
	return title

def constructTitleFromDataSummaryStat(data_ls=None):
	"""
	2012.8.17
	"""
	import numpy
	medianValue =  numpy.median(data_ls)
	meanValue = numpy.mean(data_ls)
	std = numpy.std(data_ls)
	title = 'n=%s mean=%.3g med=%.3g std %.3g'%(len(data_ls), meanValue, medianValue, std)
	return title

if __name__ == '__main__':
	#import pdb
	#pdb.set_trace()
	import pylab
	ax = pylab.gca()
	fig = pylab.gcf()
	name_ls = ['1', '2', '3']
	name2fc = assignMatPlotlibHueColorToLs(name_ls)
	drawName2FCLegend(ax, name_ls, name2fc, shape_type=1, no_face_color=True, no_edge_color=False, font_size=10)
	pylab.show()