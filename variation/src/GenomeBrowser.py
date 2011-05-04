#!/usr/bin/env python
"""
2008-01-30
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

if __name__ == '__main__':
	import pygtk
	pygtk.require('2.0')
	import gtk, gtk.glade, gobject
	from gtk import gdk
	import gnome
	import gnome.ui
	import gnomecanvas
	
	import matplotlib
	matplotlib.use('GTKAgg')  # or 'GTK'
	from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
	
	from matplotlib.figure import Figure
	
	#from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
	#2008-02-04 use a custom navigation tool bar
	from pymodule.yh_gnome import NavigationToolbar2GTKAgg_chromosome as NavigationToolbar
	
	from pymodule import yh_gnome
	from variation.src.common import get_chr_pos_from_x_axis_pos

import numpy, traceback
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle, Polygon
from matplotlib.text import Text
from matplotlib.collections import LineCollection, Collection

from pymodule.yh_matplotlib_artists import Gene, ExonIntronCollection
from pymodule.db import TableClass
from Results2DB_250k import Results2DB_250k
from pymodule import GenomeWideResults, GenomeWideResult, DataObject, getGenomeWideResultFromFile, PassingData, CNV, SNP
from pymodule.CNV import getCNVDataFromFileInGWA
from DrawSNPRegion import DrawSNPRegion	#2008-12-16 dealWithGeneAnnotation()
import Stock_250kDB
from GeneListRankTest import GeneListRankTest


class GeneModel:
	def __init__(self, gene_id=None, chromosome=None, symbol = None, description = None, type_of_gene = None, \
				start=None, stop=None, mrna_start = None, mrna_stop = None, cds_start = None, cds_stop = None, \
				strand = None, go_id_ls=None, go_evidence=None, go_description=None, go_term_type=None):
		self.gene_id = gene_id
		self.chromosome = chromosome
		self.symbol = symbol
		self.description = description
		self.type_of_gene = type_of_gene
		self.start = start
		self.stop = stop
		self.mrna_start = mrna_start
		self.mrna_stop = mrna_stop
		self.cds_start = cds_start
		self.cds_stop = cds_stop
		self.strand = strand
		self.go_id_ls = go_id_ls
		self.go_evidence = go_evidence
		self.go_description = go_description
		self.go_term_type = go_term_type

	
get1stSplitByUnderscore = lambda x: x.split('_')[0]
get1stSplitByAnything = lambda x: x.split()[0]

class GenomeBrowser(object):
	def __init__(self):
		"""
		2008-01-30
		"""
		program_path = os.path.dirname(sys.argv[0])
		xml = gtk.Builder()
		xml.add_from_file(os.path.join(program_path, 'GenomeBrowser_gtk_builder.glade'))
		
		xml.connect_signals(self)
		xml.connect_signals({ "on_window_destroy" : gtk.main_quit })		
		self.app1 = xml.get_object("window1")

		#xml = gtk.glade.XML(os.path.join(program_path, 'GenomeBrowser_gtk_builder.glade'))
		#xml.signal_autoconnect(self)
		
		self.xml = xml
		
		#self.app1 = xml.get_object("app1")
		self.app1.connect("delete_event", gtk.main_quit)
		self.app1.set_default_size(1200, 800)
		
		self.vbox_matplotlib = xml.get_object('vbox_matplotlib')
		
		# matplotlib canvas
		fig = Figure(figsize=(8,8))
		self.canvas_matplotlib = FigureCanvas(fig)  # a gtk.DrawingArea
		self.canvas_matplotlib.set_size_request(600,400)
		self.canvas_matplotlib.mpl_connect('pick_event', self.on_canvas_pick)
		self.canvas_matplotlib.props.has_tooltip = True	# 2010-3-15 to enable tooltip
		#self.canvas_matplotlib.connect("query-tooltip", )
		self.vbox_matplotlib.pack_start(self.canvas_matplotlib)
		
		#matplotlib axes
		#self.ax = fig.add_subplot(111)
		axe_y_offset1 = 0.05	#y_offset for axe_LD, axe_strain_pca, axe_phenotype, axe_map
		axe_height1 = 0.15	#height of axe_LD or axe_snp_matrix
		axe_y_offset2 = axe_y_offset1+axe_height1
		axe_height2 = 0.75	#height of axe_gene_model
		axe_y_offset3 = axe_y_offset2+axe_height2
		
		
		axe_x_offset1 = 0.1	#
		axe_width1 = 0.85
		axe_x_offset2 = axe_x_offset1 + axe_width1
		
		self.ax = fig.add_axes([axe_x_offset1, axe_y_offset2, axe_width1, axe_height2], frameon=False)
		self.ax.grid(True, alpha=0.3)
		#self.ax.set_xticklabels([])	#remove xtick labels on ax1 because axe_LD's xtick labels cover this.
		
		self.axe_gene_model = fig.add_axes([axe_x_offset1, axe_y_offset1, axe_width1, axe_height1], frameon=False, sharex=self.ax)
		#axe_gene_model.set_xticks([])	#this will set ax1's xticks off as well because the x-axis is shared.
		self.axe_gene_model.set_yticks([])
		
		
		# matplotlib toolbar
		self.toolbar = NavigationToolbar(self.canvas_matplotlib, self.app1)
		self.vbox_matplotlib.pack_start(self.toolbar, False, False)
		
		self.textview_output = xml.get_object('textview_output')
		
		self.textbuffer_output = self.textview_output.get_buffer()
		
		#redirect stdout/stderr to textbuffer_output
		t_table=self.textbuffer_output.get_tag_table()
		tag_err=gtk.TextTag("error")
		tag_err.set_property("foreground","red")
		#tag_err.set_property("font","monospace 10")
		t_table.add(tag_err)
		tag_out=gtk.TextTag("output")
		tag_out.set_property("foreground","blue")
		#tag_out.set_property("font","monospace 10")
		t_table.add(tag_out)
		
		self.dummy_out = yh_gnome.Dummy_File(self.textbuffer_output, tag_out)
		self.dummy_err = yh_gnome.Dummy_File(self.textbuffer_output, tag_err)
		
		self.app1.show_all()
		
		self.filechooserdialog1 = xml.get_object("filechooserdialog1")
		self.entry_min_value_cutoff = xml.get_object('entry_min_value_cutoff')
		self.entry_max_value_cutoff = xml.get_object('entry_max_value_cutoff')
		self.filechooserdialog1.connect("delete_event", yh_gnome.subwindow_hide)
		
		
		self.dialog_db_connect = xml.get_object("dialog_db_connect")
		self.dialog_db_connect.connect("delete_event", yh_gnome.subwindow_hide)
		self.entry_mysql_hostname = xml.get_object("entry_mysql_hostname")
		self.entry_mysql_dbname = xml.get_object("entry_mysql_dbname")
		self.entry_postgres_hostname = xml.get_object("entry_postgres_hostname")
		self.entry_postgres_dbname = xml.get_object("entry_postgres_dbname")
		self.entry_postgres_schema = xml.get_object("entry_postgres_schema")
		self.entry_gene_annotation_picklef = xml.get_object("entry_gene_annotation_picklef")
		self.filechooserbutton_gene_annot = xml.get_object("filechooserbutton_gene_annot")
		
		self.dialog_preferences = xml.get_object("dialog_preferences")
		self.dialog_preferences.connect("delete_event", yh_gnome.subwindow_hide)
		self.checkbutton_debug = xml.get_object("checkbutton_debug")
		self.checkbutton_stdout = xml.get_object("checkbutton_stdout")
		self.checkbutton_stderr = xml.get_object("checkbutton_stderr")
		self.entry_gene_width = xml.get_object("entry_gene_width")
		self.checkbutton_draw_gene_symbol = xml.get_object("checkbutton_draw_gene_symbol")
		
		self.aboutdialog1 = xml.get_object("aboutdialog1")
		self.aboutdialog1.connect("delete_event", yh_gnome.subwindow_hide)
		
		self.dialog_cnvqc_db = xml.get_object("dialog_cnvqc_db")
		
		vbox_which_cnv_data_from_db = xml.get_object("vbox_which_cnv_data_from_db")
		self.comboboxentry_which_cnv_data_from_db = gtk.ComboBoxEntry()
		self.which_cnv_data_from_db_ls = ['CNVQCCall', 'Ref-Coverage by Fragments From One Accession (1D)',\
								'Ref-Coverage by One Fragment (2D)', 'Ref-Probe-Coverage by one Fragment (2D)',\
								'CNVCall', 'CNV', 'CNVArrayCall']
		self.fillComboBox(self.comboboxentry_which_cnv_data_from_db, self.which_cnv_data_from_db_ls)
		
		hbox_seq_ref_pos_version = xml.get_object('hbox_seq_ref_pos_version')
		hbox_seq_ref_pos_version.set_sensitive(False)
		
		vbox_which_cnv_data_from_db.pack_start(self.comboboxentry_which_cnv_data_from_db)
		self.comboboxentry_which_cnv_data_from_db.connect("changed", self.on_comboboxentry_which_cnv_data_from_db_changed)
		
		self.entry_cnv_non_ecotype_id = xml.get_object("entry_cnv_non_ecotype_id")
		self.non_ecotype_id_liststore = gtk.ListStore(gobject.TYPE_STRING)
		self.entry_cnv_non_ecotype_id.set_model(self.non_ecotype_id_liststore)
		#self.entry_cnv_non_ecotype_id = gtk.combo_box_entry_new_text()	# combo_box_entry_new_text() is not available through glade.
		self.entry_cnv_ecotype_id = xml.get_object("entry_cnv_ecotype_id")
		
		self.checkbutton_within_coverage = xml.get_object("checkbutton_within_coverage")
		
		self.radiobutton_non_ecotype_id = xml.get_object("radiobutton_non_ecotype_id")
		self.radiobutton_non_ecotype_id.connect("toggled", self.on_radiobutton_cnv_db_id_toggled, "radiobutton_non_ecotype_id")
		self.radiobutton_ecotype_id = xml.get_object("radiobutton_ecotype_id")
		self.radiobutton_ecotype_id.connect("toggled", self.on_radiobutton_cnv_db_id_toggled, "radiobutton_ecotype_id")
		self.radiobutton_ecotype_id.set_group(self.radiobutton_non_ecotype_id)
		self.radiobutton_non_ecotype_id.toggled()	#2010-4-15 default, this one is active. different from .set_active(True)
		
		self.entry_cnv_intensity_fname = self.xml.get_object("entry_cnv_intensity_fname")
		self.entry_cnv_probe_extend_dist = self.xml.get_object("entry_cnv_probe_extend_dist")
		self.filechooserbutton_CNV_intensity = self.xml.get_object("filechooserbutton_CNV_intensity")
		hbox_cnv_intensity = self.xml.get_object("hbox_cnv_intensity")
		hbox_cnv_intensity.set_sensitive(False)
		
		self.filechooserdialog_cnv_gada = xml.get_object("filechooserdialog_cnv_gada")
		
		self.dialog_cnv_by_region = xml.get_object("dialog_cnv_by_region")
		self.comboboxentry_cnv_by_region = xml.get_object("comboboxentry_cnv_by_region")
		cnv_by_region_data_type_ls = ["Probe Intensity", "Probe-Blast Sequence Fragment (2D)", \
									"Ref Coverage By Sequences (1D)", "CNVQCCall", "CNVCall",\
									"Ref Coverage By Sequences (2D)", '100bp Coverage by PESolexaData']
		self.fillComboBox(self.comboboxentry_cnv_by_region, cnv_by_region_data_type_ls)
		
		self.combobox_id_type_cnv_by_region = xml.get_object("combobox_id_type_cnv_by_region")
		cnv_by_region_id_type_ls = ["array id", "ecotype id", "accession id"]
		self.fillComboBox(self.combobox_id_type_cnv_by_region, cnv_by_region_id_type_ls)
		
		# 2010-7-27 in the dialog_cnvqc_db
		self.combobox_type_of_nearby_data = xml.get_object("combobox_type_of_nearby_data")
		type_of_nearby_data_ls = ["", "tiling probe intensity", "PE genome-wide coverage"]
		self.fillComboBox(self.combobox_type_of_nearby_data, type_of_nearby_data_ls)
		
		self.combobox_gada_id_type = xml.get_object("combobox_gada_id_type")
		id_type_ls = ["array id", "ecotype id"]
		self.fillComboBox(self.combobox_gada_id_type, id_type_ls)
		
		self.mysql_conn = self.mysql_curs = self.postgres_conn = self.postgres_curs = self.db = None
		self.gene_annotation = None
		self.candidate_gene_set = None
		
		self.chr_id2size = None
		self.chr_id2cumu_size = None
		self.chr_gap = None
		self.chr_id_ls = []
		
		self.genome_wide_results = GenomeWideResults(gap=1.0)
		self.genome_wide_results.genome_wide_result_ls = []
		self.genome_wide_results.genome_wide_result_obj_id2index = {}
		self.artist_obj_id2data_obj_key = {}
		self.yticks = []
		self.yticklabels = []
		
		self.gene_id2artist_object_id = {}
		self.chr_id2gene_id_ls = {}	#chr_id here is str type (db is varchar type)
		self.gene_id2model = {}
		self.artist_obj_id2artist_gene_id_ls = {}
		
		self.gene_id2vspan_obj_id = {}	#for the axvspan()'s drawn on the canvas
		
		self.gene_width = 1.0
		
		self.draw_gene_symbol_when_clicked = 0
				
		self.debug = 0
		
		# 2010-3-18 to cache data
		self.ecotype_id2RBDictCoverageOfRef = {}
		self.cnv_intensity_fname2data = {}
		
		self.dialog_db_connect.show_all()
		
		# 2010-4-27 redirect the stdout & stderr in the end to allow user to see any error message in the __init__()
		sys.stdout = self.dummy_out
		sys.stderr = self.dummy_err
		
	def fillComboBox(self, widget, str_ls=[]):
		"""
		2010-4-27
			fill a ComboBox or ComboBoxEntry with a list of labels from str_ls
		"""
		widget.clear()	#clear previous content. For ComboBoxEntry, it also clears out the CellRendererText()
		
		liststore = gtk.ListStore(str)
		widget.set_model(liststore)

		# this is required for ComboBox. But here also on ComboBoxEntry because of widget.clear().
		cell = gtk.CellRendererText()
		widget.pack_start(cell, True)
		widget.add_attribute(cell, 'text', 0)
		
		if hasattr(widget, "set_text_column"):	# it's a ComboBoxEntry.
			widget.set_text_column(0)	# default is -1, which would render later text-appending invisible.
		
		for label in str_ls:
			widget.append_text(label)
		
	
	def load_data(self, mysql_curs, postgres_curs):
		"""
		2008-02-04 update the info related to chromosome , position in toolbar
		2008-02-01
			read the input data
			
			chr_id2cumu_size has an extra fake chromosome (0) compared to chr_id2size
			
			chr_id is all changed into str type
		"""		
		from variation.src.common import get_chr_id2size, get_chr_id2cumu_size
		chr_id_int2size = get_chr_id2size(mysql_curs)
		self.chr_id2size = {}	#change the type of chr_id into string type
		for chr_id_int, size in chr_id_int2size.iteritems():
			self.chr_id2size[str(chr_id_int)] = size
		#chr_id2cumu_size has an extra fake chromosome (0) compared to chr_id2size
		data_ls = get_chr_id2cumu_size(self.chr_id2size)
		self.chr_id2cumu_size, self.chr_gap, self.chr_id_ls = data_ls[:3]
		#2008-02-04 update the info related to chromosome , position in toolbar
		if self.debug:
			print 'self.chr_id2size', self.chr_id2size
			for chr_id in self.chr_id2size:
				print type(chr_id)
			print 'self.chr_id2cumu_size', self.chr_id2cumu_size
			for chr_id in self.chr_id2cumu_size:
				print type(chr_id)
			print 'self.chr_id_ls', self.chr_id_ls
			for chr_id in self.chr_id_ls:
				print type(chr_id)
		self.toolbar.update_chr_info(self.chr_id2size, self.chr_id2cumu_size, self.chr_gap, self.chr_id_ls)
	
	def getXposition(self, chr, pos):
		chr = str(chr)
		if getattr(self, 'chr_id2cumu_size', None) is None:
			self.load_data(self.mysql_curs, self.postgres_curs)
		this_chr_starting_pos_on_plot = self.chr_id2cumu_size[chr]-self.chr_id2size[chr]-self.chr_gap
		x = this_chr_starting_pos_on_plot + pos
		return x
	
	def getYposition(self, genome_wide_result, y_value):
		"""
		2010-4-15
			adjust for the genome_wide_result.base_value and genome_wide_result.min_value
		"""
		return genome_wide_result.base_value - genome_wide_result.min_value + y_value
		
	def plot(self, ax, canvas, genome_wide_result, draw_line_as_point=True, drawType=1):
		"""
		2010-4-15
			argument draw_cnv_block becomes drawType
				1: CNV
				2: sequence fragment
				3: line as point
		2010-3-14
			argument draw_cnv_block: each data_obj in genome_wide_result.data_obj_ls is a cnv block.
		2008-05-28
			input is genome_wide_result
			chr_id2cumu_size, chr_id2size, chr_gap hidden from arguments
		2008-02-04
			chromosome is converted to str type
		2008-02-01
			draw the p-value, snp position, and chromosome boundary
		"""
		sys.stderr.write("Plotting %s ..."%genome_wide_result.name)
		#ax.clear()
		genome_wide_result_id = id(genome_wide_result)
		x_ls = []
		y_ls = []
		plot_later = True	# a variable deciding whether to plot the each data_obj on first encounter or later
		no_of_data_objs = len(genome_wide_result.data_obj_ls)
		no_of_objs_drawn = 0
		for i in range(no_of_data_objs):
			data_obj = genome_wide_result.data_obj_ls[i]
			y_pos = self.getYposition(genome_wide_result, data_obj.value)
			x_pos = self.getXposition(data_obj.chromosome, data_obj.position)
			if data_obj.stop_position is not None:
				x_stop_pos = self.getXposition(data_obj.chromosome, data_obj.stop_position)
				if drawType==1:
					xs = [x_pos, x_stop_pos, x_stop_pos, x_pos]
					y_base_pos = genome_wide_result.base_value - genome_wide_result.min_value + 0
					ys = [y_pos, y_pos, y_base_pos, y_base_pos]
					facecolor = getattr(data_obj, 'color', 'g')
					artist_obj = Polygon(zip(xs,ys), facecolor=facecolor, alpha=0.8, linewidth=0.5, picker=True) 	#linewidth=0.7
					ax.add_patch(artist_obj)
					artist_obj_id = id(artist_obj)
					self.artist_obj_id2data_obj_key[artist_obj_id] = [genome_wide_result_id, i]
					plot_later = False
					no_of_objs_drawn += 1
				elif drawType==3 or draw_line_as_point==False:
					x_ls.append([(x_pos, y_pos), (x_stop_pos, y_pos)])
					#artist_obj = Line2D([x_pos, y_pos], [x_stop_pos, y_pos], picker=True)
				elif drawType==2:	# 2010-4-15 draw arrows to show corresponding ref genome parts for one sequence fragment
					#sys.stderr.write("drawing fragment ...")
					target_start = getattr(data_obj, 'target_start', data_obj.value)
					target_start = self.getYposition(genome_wide_result, target_start)
					target_stop = getattr(data_obj, 'target_stop', data_obj.value)
					target_stop = self.getYposition(genome_wide_result, target_stop)
					from matplotlib.patches import ConnectionPatch
					artist_obj = ConnectionPatch(xyA=(x_pos, target_start), xyB=(x_stop_pos, target_stop), coordsA="data", \
									coordsB="data", picker=True, alpha=0.8, color='b',\
									arrowstyle="->", shrinkA=0, shrinkB=0, mutation_scale=20)
					ax.add_artist(artist_obj)
					artist_obj_id = id(artist_obj)
					self.artist_obj_id2data_obj_key[artist_obj_id] = [genome_wide_result_id, i]
					#ax.add_patch(artist_obj)
					plot_later = False
					no_of_objs_drawn += 1
			else:
				if draw_line_as_point and data_obj.stop_position is not None:
					x_stop_pos = self.getXposition(data_obj.chromosome, data_obj.stop_position)
					x_pos = (x_pos+x_stop_pos)/2.0
				x_ls.append(x_pos)
				y_ls.append(y_pos)
					#artist_obj = Circle((x_pos, y_pos), picker=True)
		
		if plot_later:
			if len(y_ls)>0:
				artist_obj = ax.scatter(x_ls, y_ls, s=10, edgecolors='none', picker=True)	#2010-4-12 faceted=False is deprecated
			else:
				artist_obj = LineCollection(x_ls, picker=True)
				ax.add_artist(artist_obj)
			artist_obj_id = id(artist_obj)
			self.artist_obj_id2data_obj_key[artist_obj_id] = [genome_wide_result_id, None]
			no_of_objs_drawn += 1
		
		y_base_value = genome_wide_result.base_value
		y_top_value = genome_wide_result.base_value + genome_wide_result.max_value - genome_wide_result.min_value
		if self.debug:
			print "y_base_value", y_base_value
			print 'y_top_value', y_top_value
		self.yticks.append(y_base_value)
		self.yticks.append(y_top_value)
		ax.set_yticks(self.yticks)
		
		self.yticklabels.append('%.2f'%(genome_wide_result.min_value))
		self.yticklabels.append('%.2f'%(genome_wide_result.max_value))
		
		title_distance = abs(genome_wide_result.max_value-genome_wide_result.min_value)/10.0
		
		# 2010-3-17 add the title on top of the plot with a blended transformation
		# the x coords of this transformation are axes, and the
		# y coord are data
		trans = matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData)	
		ax.text(0.5, y_top_value+title_distance, genome_wide_result.name, horizontalalignment='center', verticalalignment='bottom', \
				transform = trans,  bbox=dict(facecolor='white', alpha=0.4))
		
		ax.set_yticklabels(self.yticklabels)
		#draw the chromosome boundary
		for chr_id, cumu_size in self.chr_id2cumu_size.iteritems():
			ax.vlines(cumu_size, y_base_value, y_top_value, color='k')
		canvas.draw()
		sys.stderr.write("%s objects drawn. Done.\n"%(no_of_objs_drawn))
	
	def respond2GeneObjPicker(self, event, artist_obj_id, artist_obj_id2artist_gene_id_ls, gene_annotation,\
							ax=None, draw_gene_symbol_when_clicked=False, canvas_matplotlib=None):
		"""
		2010-3-15
			return output_str
		2008-12-17
			a common function specifying how to respond to picker events of 
				both gene models in axe_gene_model and vertical gene spans in ax
			called by on_canvas_pick()
		"""
		output_str = ""
		
		gene_id = artist_obj_id2artist_gene_id_ls[artist_obj_id][1]
		gene_model = gene_annotation.gene_id2model[gene_id]
		
		if len(gene_model.gene_commentaries)==0:
			gene_commentary = gene_model	#fake one here
		else:
			gene_commentary = gene_model.gene_commentaries[0]
		
		protein_label = getattr(gene_commentary, 'protein_label', None)
		if not protein_label:
			protein_label = getattr(gene_commentary, 'label', '')
		
		protein_comment = getattr(gene_commentary, 'protein_comment', None)
		if not protein_comment:
			protein_comment = getattr(gene_commentary, 'comment', '')
		
		if getattr(gene_commentary, 'protein_label', None) is not None:	#true gene_commentary is available
			type_of_gene = getattr(gene_model, 'type_of_gene', '')
		else:	#it doesn't have protein, get gene_commentary_type
			type_of_gene = getattr(gene_commentary, 'gene_commentary_type', '')
		
		output_str += '%s (gene id=%s) type_of_gene: %s. chromosome: %s. start: %s. stop: %s. strand: %s.\n'%\
				(gene_model.gene_symbol, gene_id, type_of_gene, gene_model.chromosome, \
				gene_model.start, gene_model.stop, gene_model.strand)
		output_str += '\t protein_label: %s.\n'%protein_label
		output_str += '\t protein_comment: %s.\n'%protein_comment
		
		if draw_gene_symbol_when_clicked:
			if ax:
				ax.text(event.mouseevent.xdata, event.mouseevent.ydata, gene_model.gene_symbol, size=8)
			if canvas_matplotlib:
				canvas_matplotlib.draw()
		
		return output_str
	
	def on_canvas_pick(self, event):
		"""
		2010-3-15
			CNV polygon artists are stored in artist_obj_id2data_obj_key.
			use self.canvas_matplotlib.set_tooltip_text() to show what was clicked in tooltip.
		2008-11-12
			display more information (maf, genotype_var_perc, comment) of data_obj (SNP) if they exist
		2008-05-28
			pick from collection
		2008-01-31 copied from examples/pick_event_demo.py from matplotlib source code
		"""
		output_str = ""
		if self.debug:
			output_str += "x: %s\n"%event.mouseevent.x
			output_str += "y: %s\n"%event.mouseevent.y
			output_str += "xdata: %s\n"%event.mouseevent.xdata
			output_str += "ydata: %s\n"%event.mouseevent.ydata
			output_str += "dir(event): %s\n"%dir(event)
			output_str += "dir(event.guiEvent): %s\n"%dir(event.guiEvent)
			output_str += "dir(event.mouseevent): %s\n"%dir(event.mouseevent)
			output_str += "event.artist: %s\n"%event.artist
			output_str += "dir(event.artist): %s\n"%dir(event.artist)
			output_str += "type(event.artist): %s\n"%type(event.artist)
		"""
		if isinstance(event.artist, Line2D):
			thisline = event.artist
			xdata = thisline.get_xdata()
			ydata = thisline.get_ydata()
			ind = event.ind
			if self.debug:
				print "indices:", ind
				print 'onpick1 line:', zip(numpy.take(xdata, ind), numpy.take(ydata, ind))
			for i in ind:
				print "snp chromosome: %s, position: %s, pvalue: %s"%(self.snp_pos_ls[i][0], self.snp_pos_ls[i][1], self.pvalue_ls[i])
		"""
		artist_obj_id = id(event.artist)
		if isinstance(event.artist, ExonIntronCollection):	#ExonIntronCollection is also a kind of Collection
			if artist_obj_id in self.artist_obj_id2artist_gene_id_ls:
				output_str += self.respond2GeneObjPicker(event, artist_obj_id, self.artist_obj_id2artist_gene_id_ls, self.gene_annotation,\
							ax=self.axe_gene_model, draw_gene_symbol_when_clicked=self.draw_gene_symbol_when_clicked, \
							canvas_matplotlib=self.canvas_matplotlib)
			else:
				sys.stderr.write("%s not in artist_obj_id2artist_gene_id_ls.\n"%(artist_obj_id))
		elif isinstance(event.artist, Collection) or isinstance(event.artist, LineCollection):	#
			if artist_obj_id in self.artist_obj_id2data_obj_key:
				genome_wide_result_id, data_obj_index = self.artist_obj_id2data_obj_key[artist_obj_id]
				genome_wide_result = self.genome_wide_results.get_genome_wide_result_by_obj_id(genome_wide_result_id)
				for obj_index in event.ind:
					if isinstance(obj_index, tuple) or isinstance(obj_index, list):
						obj_index = obj_index[0]
					data_obj = genome_wide_result.get_data_obj_by_obj_index(obj_index)
					output_str += str(data_obj)
			else:
				sys.stderr.write("%s not in artist_obj_id2data_obj_key.\n"%(artist_obj_id))
		elif isinstance(event.artist, Polygon):	#ExonIntronCollection is also a kind of Collection			
			patch = event.artist
			if self.debug:
				output_str += 'artist ID: %s\n'%artist_obj_id
				output_str += 'onpick1 patch: %s\n'%patch.get_verts()
			
			if artist_obj_id in self.artist_obj_id2artist_gene_id_ls:
				output_str += self.respond2GeneObjPicker(event, artist_obj_id, self.artist_obj_id2artist_gene_id_ls, self.gene_annotation,\
										ax=self.ax, draw_gene_symbol_when_clicked=self.draw_gene_symbol_when_clicked, \
										canvas_matplotlib=self.canvas_matplotlib)
			elif artist_obj_id in self.artist_obj_id2data_obj_key:	# 2010-3-15 CNV artists are stored in artist_obj_id2data_obj_key.
				genome_wide_result_id, data_obj_index = self.artist_obj_id2data_obj_key[artist_obj_id]
				genome_wide_result = self.genome_wide_results.get_genome_wide_result_by_obj_id(genome_wide_result_id)
				data_obj = genome_wide_result.get_data_obj_by_obj_index(data_obj_index)
				output_str += str(data_obj)
			else:
				sys.stderr.write("%s not in artist_obj_id2artist_gene_id_ls.\n"%(artist_obj_id))
		else:	#2010-4-15 for everything else, just check if it's in artist_obj_id2data_obj_key
			if artist_obj_id in self.artist_obj_id2data_obj_key:
				genome_wide_result_id, data_obj_index = self.artist_obj_id2data_obj_key[artist_obj_id]
				genome_wide_result = self.genome_wide_results.get_genome_wide_result_by_obj_id(genome_wide_result_id)
				data_obj = genome_wide_result.get_data_obj_by_obj_index(data_obj_index)
				output_str += str(data_obj)
			else:
				sys.stderr.write("%s not in artist_obj_id2data_obj_key.\n"%(artist_obj_id))
			
			if isinstance(event.artist, Rectangle):
				patch = event.artist
				output_str += 'onpick1 rectangle: %s\n'%patch.get_verts()
			elif isinstance(event.artist, Text):
				text = event.artist
				output_str += 'onpick1 text: %s\n'%text.get_text()
		print output_str
		self.canvas_matplotlib.set_tooltip_text(output_str)	# 2010-3

	
	def on_imagemenuitem_quit_activate(self, data=None):
		"""
		2008-02-01
			program quits
		"""
		gtk.main_quit()
	
	def on_imagemenuitem_open_activate(self, event, data=None):
		self.filechooserdialog1.show_all()
	
	def on_imagemenuitem_db_connect_activate(self, event, data=None):
		self.dialog_db_connect.show_all()
	
	def on_imagemenuitem_cnvqc_db_activate(self, event, data=None):
		"""
		2009-10-30
		"""
		self.dialog_cnvqc_db.show_all()
	
	def on_imagemenuitem_cnv_file_activate(self, event, data=None):
		"""
		2009-10-30
		"""
		self.filechooserdialog_cnv_gada.show_all()
	
	def on_imagemenuitem_cnv_by_region_activate(self, event, data=None):
		"""
		2010-4-27
			show the cnv_by_region dialog
		"""
		self.dialog_cnv_by_region.show_all()
		
		
	def on_button_filechooser_ok_clicked(self, widget, data=None):
		"""
		2008-12-16
			allow gwr name to be specified
			add function to get gwr from db based on call_method_id, analysis_method_id, phenotype_method_id
		2008-10-12
			add checkbutton_draw_line_as_point
			add checkbutton_4th_col_stop_pos
		2008-08-03
			restrict the data by (chromosome, start, stop)
		2008-05-31
			add check button to handle log10 transformation
		2008-05-28
			use GenomeWideResult and etc
		2008-02-14
			set the window title by the input filename
		"""
		input_fname = self.filechooserdialog1.get_filename()
		self.filechooserdialog1.hide()
		if not self.mysql_conn or not self.mysql_curs:
			self.db_connect()
		self.app1.set_title("Genome Browser: %s"%input_fname)
		
		checkbutton_log10_transformation = self.xml.get_object("checkbutton_log10_transformation")
		if checkbutton_log10_transformation.get_active():
			do_log10_transformation = True
		else:
			do_log10_transformation = False
		
		if self.entry_min_value_cutoff.get_text():
			min_value_cutoff = float(self.entry_min_value_cutoff.get_text())
		else:
			min_value_cutoff = None
		
		#2008-08-03
		pdata = PassingData()
		entry_chromosome = self.xml.get_object("entry_chromosome")
		if entry_chromosome.get_text():
			pdata.chromosome = int(entry_chromosome.get_text())
		entry_start = self.xml.get_object("entry_start")
		if entry_start.get_text():
			pdata.start = int(entry_start.get_text())
		entry_stop = self.xml.get_object("entry_stop")
		if entry_stop.get_text():
			pdata.stop = int(entry_stop.get_text())
		
		# 2009-10-27
		if self.entry_max_value_cutoff.get_text():
			pdata.max_value_cutoff = float(self.entry_max_value_cutoff.get_text())
		else:
			pdata.max_value_cutoff = None
		# 2009-10-27
		checkbutton_OR_min_max = self.xml.get_object("checkbutton_OR_min_max")
		if checkbutton_OR_min_max.get_active():
			pdata.OR_min_max = True
		else:
			pdata.OR_min_max = False
		
		checkbutton_4th_col_stop_pos = self.xml.get_object("checkbutton_4th_col_stop_pos")
		if checkbutton_4th_col_stop_pos.get_active():
			pdata.is_4th_col_stop_pos = True
		else:
			pdata.is_4th_col_stop_pos = False
		
		checkbutton_draw_line_as_point = self.xml.get_object("checkbutton_draw_line_as_point")
		if checkbutton_draw_line_as_point.get_active():
			draw_line_as_point= True
		else:
			draw_line_as_point = False
		
		entry_gwr_name = self.xml.get_object("entry_gwr_name")
		if entry_gwr_name.get_text():
			pdata.gwr_name = entry_gwr_name.get_text()
		else:
			pdata.gwr_name = None
		
		entry_call_method_id = self.xml.get_object("entry_call_method_id")
		call_method_id = entry_call_method_id.get_text()
		entry_analysis_method_id = self.xml.get_object("entry_analysis_method_id")
		analysis_method_id = entry_analysis_method_id.get_text()
		entry_phenotype_method_id = self.xml.get_object("entry_phenotype_method_id")
		phenotype_method_id = entry_phenotype_method_id.get_text()
		
		if call_method_id and analysis_method_id and phenotype_method_id:
			call_method_id = int(call_method_id)
			analysis_method_id = int(analysis_method_id)
			phenotype_method_id = int(phenotype_method_id)
			rows = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).filter_by(analysis_method_id=analysis_method_id).\
					filter_by(phenotype_method_id=phenotype_method_id).filter_by(results_method_type_id=1)
			if rows.count()==1:
				rm = rows.first()
			elif rows.count()==0:
				sys.stderr.write("No result fetched from db based on call_method_id=%s, analysis_method_id=%s, phenotype_method_id=%s.\n"%\
								(call_method_id, analysis_method_id, phenotype_method_id))
				rm = None
			else:
				sys.stderr.write("First result out of %s results fetched from db based on call_method_id=%s, analysis_method_id=%s, phenotype_method_id=%s.\n"%\
								(rows.count(), call_method_id, analysis_method_id, phenotype_method_id))
				rm = rows.first()
			if rm:
				input_fname = rm.filename
				pdata.gwr_name = '%s_%s_%s_call_%s'%(rm.analysis_method.short_name, \
													rm.phenotype_method_id, rm.phenotype_method.short_name,\
													rm.call_method_id)
		
		genome_wide_result = getGenomeWideResultFromFile(input_fname, min_value_cutoff, do_log10_transformation, pdata)
		if len(genome_wide_result.data_obj_ls)>0:
			self.genome_wide_results.add_genome_wide_result(genome_wide_result)
			#self.load_data(input_fname, self.mysql_curs, self.postgres_curs)
			self.plot(self.ax, self.canvas_matplotlib, self.genome_wide_results.genome_wide_result_ls[-1], draw_line_as_point=draw_line_as_point)
		else:
			sys.stderr.write("No data in %s under min_value_cutoff=%s. Maybe min_value_cutoff is too high.\n"%(input_fname, min_value_cutoff))
	
	def on_button_filechooser_cancel_clicked(self, widget, data=None):
		self.filechooserdialog1.hide()
	
	def on_button_dialog_db_connect_cancel_clicked(self, widget, data=None):
		self.dialog_db_connect.hide()
	
	def on_button_dialog_db_connect_clicked(self, widget, data=None):
		self.dialog_db_connect.hide()
		self.db_connect()
	
	def on_button_cnvqc_ok_clicked(self, widget, data=None):
		"""
		2010-3-14
			add functionality to deal with checkbutton_sequence_refpos, which allows getting sequence coverage data from db
		2009-10-30
			get CNV QC data as a genome_wide_result from db
		"""
		if not self.mysql_conn or not self.mysql_curs:
			self.db_connect()
		
		entry_cnv_non_ecotype_id = self.entry_cnv_non_ecotype_id
		if entry_cnv_non_ecotype_id.get_property('sensitive'):
			cnv_non_ecotype_id = yh_gnome.getDataOutOfTextEntry(entry_cnv_non_ecotype_id, data_type=int,\
															filter_func = get1stSplitByUnderscore)
		else:
			cnv_non_ecotype_id = None
		
		entry_cnv_ecotype_id = self.xml.get_object("entry_cnv_ecotype_id")
		if entry_cnv_ecotype_id.get_property('sensitive'):
			cnv_ecotype_id = yh_gnome.getDataOutOfTextEntry(entry_cnv_ecotype_id, data_type=int)
		else:
			cnv_ecotype_id = None
		
		if cnv_non_ecotype_id is None and cnv_ecotype_id is None:
			sys.stderr.write("Neither non-ecotype nor ecotype id is given.\n")
			return
		
		cnv_method_id = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_cnv_method_id"), data_type=int,\
													filter_func=get1stSplitByAnything, default=1)
		
		cnv_type_id = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_cnv_type_id"), data_type=int,\
													filter_func=get1stSplitByAnything, default=1)
		
		cnv_qc_min_size = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_cnv_qc_min_size"), data_type=int)
		
		cnv_qc_min_no_of_probes = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_cnv_qc_min_no_of_probes"), \
																data_type=int)
		
		min_reciprocal_overlap = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_min_reciprocal_overlap"), \
															data_type=float, default=0.6)
		
		chr = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_chr_cnv_from_db"), data_type=int)
		start = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_start_cnv_from_db"), data_type=int)
		stop = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_stop_cnv_from_db"), data_type=int)
		
		getFstSplitFunc = lambda x: x.split()[0]
		seq_ref_pos_version = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("comboboxentry_seq_ref_pos_version"), \
											data_type=int, filter_func=getFstSplitFunc, default=1)
		
		type_of_cnv_data = self.comboboxentry_which_cnv_data_from_db.get_active()
		drawType = 1	# default is to draw CNV (box)
		local_genome_wide_result_ls = []
		if type_of_cnv_data==0:	#CNVQCCall
			genome_wide_result = self.db.getCNVQCInGWA(cnv_non_ecotype_id, cnv_type_id=cnv_type_id, chr=chr, start=start, \
									stop=stop, \
									min_size=cnv_qc_min_size, cnv_method_id=cnv_method_id,\
									min_no_of_probes=cnv_qc_min_no_of_probes)
			local_genome_wide_result_ls.append(genome_wide_result)
		elif type_of_cnv_data==1:	# Ref-Coverage by Fragments From One Accession (1D)
			genome_wide_result = self.db.getSequenceRefPosInGWA(accession_id=cnv_non_ecotype_id, min_size=cnv_qc_min_size, \
									min_no_of_probes=cnv_qc_min_no_of_probes, version=seq_ref_pos_version)
			local_genome_wide_result_ls.append(genome_wide_result)
		elif type_of_cnv_data==2:	# Ref-Coverage by One Fragment (2D)
			genome_wide_result = self.db.getOneSequenceRefPosInGWA(sequence_fragment_id=cnv_non_ecotype_id, \
										min_size=cnv_qc_min_size, \
										min_no_of_probes=cnv_qc_min_no_of_probes, version=seq_ref_pos_version)
			local_genome_wide_result_ls.append(genome_wide_result)
			drawType = 2	# draw sequence fragment (arrow)
		elif type_of_cnv_data==3:	# Ref-Probe-Coverage by one Fragment (2D)
			genome_wide_result = self.db.getOneSequence2ProbeInGWA(sequence_fragment_id=cnv_non_ecotype_id, \
									min_no_of_identities=25)
			drawType = 2	# draw sequence fragment (arrow)
			local_genome_wide_result_ls.append(genome_wide_result)
		elif type_of_cnv_data == 4:	# CNVCall
			array_id_ls = []
			if cnv_non_ecotype_id:
				array_id_ls = [cnv_non_ecotype_id]
			elif cnv_ecotype_id:
				array_id_ls = self.db.getArrayIDLsGivenEcotypeID(cnv_ecotype_id)
			for array_id in array_id_ls:
				genome_wide_result = self.db.getCNVCallInGWA(array_id=array_id, cnv_type_id=cnv_type_id, \
												chr=chr, start=start, \
												stop=stop, \
												cnv_method_id=cnv_method_id, min_size=cnv_qc_min_size, \
												min_no_of_probes=cnv_qc_min_no_of_probes)
				if genome_wide_result:
					local_genome_wide_result_ls.append(genome_wide_result)
		elif type_of_cnv_data==5: #2010-8-5 CNV
			genome_wide_result = self.db.getCNVInGWA(cnv_method_id=cnv_method_id,cnv_type_id=cnv_type_id, chr=chr, start=start, \
										stop=stop, \
										min_size=cnv_qc_min_size,\
										min_no_of_probes=cnv_qc_min_no_of_probes)
			local_genome_wide_result_ls.append(genome_wide_result)
			
		elif type_of_cnv_data==6: #2010-8-5 CNVArrayCall
			array_id_ls = []
			if cnv_non_ecotype_id:
				array_id_ls = [cnv_non_ecotype_id]
			elif cnv_ecotype_id:
				array_id_ls = self.db.getArrayIDLsGivenEcotypeID(cnv_ecotype_id)
			for array_id in array_id_ls:
				genome_wide_result = self.db.getCNVArrayCallInGWA(array_id=array_id, cnv_method_id=cnv_method_id,cnv_type_id=cnv_type_id, \
										chr=chr, start=start, stop=stop, \
										min_size=cnv_qc_min_size,\
										min_no_of_probes=cnv_qc_min_no_of_probes)
				if genome_wide_result:
					local_genome_wide_result_ls.append(genome_wide_result)
		else:
			sys.stderr.write("Choose which type of CNV data first!\n")
			return
		
		for genome_wide_result in local_genome_wide_result_ls:
			checkbutton_within_coverage_sensitive_state = self.checkbutton_within_coverage.get_property('sensitive')
			if checkbutton_within_coverage_sensitive_state and self.checkbutton_within_coverage.get_active():
				ecotype_id = genome_wide_result.ecotype_id 
				if ecotype_id not in self.ecotype_id2RBDictCoverageOfRef:
					self.ecotype_id2RBDictCoverageOfRef.update(self.db.getSequenceFragmentRefPosFromDBInRBDict(data_source_id=None, \
																			ecotype_id=ecotype_id, \
									min_QC_segment_size=cnv_qc_min_size, min_no_of_probes=cnv_qc_min_no_of_probes,\
									min_reciprocal_overlap=min_reciprocal_overlap))
				RBDictCoverageOfRef = self.ecotype_id2RBDictCoverageOfRef.get(ecotype_id)
				if RBDictCoverageOfRef:
					# remove objects from genome_wide_result which fall outside of RBDictCoverageOfRef
					genome_wide_result.keepGWRObjectsWithinGivenRBDict(RBDictCoverageOfRef, min_reciprocal_overlap=min_reciprocal_overlap)
			
			# plot genome_wide_result
			if len(genome_wide_result.data_obj_ls)>0:
				self.genome_wide_results.add_genome_wide_result(genome_wide_result)
				self.plot(self.ax, self.canvas_matplotlib, self.genome_wide_results.genome_wide_result_ls[-1], drawType=drawType)
			else:
				sys.stderr.write("No CNV QC data for accession_id=%s, cnv_type=%s, min_size=%s, min_no_of_probes=%s.\n"%\
								(cnv_non_ecotype_id, cnv_type_id, cnv_qc_min_size, cnv_qc_min_no_of_probes))
				continue
			combobox_type_of_nearby_data = self.xml.get_object("combobox_type_of_nearby_data")
			type_of_nearby_data = combobox_type_of_nearby_data.get_active()
			
			input_fname = self.entry_cnv_intensity_fname.get_text()
			
			extend_dist = self.entry_cnv_probe_extend_dist.get_text()
			if extend_dist:
				extend_dist = int(extend_dist)*1000
			else:
				extend_dist = 20000
			rbDict = CNV.turnSegmentGWRIntoRBDict(genome_wide_result, extend_dist=extend_dist, \
										min_reciprocal_overlap=min_reciprocal_overlap)
			if type_of_nearby_data==1:	#intensity
				self.plotProbeIntensityNearGWR(genome_wide_result, rbDict=rbDict, cnv_intensity_fname=input_fname, \
											min_reciprocal_overlap=min_reciprocal_overlap)
			elif type_of_nearby_data==2:	#coverage
				self.plotQuanLongPECoverageWithinRBDict(rbDict, input_fname, \
							min_reciprocal_overlap=min_reciprocal_overlap, additionalTitle=None)
		self.dialog_cnvqc_db.hide()
	
	def plotProbeIntensityNearGWR(self, genome_wide_result, rbDict=None, cnv_intensity_fname=None, \
								min_reciprocal_overlap=None):
		"""
		2010-4-27
			split out of on_button_cnvqc_ok_clicked()
			This function draws intensity of probes which are near or within the genome_wide_result.
		"""
		if not os.path.isfile(cnv_intensity_fname):
			sys.stderr.write("Either CNV probe intensity filename not specified or it doesn't exist.\n")
			return
		ecotype_id = getattr(genome_wide_result, 'ecotype_id', None)
		array_id = getattr(genome_wide_result, 'array_id', None)
		checkbutton_match_intensity_by_ecotypeid = self.xml.get_object("checkbutton_match_intensity_by_ecotypeid")
		array_id_ls = []
		if checkbutton_match_intensity_by_ecotypeid.get_active():	# try ecotype ID first
			if ecotype_id:
				array_id_ls = self.db.getArrayIDLsGivenEcotypeID(ecotype_id)
			elif array_id:	#
				array_id_ls = [array_id]
		else:
			if array_id:
				array_id_ls = [array_id]
			elif ecotype_id:
				array_id_ls = self.db.getArrayIDLsGivenEcotypeID(ecotype_id)
		if array_id_ls:
			self.plotIntensityOfProbesWithinRBDict(rbDict, cnv_intensity_fname, array_id_ls, \
										min_reciprocal_overlap=min_reciprocal_overlap, \
										additionalTitle=genome_wide_result.name)
	
	def plotIntensityOfProbesWithinRBDict(self, rbDict, cnv_intensity_fname, array_id_ls, \
										min_reciprocal_overlap=0.6, additionalTitle=None):
		"""
		2010-4-28
			split out of plotProbeIntensityNearGWR()
		"""
		tilingIntensityData = self.cnv_intensity_fname2data.get(cnv_intensity_fname)
		if not tilingIntensityData:
			tilingIntensityData = CNV.TilingProbeIntensityData(cnv_intensity_fname, \
												min_reciprocal_overlap=min_reciprocal_overlap)
			"""
			data_matrix, probe_id_ls, chr_pos_ls, header = CNV.getProbeIntensityData(cnv_intensity_fname)
			col_id_ls = header[1:-2]
			col_id_ls = map(int, col_id_ls)
			tilingIntensityData = SNP.SNPData(row_id_ls=chr_pos_ls, col_id_ls=col_id_ls, data_matrix=data_matrix)
			tilingIntensityData.probe_id_ls = probe_id_ls
			"""
			self.cnv_intensity_fname2data[cnv_intensity_fname] = tilingIntensityData
		
		for array_id in array_id_ls:
			intensity_gwr = tilingIntensityData.getIntensityForOneArrayInGWRGivenRBDict(array_id, rbDict=rbDict,\
																additionalTitle=additionalTitle)
			"""
			intensity_gwr = CNV.fetchIntensityInGWAWithinRBDictGivenArrayIDFromTilingIntensity(tilingIntensityData, \
																		array_id, rbDict,\
																		gwr_name=gwr_name,\
																		min_reciprocal_overlap=min_reciprocal_overlap)
			"""
			if intensity_gwr and len(intensity_gwr.data_obj_ls)>0:
				self.genome_wide_results.add_genome_wide_result(intensity_gwr)
				self.plot(self.ax, self.canvas_matplotlib, self.genome_wide_results.genome_wide_result_ls[-1], drawType=1)
			else:
				sys.stderr.write("No intensity data for array %s.\n"%array_id)
	
	coverageFname2gwr = {}
	def plotQuanLongPECoverageWithinRBDict(self, rbDict, input_fname, \
							min_reciprocal_overlap=0.6, additionalTitle=None, windowSize=100):
		"""
		2010-7-28
			similar to plotIntensityOfProbesWithinRBDict, but replace probe intensity with PE coverage
		"""
		coverage_gwr = self.coverageFname2gwr.get(input_fname)
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey
		from pymodule.SNP import GenomeWideResult, DataObject
		from pymodule.CNV import readQuanLongPECoverageIntoGWR
		if coverage_gwr is None:
			coverage_gwr = readQuanLongPECoverageIntoGWR(input_fname, additionalTitle=additionalTitle, \
											windowSize=windowSize)
			self.coverageFname2gwr[input_fname] = coverage_gwr
		
		# restrict the coverage_gwr within rbDict
		inRBDictCoverageGWR = GenomeWideResult(name=coverage_gwr.name)
		genome_wide_result_id = id(inRBDictCoverageGWR)
		for data_obj in coverage_gwr.data_obj_ls:
			cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=data_obj.chromosome, span_ls=[data_obj.position],\
										min_reciprocal_overlap=min_reciprocal_overlap)
			if cnvSegmentKey in rbDict:
				inRBDictCoverageGWR.add_one_data_obj(data_obj)
		
		if inRBDictCoverageGWR and len(inRBDictCoverageGWR.data_obj_ls)>0:
			self.genome_wide_results.add_genome_wide_result(inRBDictCoverageGWR)
			self.plot(self.ax, self.canvas_matplotlib, self.genome_wide_results.genome_wide_result_ls[-1], drawType=1)
		else:
			sys.stderr.write("No intensity data for array %s.\n"%array_id)
	
	def on_button_cnvqc_cancel_clicked(self, widget, data=None):
		"""
		2009-10-30
		"""
		self.dialog_cnvqc_db.hide()
		
	def on_button_cnv_gada_ok_clicked(self, widget, data=None):
		"""
		2010-6-3
			deal with ecotype-id/array-id, chr, start, stop
		2009-10-30
			get CNVs from GADA output file
		"""
		input_fname = self.filechooserdialog_cnv_gada.get_filename()
		self.filechooserdialog_cnv_gada.hide()
		if not self.mysql_conn or not self.mysql_curs:
			self.db_connect()
		
		id_type = yh_gnome.getDataOutOfTextEntry(self.combobox_gada_id_type)
		gada_id = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_gada_id"), data_type=int)
		array_id_ls = []
		if id_type=='ecotype id':
			array_id_ls = self.db.getArrayIDLsGivenEcotypeID(gada_id)
			if not array_id_ls:
				sys.stderr.write("no arrays found for ecotype id %s.\n"%(gada_id))
				return
		elif gada_id:
			array_id_ls.append(gada_id)
		
		if not array_id_ls:
			sys.stderr.write("no arrays specified for gada id %s.\n"%(gada_id))
			return
		
		max_amp = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_cnv_max_amp"), data_type=float)
		min_amp = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_cnv_min_amp"), data_type=float)
		min_size = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_cnv_min_size"), data_type=int)
		
		chr = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_gada_chromosome"), \
											data_type=int)
		start = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_gada_start"), data_type=int)
		stop = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_gada_stop"), data_type=int)
		
		min_no_of_probes = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_cnv_min_no_of_probes"), data_type=int)
		
		checkbutton_filter_gada_by_cutoff = self.xml.get_object("checkbutton_filter_gada_by_cutoff")
		if checkbutton_filter_gada_by_cutoff.get_active():
			filter_gada_by_cutoff = True
		else:
			filter_gada_by_cutoff = False
		
		input_fname_ls = [input_fname]
		local_genome_wide_result_ls = []
		for array_id in array_id_ls:
			genome_wide_result = getCNVDataFromFileInGWA(input_fname_ls, array_id, max_amp=max_amp, \
									min_amp=min_amp, min_size=min_size,\
									min_no_of_probes=min_no_of_probes, chr=chr, start=start, stop=stop,\
									filter_gada_by_cutoff=filter_gada_by_cutoff)
			if len(genome_wide_result.data_obj_ls)>0:
				self.genome_wide_results.add_genome_wide_result(genome_wide_result)
				#self.load_data(input_fname, self.mysql_curs, self.postgres_curs)
				self.plot(self.ax, self.canvas_matplotlib, self.genome_wide_results.genome_wide_result_ls[-1], drawType=1)
			else:
				sys.stderr.write("No CNV data for array_id=%s, max_amp=%s, min_amp=%s, min_size=%s, min_no_of_probes=%s.\n"%\
								(array_id, max_amp, min_amp, min_size, min_no_of_probes))
	
	
	def on_button_cnv_gada_cancel_clicked(self, widget, data=None):
		"""
		2009-10-30
		"""
		self.filechooserdialog_cnv_gada.hide()

	def on_button_ok_cnv_by_region_clicked(self, widget, data=None):
		"""
		2010-4-27
		"""
		cnv_by_region_id = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("comboboxentry_id_cnv_by_region"), \
												data_type=int, \
												filter_func=get1stSplitByUnderscore)
		
		chr = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_chr_cnv_by_region"), data_type=int)
		start = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_start_cnv_by_region"), data_type=int)
		stop = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_stop_cnv_by_region"), data_type=int)	
		
		#2010-7-22
		cnv_method_id = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("comboboxentry_cnv_method_id_by_region"), \
											data_type=int, \
											filter_func=get1stSplitByAnything, default=1)
		cnv_type_id = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("comboboxentry_cnv_type_id_by_region"), \
											data_type=int, \
											filter_func=get1stSplitByAnything, default=1)
		cnv_min_size = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_min_size_by_region"), data_type=int)
		cnv_min_no_of_probes = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("entry_min_no_of_probes_by_region"), \
													data_type=int)
		
		cnv_by_region_id_type = yh_gnome.getDataOutOfTextEntry(self.combobox_id_type_cnv_by_region)
		cnv_by_region_data_type = self.comboboxentry_cnv_by_region.get_active()
		
		seq_ref_pos_version = yh_gnome.getDataOutOfTextEntry(self.xml.get_object("comboboxentry_seq_ref_pos_version_by_region"), \
								data_type=int, filter_func=get1stSplitByAnything, default=1)
		
		entry_intensity_fname = self.xml.get_object("entry_intensity_fname")
		cnv_intensity_fname = yh_gnome.getDataOutOfTextEntry(entry_intensity_fname)
		
		drawType = 1	# default is to draw CNV (box)
		local_genome_wide_result_ls = []
		if cnv_by_region_data_type==0:	#Probe Intensity
			array_id_ls = []
			additionalTitle = None
			if cnv_by_region_id_type=='ecotype id':
				additionalTitle = 'ecotype %s'%cnv_by_region_id
				array_id_ls = self.db.getArrayIDLsGivenEcotypeID(cnv_by_region_id)
				if not array_id_ls:
					sys.stderr.write("no arrays found for ecotype id %s.\n"%(cnv_by_region_id))
					return
			elif cnv_by_region_id:
				array_id_ls.append(cnv_by_region_id)
			if not array_id_ls:
				sys.stderr.write("no array id set.\n")
				return
			if not cnv_intensity_fname:
				sys.stderr.write("Please specify cnv_intensity_fname.\n")
				return
			if chr and start and stop:
				from pymodule.RBTree import RBDict
				rbDict = RBDict(cmpfn=CNV.leftWithinRightAlsoEqualCmp)
				segmentKey = CNV.CNVSegmentBinarySearchTreeKey(chromosome=chr, span_ls=[start, stop])
				rbDict[segmentKey] = (chr, start, stop)
			else:
				sys.stderr.write("at least one of chr, start, stop not specified.\n")
				return
			if rbDict and cnv_intensity_fname and array_id_ls:
				self.plotIntensityOfProbesWithinRBDict(rbDict, cnv_intensity_fname, array_id_ls, \
										additionalTitle=additionalTitle)
		elif cnv_by_region_data_type==1:	# SequenceFragment Probe-Blast Result	(2D)
			local_genome_wide_result_ls.extend(self.db.getSequence2ProbeInMultiGWA(accession_id=cnv_by_region_id, \
														min_no_of_identities=25, chr=chr, start=start, stop=stop))
			drawType = 2	# draw sequence fragment (arrow)
		elif cnv_by_region_data_type==2:	# Reference Coverage By SequenceFragment (1D)
			genome_wide_result = self.db.getSequenceRefPosInGWA(accession_id=cnv_by_region_id, chr=chr, \
															start=start, stop=stop, version=seq_ref_pos_version)
			local_genome_wide_result_ls.append(genome_wide_result)
		elif cnv_by_region_data_type == 3:	# CNVQCCall
			genome_wide_result = self.db.getCNVQCInGWA(accession_id=cnv_by_region_id, chr=chr, start=start, stop=stop,\
											cnv_type_id=cnv_type_id, cnv_method_id=cnv_method_id, min_size=cnv_min_size,\
											min_no_of_probes=cnv_min_no_of_probes)
			local_genome_wide_result_ls.append(genome_wide_result)
		elif cnv_by_region_data_type==4:	# CNVCall
			array_id_ls = []
			if cnv_by_region_id_type=='ecotype id':
				array_id_ls = self.db.getArrayIDLsGivenEcotypeID(cnv_by_region_id)
			elif cnv_by_region_id:
				array_id_ls.append(cnv_by_region_id)
			for array_id in array_id_ls:
				genome_wide_result = self.db.getCNVCallInGWA(array_id=array_id, chr=chr, start=start, stop=stop,\
											cnv_type_id=cnv_type_id, cnv_method_id=cnv_method_id, min_size=cnv_min_size,\
											min_no_of_probes=cnv_min_no_of_probes)
				local_genome_wide_result_ls.append(genome_wide_result)
		
		elif cnv_by_region_data_type == 5:	#Ref Coverage By Sequences (2D)
			local_genome_wide_result_ls.extend(self.db.getSequenceRefPosInMultiGWA(accession_id=cnv_by_region_id, \
							chr=chr, start=start, stop=stop, version=seq_ref_pos_version))
			drawType = 2	# draw sequence fragment (arrow)
		elif cnv_by_region_data_type == 6:	#100bp Coverage by PESolexaData
			from pymodule.CNV import readQuanLongPECoverageIntoGWR
			if not cnv_intensity_fname:
				sys.stderr.write("Please specify cnv_intensity_fname.\n")
				return
			coverage_fname = cnv_intensity_fname
			#the coverage data is split into different chromosomes. so no need to specify chromosome here.
			genome_wide_result = readQuanLongPECoverageIntoGWR(coverage_fname, \
											additionalTitle=None, windowSize=100, start=start, stop=stop)
			local_genome_wide_result_ls.append(genome_wide_result)
		else:
			sys.stderr.write("Type of CNV data %s unknown.\n"%cnv_by_region_data_type)
			return
		
		
		for genome_wide_result in local_genome_wide_result_ls:
			# plot genome_wide_result
			if len(genome_wide_result.data_obj_ls)>0:
				self.genome_wide_results.add_genome_wide_result(genome_wide_result)
				self.plot(self.ax, self.canvas_matplotlib, self.genome_wide_results.genome_wide_result_ls[-1], \
						drawType=drawType)
			else:
				sys.stderr.write("No data for id=%s, chr=%s, start=%s, stop=%s.\n"%\
								(cnv_by_region_id, chr, start, stop))
				continue
		self.dialog_cnv_by_region.hide()
		
	def on_button_cancel_cnv_by_region_clicked(self, widget, data=None):
		"""
		2010-4-27
		"""
		self.dialog_cnv_by_region.hide()
	
	def db_connect(self):
		"""
		2010-1-15
			pass "cls_with_db_args=self" to DrawSNPRegion.dealWithGeneAnnotation()
		2009-12-09
			add db_user, db_passwd to MySQLdb.connect()
		2008-12-16
			add gene_annotation_picklef
		2008-02-01
			read the data in dialog_db_connect and establish the connections to two databases
		"""
		sys.stderr.write("Database Connecting ...")
		self.drivername = 'mysql'
		self.hostname = self.entry_mysql_hostname.get_text()
		self.dbname = self.entry_mysql_dbname.get_text()
		self.db_user = self.xml.get_object("entry_db_user").get_text()
		self.db_passwd = self.xml.get_object("entry_db_passwd").get_text()
		
		import MySQLdb
		try:
			self.mysql_conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user=self.db_user, passwd=self.db_passwd)
			self.mysql_curs = self.mysql_conn.cursor()
			self.db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
					   password=self.db_passwd, hostname=self.hostname, database=self.dbname)
			self.db.setup(create_tables=False)
			self.session = self.db.session
		except:
			sys.stderr.write('DB connection error: %s\n'%repr(sys.exc_info()))
			traceback.print_exc()
		
		if not self.gene_annotation:
			gene_annotation_picklef = self.entry_gene_annotation_picklef.get_text()
			self.gene_annotation = DrawSNPRegion.dealWithGeneAnnotation(gene_annotation_picklef, cls_with_db_args=self)
		
		#2010-1-13 for postgresql. commented out
		#hostname = self.entry_postgres_hostname.get_text()
		#dbname = self.entry_postgres_dbname.get_text()
		#schema = self.entry_postgres_schema.get_text()
		
		#from annot.bin.codense.common import db_connect			#2008-12-16 don't need postgres conn anymore
		#self.postgres_conn, self.postgres_curs = db_connect(hostname, dbname, schema)
		
		sys.stderr.write("Done.\n")
	
	def get_gene_id2model(cls, curs, entrezgene_mapping_table='genome.entrezgene_mapping', \
						annot_assembly_table = 'genome.annot_assembly', gene_table='genome.gene', \
						gene2go_table='genome.gene2go', tax_id=3702):
		"""
		2008-09-24
			turn gene_id into integer
		2008-08-03
			schema where tables about genes are from is renamed from 'sequence' to 'genome'
		2008-02-02
			get all the necessary info for genes.
			watch, chromosome here is varchar type (because of chromosome X, Y etc)
		"""
		sys.stderr.write("Getting gene_id2model and chr_id2gene_id_ls...")
		from annot.bin.codense.common import pg_1d_array2python_ls
		gene_id2model = {}
		chr_id2gene_id_ls = {}
		curs.execute("DECLARE gene_crs CURSOR FOR select e.gene_id, a.chromosome, e.start, e.stop, e.mrna_start, e.mrna_stop, e.cds_start, e.cds_stop, e.strand, g.gene_symbol, g.description, g.type_of_gene \
					from %s e, %s a, %s g where e.gene_id=g.gene_id and e.genomic_gi=a.gi and e.tax_id=%s order by chromosome, start, stop"%\
					(entrezgene_mapping_table, annot_assembly_table, gene_table, tax_id))
		curs.execute("fetch 5000 from gene_crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				#gene_id is integer. chromosome is varchar.
				gene_id, chromosome, start, stop, mrna_start, mrna_stop, cds_start, cds_stop, strand, symbol, description, type_of_gene = row
				gene_id = int(gene_id)	#2008-09-24
				if cds_start and cds_stop:
					if type(cds_start)!=list:
						cds_start = pg_1d_array2python_ls(cds_start, int)
						cds_stop = pg_1d_array2python_ls(cds_stop, int)
				else:
					cds_start = cds_stop = None
				
				if mrna_start and mrna_stop:
					if type(mrna_stop)!=list:
						mrna_start = pg_1d_array2python_ls(mrna_start, int)
						mrna_stop = pg_1d_array2python_ls(mrna_stop, int)
				else:
					mrna_start = mrna_stop = None
				
				if chromosome not in chr_id2gene_id_ls:
					chr_id2gene_id_ls[chromosome] = []
				chr_id2gene_id_ls[chromosome].append(gene_id)
				if gene_id not in gene_id2model:
					gene_id2model[gene_id] = GeneModel(gene_id, chromosome, symbol, description, type_of_gene, \
														start, stop, mrna_start, mrna_stop, cds_start, cds_stop, strand)
			curs.execute("fetch 5000 from gene_crs")
			rows = curs.fetchall()
		curs.execute("close gene_crs")
		sys.stderr.write("Done.\n")
		return gene_id2model, chr_id2gene_id_ls
	
	get_gene_id2model = classmethod(get_gene_id2model)
	
	def plot_one_gene(self, ax, gene_id, gene_id2model, chr_id2cumu_size, chr_id2size, chr_gap, y_value=1, gene_width=1.0):
		"""
		2008-12-16
			defunct. DrawSNPRegion.drawGeneModel() is used in on_button_draw_annotation_clicked()
		2008-02-02
			draw a single gene on the canvas, 
		"""
		gene_model = gene_id2model.get(gene_id)
		if gene_model:
			c_start_ls = None
			c_end_ls = None
			if gene_model.cds_start!=None and gene_model.cds_stop!=None:
				c_start_ls = gene_model.cds_start
				c_end_ls = gene_model.cds_stop
			elif gene_model.mrna_start!=None and gene_model.mrna_stop!=None:
				c_start_ls = gene_model.mrna_start
				c_end_ls = gene_model.mrna_stop
			elif gene_model.start!=None and gene_model.stop!=None:
				c_start_ls = [gene_model.start]
				c_end_ls = [gene_model.stop]
			if c_start_ls and c_end_ls:
				chromosome = gene_model.chromosome
				this_chr_starting_pos_on_plot = chr_id2cumu_size[chromosome]-chr_id2size[chromosome]-chr_gap
				if gene_model.strand=="1":
					g_artist = Gene(c_start_ls, c_end_ls, y=y_value, x_offset=this_chr_starting_pos_on_plot, \
								width=gene_width, alpha=0.3, facecolor='r', picker=True)
				elif gene_model.strand=="-1":	#to draw opposite strand, 1st is to order c_start_ls and c_end_ls 
					# in descending order. 2nd is to swap c_start_ls and c_end_ls.
					#c_start_ls.reverse()	#2008-02-04 it's already in descending order in db.
					#c_end_ls.reverse()	#2008-02-04 it's already in descending order in db.
					g_artist = Gene(c_end_ls, c_start_ls, y=y_value, x_offset=this_chr_starting_pos_on_plot, \
								width=gene_width, alpha=0.3, facecolor='r', picker=True)
				else:	#no arrow
					g_artist = Gene(c_start_ls, c_end_ls, y=y_value, is_arrow=False, x_offset=this_chr_starting_pos_on_plot, \
								width=gene_width, alpha=0.3, facecolor='r', picker=True)
				ax.add_artist(g_artist)
				artist_obj_id = id(g_artist)
				self.artist_obj_id2artist_gene_id_ls[artist_obj_id] = [g_artist, gene_id]
				self.gene_id2artist_object_id[gene_id] = artist_obj_id
	
	def on_button_draw_annotation_clicked(self, widget, data=None):
		"""
		2008-12-16
			use DrawSNPRegion.drawGeneModel() to draw gene models
		2008-02-02
		"""
		if not self.chr_id2size:
			sys.stderr.write("No genome-wide pvalue plot has been drawn yet. Do it first!\n")
			return
		#if not self.gene_id2model:
		#	self.gene_id2model, self.chr_id2gene_id_ls = self.get_gene_id2model(self.postgres_curs, tax_id=3702)
		if not self.gene_annotation:
			self.db_connect()
		
		xlim = self.axe_gene_model.get_xlim()
		left_chr, left_pos = get_chr_pos_from_x_axis_pos(xlim[0], self.chr_gap, self.chr_id2cumu_size, self.chr_id_ls)
		right_chr, right_pos = get_chr_pos_from_x_axis_pos(xlim[1], self.chr_gap, self.chr_id2cumu_size, self.chr_id_ls)
		
		#fake a snps_within_this_region for drawGeneModel()
		snps_within_this_region = PassingData(chr_pos_ls=[[left_chr, left_pos],[right_chr, right_pos]])
		base_y_value = 1
		gene_width = 0.8
		gene_position_cycle = 5
		
		return_data = DrawSNPRegion.drawGeneModel(self.axe_gene_model, snps_within_this_region, self.gene_annotation, candidate_gene_set=None,\
								gene_width=gene_width, gene_position_cycle=gene_position_cycle, base_y_value=base_y_value, \
								gene_box_text_gap=20, label_gene=0, rotate_xy=False,\
								chr_id2cumu_size=self.chr_id2cumu_size, chr_id2size=self.chr_id2size, chr_gap=self.chr_gap,\
								artist_obj_id2artist_gene_id_ls=self.artist_obj_id2artist_gene_id_ls, \
								gene_id2artist_object_id=self.gene_id2artist_object_id, drawGeneOnTheBoundary=False)
					#set drawGeneOnTheBoundary to False because later adding text to these genes would corrupt the running program.
		self.axe_gene_model.set_ylim([base_y_value-gene_width, gene_position_cycle+gene_width*2])
		
		"""
		for gene_id in self.chr_id2gene_id_ls[left_chr]:
			gene_model = self.gene_id2model[gene_id]
			if gene_model.start!=None and gene_model.stop!=None and gene_model.stop>left_pos and gene_id not in self.gene_id2artist_object_id:
				if left_chr==right_chr:	#same chromosome
					if gene_model.start>right_pos:	#totally out of range, skip it
						continue
				y_value = len(self.gene_id2artist_object_id)%4	#cycling through the y position to avoid clogging
				self.plot_one_gene(self.ax, gene_id, self.gene_id2model, self.chr_id2cumu_size, self.chr_id2size, self.chr_gap, y_value=-1-y_value, gene_width=self.gene_width)
		if left_chr!=right_chr:
			for gene_id in self.chr_id2gene_id_ls[right_chr]:
				gene_model = self.gene_id2model[gene_id]
				if gene_model.start!=None and gene_model.stop!=None and gene_model.start<right_pos and gene_id not in self.gene_id2artist_object_id:
					y_value = len(self.gene_id2artist_object_id)%4	#cycling through the y position to avoid clogging
					self.plot_one_gene(self.ax, gene_id, self.gene_id2model, self.chr_id2cumu_size, self.chr_id2size, self.chr_gap, y_value=-1-y_value, gene_width=self.gene_width)
		"""
		self.canvas_matplotlib.draw()
	
	def on_imagemenuitem_preferences_activate(self, event, data=None):
		"""
		2008-02-04
		"""
		self.dialog_preferences.show_all()
	
	def on_button_dialog_preferences_ok_clicked(self, widget, data=None):
		"""
		2008-02-04
			change some preferences
		"""
		self.dialog_preferences.hide()
		if self.checkbutton_debug.get_active():
			self.debug = 1
		else:
			self.debug = 0
		if self.checkbutton_stderr.get_active():
			sys.stderr = self.dummy_err
		else:
			sys.stderr = sys.__stderr__
		if self.checkbutton_stdout.get_active():
			sys.stdout = self.dummy_out
		else:
			sys.stdout = sys.__stdout__
		if self.checkbutton_draw_gene_symbol.get_active():
			self.draw_gene_symbol_when_clicked = 1
		else:
			self.draw_gene_symbol_when_clicked = 0
		self.gene_width = float(self.entry_gene_width.get_text())
	
	def on_button_dialog_preferences_cancel_clicked(self, widget, data=None):
		"""
		2008-02-04
			don't change any preferences
		"""
		self.dialog_preferences.hide()
	
	def on_imagemenuitem_about_activate(self, widget):
		"""
		2008-02-04
		"""
		self.aboutdialog1.show_all()
	
	def on_imagemenuitem_cleanup_output_activate(self, widget):
		"""
		2008-02-04
			clean up output buffer
		"""
		self.textbuffer_output.set_text('')
	
	def on_imagemenuitem_clear_plot_activate(self, widget):
		"""
		2010-11-22
			clear the axes (not the axe_gene_model)
		"""
		self.ax.clear()
		self.genome_wide_results.clear()
	
	def on_checkbutton_debug_toggled(self, widget):
		"""
		2008-05-28
		"""
		if self.checkbutton_debug.get_active():
			self.debug = 1
		else:
			self.debug = 0
	
	def on_checkbutton_stdout_toggled(self, widget):
		"""
		2008-05-28
		"""
		if self.checkbutton_stdout.get_active():
			sys.stdout = self.dummy_out
		else:
			sys.stdout = sys.__stdout__
	
	def on_checkbutton_stderr_toggled(self, widget):
		"""
		2008-05-28
		"""
		if self.checkbutton_stderr.get_active():
			sys.stderr = self.dummy_err
		else:
			sys.stderr = sys.__stderr__
	
	def on_checkbutton_draw_gene_symbol_toggled(self, widget):
		"""
		2008-05-28
		"""
		if self.checkbutton_draw_gene_symbol.get_active():
			self.draw_gene_symbol_when_clicked = 1
		else:
			self.draw_gene_symbol_when_clicked = 0
	
	def on_entry_gene_width_changed(self, widget):
		"""
		2008-05-28
		"""
		self.gene_width = float(self.entry_gene_width.get_text())
	
	def on_filechooserbutton_gene_annot_file_set(self, widget):
		"""
		2008-12-16
		"""
		self.entry_gene_annotation_picklef.set_text(self.filechooserbutton_gene_annot.get_filename())
	
	def on_button_draw_gene_list_bars_clicked(self, widget):
		"""
		2008-12-16
			draw vertical spans to denote the locations of genes from a candidate list
		"""
		if self.db is None:
			self.db_connect()
		if not self.chr_id2size:
			sys.stderr.write("No genome-wide pvalue plot has been drawn yet. Do it first!\n")
			return
		entry_gene_list_id = self.xml.get_object("entry_gene_list_id")
		list_type_id = entry_gene_list_id.get_text()
		comboboxentry_bar_color = self.xml.get_object("comboboxentry_bar_color")
		bar_color = comboboxentry_bar_color.get_active_text()
		if not bar_color:	#default is black
			bar_color = 'k'
		if list_type_id:
			list_type_id = int(list_type_id)
			self.candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
			for gene_id in self.candidate_gene_set:
				gene_model = self.gene_annotation.gene_id2model[gene_id]
				if gene_id in self.gene_id2vspan_obj_id:
					artist_obj_id = self.gene_id2vspan_obj_id[gene_id]
					artist = self.artist_obj_id2artist_gene_id_ls[artist_obj_id][0]
					if artist.get_edgecolor()!=bar_color:
						artist.set_edgecolor(bar_color)
					if artist.get_facecolor()!=bar_color:
						artist.set_facecolor(bar_color)
					#artist.remove()
				else:
					this_chr_starting_pos_on_plot = self.chr_id2cumu_size[gene_model.chromosome]-\
							self.chr_id2size[gene_model.chromosome]-self.chr_gap
					xmin = this_chr_starting_pos_on_plot + gene_model.start
					xmax = this_chr_starting_pos_on_plot + gene_model.stop
					artist = self.ax.axvspan(xmin, xmax, edgecolor=bar_color, facecolor=bar_color, alpha=0.3, picker=6)
					artist_obj_id = id(artist)
					self.artist_obj_id2artist_gene_id_ls[artist_obj_id] = [artist, gene_id]
					self.gene_id2vspan_obj_id[gene_id] = artist_obj_id
			self.canvas_matplotlib.draw()
	
	def on_button_adjust_gene_axis_clicked(self, widget):
		"""
		2008-12-19
			sometimes after zoom-in/out, axe_gene_model loses track of its y-range and the gene models in it float into ax.
			this function would bring the y-range of axe_gene_model into normal range.
		"""
		base_y_value = 1
		gene_width = 0.8
		gene_position_cycle = 5
		self.axe_gene_model.set_ylim([base_y_value-gene_width, gene_position_cycle+gene_width*2])
		self.canvas_matplotlib.draw()
	
	def on_comboboxentry_which_cnv_data_from_db_changed(self, widget, data=None):
		"""
		2010-6-14
			deal with hbox_seq_ref_pos_version
		2010-3-17
		"""
		type_of_cnv_data = self.comboboxentry_which_cnv_data_from_db.get_active()
		entry_cnv_type_id = self.xml.get_object("entry_cnv_type_id")
		hbox_seq_ref_pos_version = self.xml.get_object('hbox_seq_ref_pos_version')
		#2010-7-26
		hbox_cnv_type_cnvqc_db = self.xml.get_object('hbox_cnv_type_cnvqc_db')
		hbox_cnv_method_cnvqc_db = self.xml.get_object('hbox_cnv_method_cnvqc_db')
		entry_cnv_method_id = self.xml.get_object("entry_cnv_method_id")
		
		if type_of_cnv_data==0:	# CNVQCCall
			self.checkbutton_within_coverage.set_sensitive(True)
			self.radiobutton_non_ecotype_id.set_label("accession id:")
		elif type_of_cnv_data==1:	# sequence fragments of one accession
			self.checkbutton_within_coverage.set_sensitive(False)
			self.radiobutton_non_ecotype_id.set_label("accession id:")
			self.radiobutton_ecotype_id.set_sensitive(False)
		elif type_of_cnv_data == 2 or type_of_cnv_data == 3:	# one sequence fragment
			self.checkbutton_within_coverage.set_sensitive(False)
			self.radiobutton_non_ecotype_id.set_label("sequence fragment id:")
			self.radiobutton_ecotype_id.set_sensitive(False)
			"""
			self.entry_cnv_non_ecotype_id.append_text("16719")
			self.entry_cnv_non_ecotype_id.append_text("16720")
			"""
		elif type_of_cnv_data == 4:	# CNVCall
			self.checkbutton_within_coverage.set_sensitive(True)
			self.radiobutton_non_ecotype_id.set_label("array id:")
			self.radiobutton_ecotype_id.set_sensitive(True)
		elif type_of_cnv_data == 5:	# CNV	#2010-8-5
			self.checkbutton_within_coverage.set_sensitive(True)
			self.radiobutton_ecotype_id.set_sensitive(False)
		elif type_of_cnv_data == 6:	# CNVArrayCall
			self.checkbutton_within_coverage.set_sensitive(True)
			self.radiobutton_non_ecotype_id.set_label("array id:")
			self.radiobutton_ecotype_id.set_sensitive(False)
		
		if type_of_cnv_data==5:	#2010-8-5 CNV doesn't need any ID
			self.radiobutton_non_ecotype_id.set_sensitive(False)
		else:
			self.radiobutton_non_ecotype_id.set_sensitive(True)
		
		data_type_id2CNVTable = {0: Stock_250kDB.CNVQCCall, 4:Stock_250kDB.CNVCall, 5: Stock_250kDB.CNV, 6:Stock_250kDB.CNVArrayCall}
		if type_of_cnv_data in [0, 4, 5, 6]:	# CNVQCCall or CNVCall
			CNVTableClass = data_type_id2CNVTable[type_of_cnv_data]
			
			hbox_cnv_method_cnvqc_db.set_sensitive(True)
			cnv_method_info = self.db.getCNVMethodOrTypeInfoDataInCNVCallOrQC(CNVTableClass=CNVTableClass)
			self.fillComboBox(entry_cnv_method_id, cnv_method_info.list_label_ls)
			
			if type_of_cnv_data==6:	#2010-8-5 CNVArrayCall doesn't have cnv type info 
				hbox_cnv_type_cnvqc_db.set_sensitive(False)
				self.fillComboBox(entry_cnv_type_id, [])
			else:
				hbox_cnv_type_cnvqc_db.set_sensitive(True)
				cnv_type_info = self.db.getCNVMethodOrTypeInfoDataInCNVCallOrQC(TableClass=Stock_250kDB.CNVType, \
															CNVTableClass=CNVTableClass)
				self.fillComboBox(entry_cnv_type_id, cnv_type_info.list_label_ls)
			
			if type_of_cnv_data==4:	#CNVCall
				# fill the comboboxentry_id_cnv_by_region with ids from array_info that have data in CNVCall
				array_id_info = self.db.getArrayIDInfoWithDataInGivenTable()
				self.fillComboBox(self.entry_cnv_non_ecotype_id, array_id_info.list_label_ls)
			elif type_of_cnv_data==0:	#CNVQCCall
				accession_id_info = self.db.getAccessionIDInfoWithDataInCNVQCCall()
				self.fillComboBox(self.entry_cnv_non_ecotype_id, accession_id_info.list_label_ls)
			elif type_of_cnv_data==6:	#CNVArrayCall
				accession_id_info = self.db.getArrayIDInfoWithDataInGivenTable(TableClass=Stock_250kDB.CNVArrayCall)
				self.fillComboBox(self.entry_cnv_non_ecotype_id, accession_id_info.list_label_ls)
		else:
			# fill the comboboxentry_id_cnv_by_region with nothing
			self.fillComboBox(self.entry_cnv_non_ecotype_id, [])
			
			hbox_cnv_type_cnvqc_db.set_sensitive(False)
			self.fillComboBox(entry_cnv_type_id, [])
			hbox_cnv_method_cnvqc_db.set_sensitive(False)
			self.fillComboBox(entry_cnv_method_id, [])
			
			
		#2010-6-14
		comboboxentry_seq_ref_pos_version = self.xml.get_object("comboboxentry_seq_ref_pos_version")
		if type_of_cnv_data == 1 or type_of_cnv_data == 2:	# SequenceFragmentRefPos
			hbox_seq_ref_pos_version.set_sensitive(True)
			version_info = self.db.getSeqFragmentRefPosVersionInfo()
			self.fillComboBox(comboboxentry_seq_ref_pos_version, version_info.list_label_ls)
		else:
			hbox_seq_ref_pos_version.set_sensitive(False)
			self.fillComboBox(comboboxentry_seq_ref_pos_version, [])

	def on_comboboxentry_cnv_by_region_changed(self, widget, data=None):
		"""
		2010-4-27
			cnv_by_region_data_type_ls = ["Probe Intensity", "SequenceFragment Probe-Blast Result", \
									"Reference Coverage By SequenceFragment", "CNVQCCall", "CNVCall"]
		"""
		cnv_by_region_data_type = widget.get_active()
		comboboxentry_id_cnv_by_region = self.xml.get_object("comboboxentry_id_cnv_by_region")
		
		hbox_intensity_fname = self.xml.get_object("hbox_intensity_fname")
		if cnv_by_region_data_type==0 or cnv_by_region_data_type==6:
			hbox_intensity_fname.set_sensitive(True)
		else:
			hbox_intensity_fname.set_sensitive(False)
		
		# 2010-11-23 determine which widgets to turn off when type 6 is chosen
		hbox_combox_id_entry = self.xml.get_object("hbox_combox_id_entry")
		hbox_chr_cnv_by_region = self.xml.get_object("hbox_chr_cnv_by_region")
		hbox_entry_min_size_by_region = self.xml.get_object("hbox_entry_min_size_by_region")
		hbox_min_no_of_probes_by_region = self.xml.get_object("hbox_min_no_of_probes_by_region")
		label_intensity_fname = self.xml.get_object("label_intensity_fname")
		if cnv_by_region_data_type == 6:	#100bp Coverage by PESolexaData
			hbox_combox_id_entry.set_sensitive(False)
			hbox_chr_cnv_by_region.set_sensitive(False)
			hbox_entry_min_size_by_region.set_sensitive(False)
			hbox_min_no_of_probes_by_region.set_sensitive(False)
			label_intensity_fname.set_text("coverage file: ")
		else:
			hbox_combox_id_entry.set_sensitive(True)
			hbox_chr_cnv_by_region.set_sensitive(True)
			hbox_entry_min_size_by_region.set_sensitive(True)
			hbox_min_no_of_probes_by_region.set_sensitive(True)
			label_intensity_fname.set_text("intensity file: ")
		
		if cnv_by_region_data_type==0:	#Probe Intensity
			self.combobox_id_type_cnv_by_region.set_active(0)
			# fill the comboboxentry_id_cnv_by_region with nothing
			self.fillComboBox(comboboxentry_id_cnv_by_region, [])
		elif cnv_by_region_data_type==1:	# SequenceFragment Probe-Blast Result
			self.combobox_id_type_cnv_by_region.set_active(2)
			# fill the comboboxentry_id_cnv_by_region with ids from CNVQCAccession that have data in SequenceFragment2Probe
			accession_id_info = self.db.getAccessionIDInfoWithDataInSequenceFragment2Probe()
			self.fillComboBox(comboboxentry_id_cnv_by_region, accession_id_info.list_label_ls)
		elif cnv_by_region_data_type==2:	# Reference Coverage By SequenceFragment
			self.combobox_id_type_cnv_by_region.set_active(2)
			# fill the comboboxentry_id_cnv_by_region with ids from CNVQCAccession that have data in SequenceFragmentRefPos
			accession_id_info = self.db.getAccessionIDInfoWithDataInSequenceFragmentRefPos()
			self.fillComboBox(comboboxentry_id_cnv_by_region, accession_id_info.list_label_ls)
		elif cnv_by_region_data_type == 3:	# CNVQCCall
			self.combobox_id_type_cnv_by_region.set_active(2)
			# fill the comboboxentry_id_cnv_by_region with ids from CNVQCAccession that have data in CNVQCCall
			accession_id_info = self.db.getAccessionIDInfoWithDataInCNVQCCall()
			self.fillComboBox(comboboxentry_id_cnv_by_region, accession_id_info.list_label_ls)
		elif cnv_by_region_data_type==4:	# CNVCall
			self.combobox_id_type_cnv_by_region.set_active(0)
			# fill the comboboxentry_id_cnv_by_region with ids from array_info that have data in CNVCall
			array_id_info = self.db.getArrayIDInfoWithDataInGivenTable()
			self.fillComboBox(comboboxentry_id_cnv_by_region, array_id_info.list_label_ls)
		
		elif cnv_by_region_data_type == 5:	#Ref Coverage By Sequences (2D)
			self.combobox_id_type_cnv_by_region.set_active(2)
			accession_id_info = self.db.getAccessionIDInfoWithDataInSequenceFragmentRefPos()
			self.fillComboBox(comboboxentry_id_cnv_by_region, accession_id_info.list_label_ls)
		
		#2010-6-17
		hbox_seq_ref_pos_version_by_region = self.xml.get_object('hbox_seq_ref_pos_version_by_region')
		comboboxentry_seq_ref_pos_version_by_region = self.xml.get_object("comboboxentry_seq_ref_pos_version_by_region")
		if cnv_by_region_data_type == 2 or cnv_by_region_data_type==5:	# SequenceFragmentRefPos
			hbox_seq_ref_pos_version_by_region.set_sensitive(True)
			version_info = self.db.getSeqFragmentRefPosVersionInfo()
			self.fillComboBox(comboboxentry_seq_ref_pos_version_by_region, version_info.list_label_ls)
		else:
			hbox_seq_ref_pos_version_by_region.set_sensitive(False)
			self.fillComboBox(comboboxentry_seq_ref_pos_version_by_region, [])
	
		#2010-7-23
		hbox_cnv_type_id_by_region = self.xml.get_object('hbox_cnv_type_id_by_region')
		comboboxentry_cnv_type_id_by_region = self.xml.get_object("comboboxentry_cnv_type_id_by_region")
		hbox_cnv_method_id_by_region = self.xml.get_object('hbox_cnv_method_id_by_region')
		comboboxentry_cnv_method_id_by_region = self.xml.get_object("comboboxentry_cnv_method_id_by_region")
		
		data_type_id2CNVTable = {3: Stock_250kDB.CNVQCCall, 4:Stock_250kDB.CNVCall}
		if cnv_by_region_data_type == 3 or cnv_by_region_data_type==4:	# CNVQCCall or CNVCall
			CNVTableClass = data_type_id2CNVTable[cnv_by_region_data_type]
			hbox_cnv_type_id_by_region.set_sensitive(True)
			cnv_type_info = self.db.getCNVMethodOrTypeInfoDataInCNVCallOrQC(TableClass=Stock_250kDB.CNVType, \
															CNVTableClass=CNVTableClass)
			self.fillComboBox(comboboxentry_cnv_type_id_by_region, cnv_type_info.list_label_ls)
			
			hbox_cnv_method_id_by_region.set_sensitive(True)
			cnv_method_info = self.db.getCNVMethodOrTypeInfoDataInCNVCallOrQC(CNVTableClass=CNVTableClass)
			self.fillComboBox(comboboxentry_cnv_method_id_by_region, cnv_method_info.list_label_ls)
		else:
			hbox_cnv_type_id_by_region.set_sensitive(False)
			self.fillComboBox(comboboxentry_cnv_type_id_by_region, [])
			hbox_cnv_method_id_by_region.set_sensitive(False)
			self.fillComboBox(comboboxentry_cnv_method_id_by_region, [])
	
	def on_radiobutton_cnv_db_id_toggled(self, widget, data=None):
		"""
		2010-3-17
		"""
		if widget.get_active() == 1:	#only change self.submit_option to the active radiobutton
			if data=='radiobutton_non_ecotype_id':
				self.entry_cnv_non_ecotype_id.set_sensitive(True)
				self.entry_cnv_ecotype_id.set_sensitive(False)
			elif data=='radiobutton_ecotype_id':
				self.entry_cnv_non_ecotype_id.set_sensitive(False)
				self.entry_cnv_ecotype_id.set_sensitive(True)
	
	def on_filechooserbutton_CNV_intensity_file_set(self, widget, data=None):
		"""
		2010-3-17
		"""
		entry_cnv_intensity_fname = self.xml.get_object("entry_cnv_intensity_fname")
		cnv_intensity_fname = widget.get_filename()
		entry_cnv_intensity_fname.set_text(cnv_intensity_fname)
	
	def on_filechooserbutton_intensity_fname_file_set(self, widget, data=None):
		"""
		2010-4-27
		"""
		entry_intensity_fname = self.xml.get_object("entry_intensity_fname")
		cnv_intensity_fname = widget.get_filename()
		entry_intensity_fname.set_text(cnv_intensity_fname)
		
	def on_combobox_type_of_nearby_data_changed(self, widget, data=None):
		"""
		2010-7-28
			renamed from on_checkbutton_plot_cnv_intensity_toggled().
			check button got replaced by a combo box
		2010-3-17
		"""
		hbox_cnv_intensity = self.xml.get_object("hbox_cnv_intensity")
		if widget.get_active()>0:
			hbox_cnv_intensity.set_sensitive(True)
		else:
			hbox_cnv_intensity.set_sensitive(False)
		
		
if __name__ == '__main__':
	#prog = gnome.program_init('GenomeBrowser', '0.1')	# 2010-3-14 no more gnome part.
	instance = GenomeBrowser()
	gtk.main()
