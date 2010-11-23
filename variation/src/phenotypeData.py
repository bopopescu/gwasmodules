import pdb
from env import *
import math
import scipy as sp
import itertools as it

# Phenotype categories: (category,order)
#phenotypeCategories = {
#				1: (1,1), 2: (1,2), 3: (1,3), 4: (1,4), 39: (1,5), 40: (1,6), 41: (1,7), 42: (1,8), 
#				5: (1,9), 6: (1,10), 7: (1,11), 80: (1,12), 81: (1,13), 82: (1,14),	45: (1,15), 
#				46: (1,16), 47: (1,17), 48: (1,18), 57: (1,19), 59: (1,20), 58: (1,21), 43: (1,22), 44: (1,23),
#					
#				60: (2,1), 61: (2,2), 62: (2,3), 63: (2,4), 64: (2,5), 8: (2,6), 161: (2,7), 162: (2,8), 163: (2,9), 182: (2,10),
#				164: (2,11), 165: (2,12), 166: (2,13), 173: (2,14), 174: (2,15), 175: (2,16), 176: (2,17), 177: (2,18), 178: (2,19), 179: (2,19),
#				158: (2,19), 159: (2,19), 75: (2,20), 76: (2,21), 77: (2,22), 78: (2,23), 79: (2,24), 
#				167: (2,25), 168: (2,26), 169: (2,27), 170: (2,28), 171: (2,29), 172: (2,30), 183: (2,31), 184: (2,32),
#				272: (2,33),273: (2,34),274: (2,35),277: (2,36),278: (2,37),279: (2,38),280: (2,39),281: (2,40),282: (2,41),
#					
#				32: (3,1), 33: (3,2), 35: (3,3), 34: (3,4),
#				9: (3,5), 10: (3,6), 11: (3,7), 12: (3,8), 13: (3,9),
#				65: (3,10), 67: (3,11), 69: (3,12), 71: (3,13), 73: (3,14), 
#				66: (3,15), 68: (3,16), 70: (3,17), 72: (3,18), 74: (3,19),
#				186: (3,20), 185: (3,21),
#					
#				14: (4,1), 15: (4,2), 16: (4,3), 17: (4,4), 18: (4,5), 19: (4,6), 20: (4,7), 21: (4,8), 22: (4,9), 23: (4,10), 
#				24: (4,11), 25: (4,12), 26: (4,13), 27: (4,14), 28: (4,15), 29: (4,16), 30: (4,17), 31: (4,18)
#			}
#
#categories_2_phenotypes = {1:[1,2,3,4,39,40,41,42,5,6,7,80,81,82,45,46,47,48,57,59,58,43,44],
#				2:[60,61,62,63,64,8,161,162,163,182,164,165,166,173,174,175,176,177,178,179,158,159,
#					75,76,77,78,79,167,168,169,170,171,172,272,273,274,277,278,279,280,281,282],
#				3:[32,33,35,34,9,10,11,12,13,65,67,69,71,73,66,68,70,72,74,186,185,183,184],
#				4:[14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]}
#
#categories_2_phenotypes = {1:[5,80,6,81,7,82,1,2,47,48,45,46,59,57,39,40,41,42,3,4,43,44,58],
#categories_2_phenotypes = {1:[1,2,3,4,39,40,41,42,5,6,7,80,81,82,45,46,47,48,57,59,58,43,44],
#				2:[60,61,62,63,64,8,161,162,163,182,164,165,166,173,174,175,176,177,178,179,158,159,
#					75,76,77,78,79,167,168,169,170,171,172,272,273,274,277,278,279,280,281,282,283],
#				3:[32,33,35,34,9,10,11,12,13,65,67,69,71,73,66,68,70,72,74,186,185,183,184],
#				4:[14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]}
#
#Log 02/25/09 - Chlor_16,Chlor_22 were removed from cat. 2: 168,169, { 168: (2,26), 169: (2,27),} 
#new_name_dict = {11:"Emoy2", 32:"Emoy2", }


class phenotype_data:
	"""
	A class that encapsulates phenotype values and provides basic functionality for these.
	
	This is an update of an older class.
	"""

	def __init__(self, phen_dict=None, phenotype_names=[], phen_ids=None):
		if phen_dict:
			self.phen_dict = phen_dict
			self.phen_ids = phen_dict.keys()
		else:
			if phen_ids:
				self.phen_ids = phen_ids
			else:
				self.phen_ids = range(len(phenotype_names))
			self.phen_dict = {}
			for i, pid in enumerate(self.phen_ids):
				self.phen_dict[pid] = {'name':phenotype_names[i], 'ecotypes':[], 'values':[]}


	def log_transform(self, pid, method='standard'):
		a = sp.array(self.phen_dict[pid]['values'])
		if method == 'standard':
			vals = sp.log((a - min(a)) + 0.1 * sp.std(a))
		else:
			vals = sp.log(a)
		self.phen_dict[pid]['values'] = vals.tolist()

	def sqrt_transform(self, pid, method='standard'):
		a = sp.array(self.phen_dict[pid]['values'])
		if method == 'standard':
			vals = sp.sqrt((a - min(a)) + 0.1 * sp.std(a))
		else:
			vals = sp.sqrt(a)
		self.phen_dict[pid]['values'] = vals.tolist()


	def na_outliers(self, pids, iqr_threshold):
		raise NotImplementedError


	def filter_phenotypes(self, pids_to_keep):
		"""
		Removes phenotypes.
		"""
		self.phen_ids = pids_to_keep
		d = {}
		for pid in pids_to_keep:
			d[pid] = self.phen_dict[pid]
		self.phen_dict = d


	def filter_ecotypes(self, indices_to_keep, pids=None):
		"""
		Removes the ecotypes from all data.
		"""
		if not pids:
			pids = self.phen_ids
		for pid in pids:
			el = []
			vl = []
			d = self.phen_dict[pid]
			for i in indices_to_keep:
				el.append(d['ecotypes'][i])
				vl.append(d['values'][i])
			self.phen_dict[pid]['ecotypes'] = el
			self.phen_dict[pid]['values'] = vl


	def order_ecotypes(self, ets_map, pids=None):
		if not pids:
			pids = self.phen_ids
		for pid in pids:
			d = self.phen_dict[pid]
			ets = []
			vals = []
			for i in ets_map:
				ets.append(d['ecotypes'][i])
				vals.append(d['values'][i])

	def get_name(self, pid):
		return self.phen_dict[pid]['name']

	def get_values(self, pid):
		return self.phen_dict[pid]['values']

	def get_ecotypes(self, pid):
		return self.phen_dict[pid]['ecotypes']

	def get_incidence_matrix(self, pid):
		ets = sp.array(map(int, self.phen_dict[pid]['ecotypes']))
		unique_ets = []
		i = 0
		while i < len(ets):
			et = ets[i]
			unique_ets.append(et)
			while i < len(ets) and ets[i] == et:
				i += 1
#		unique_ets = sp.mat(sp.unique(ets))
		Z = sp.int16(sp.mat(ets).T == sp.mat(unique_ets))
		#print Z
		return Z


	def _get_ecotype_value_dict_(self, pid):
		ecotypes = self.get_ecotypes(pid)
		values = self.get_values(pid)
		d = {}
		for et in set(ecotypes):
			d[et] = {'values':[], 'rep_num':0}

		for et, val in it.izip(ecotypes, values):
			d[et]['values'].append(val)
			d[et]['rep_num'] += 1
		return d

	def get_avg_values(self, pid):
		"""
		Returns the average values, along with the ecotypes, and rep_number
		"""
		d = self._get_ecotype_value_dict_(pid)
		ecotypes = d.keys()
		avg_vals = []
		rep_nums = []
		for et in d:
			avg_vals.append(sp.mean(d[et]['values']))
			rep_nums.append(d[et]['rep_num'])
		return {'name':self.get_name(pid) + '_avg', 'ecotypes':ecotypes, 'values':avg_vals, 'rep_nums': rep_nums}




	def get_db_pid(self, pid, conn=None):
		"""
		Retrieves the DB pid, using the phenotype name.
		"""
		phen_name = self.get_name(pid)
		import dbutils
		if not conn:
			conn = dbutils.connect_to_default_lookup()
		cursor = conn.cursor()
		sql_statement = "SELECT id FROM stock_250k.phenotype_method WHERE short_name='%s'" % (phen_name)
		print sql_statement
		numRows = int(cursor.execute(sql_statement))
		row = cursor.fetchone()
		if row:
			db_pid = int(row[0])
		else:
			print "No id found in DB for phenotype:%s" % (phen_name)
			db_pid = None
		cursor.close()
		if not conn:
			conn.close()
		return db_pid


	def convert_to_averages(self, pids=None):
		"""
		Replaces phenotypes with replicates with their averages.
		"""
		if not pids:
			pids = self.phen_dict.keys()
		for pid in pids:
			phen_name = self.get_name(pid)
			self.phen_dict[pid] = self.get_avg_values(pid)
			self.phen_dict[pid]['name'] = phen_name



	def plot_histogram(self, pid, title=None , pdf_file=None, png_file=None, x_label=None):
		import pylab as plt
		plt.figure(figsize=(5, 4))
		plt.axes([0.14, 0.13, 0.81, 0.83])
		if x_label:
			plt.xlabel(x_label)
		phen_vals = self.get_values(pid)

		minVal = min(phen_vals)
		maxVal = max(phen_vals)
		x_range = maxVal - minVal
		histRes = plt.hist(phen_vals, bins=round(8 + 2 * sp.log(len(phen_vals))), alpha=0.8)
		y_max = max(histRes[0])
		plt.axis([minVal - 0.035 * x_range, maxVal + 0.035 * x_range, -0.035 * y_max, 1.18 * y_max])
		num_phen_vals = len(phen_vals)
		shapiro_pval = sp.stats.shapiro(phen_vals)[1]
		plt.text(maxVal - 0.7 * x_range, 1.1 * y_max, "Number of values: " + str(num_phen_vals), size="x-small")
		plt.text(maxVal - 0.85 * x_range, 1.02 * y_max, "Shapiro-Wilk normality $p$-value: %0.6f" % shapiro_pval , size="x-small")
		print max(histRes[0])
		plt.ylabel("Frequency")
		if title:
			plt.title(title)
		if pdf_file:
			plt.savefig(pdf_file, format="pdf")
		if png_file:
			plt.savefig(png_file, format="png", dpi=300)
		elif not pdf_file or png_file:
			plt.show()
		plt.clf()


	def write_to_file(self, file_name, delim=','):
		"""
		Writes the object to a file.. (in the new format)
		"""
		f = open(file_name, 'w')
		header = ['phenotype_id', 'phenotype_name', 'ecotype_id', 'value', 'replicate_id', ]
		f.write(delim.join(header) + '\n')
		for pid in self.phen_dict:
			d = self.phen_dict[pid]
			phen_name = d['name']
			ets_vals = zip(d['ecotypes'], d['values'])
			ets_vals.sort()
			last_et = -1
			for (et, val) in ets_vals:
				if et != last_et:
					repl_id = 1
				else:
					repl_id += 1
				l = map(str, [pid, phen_name, et, val, repl_id])
				f.write(delim.join(l) + '\n')
				last_et = et
		f.close()


	def is_binary(self, pid):
		return len(sp.unique(self.phen_dict[pid]['values'])) == 2


	def plot_accession_map(self, pid, ecotypes=None, pdf_file=None, png_file=None):
		"""
		Plot accessions on a map.
		
		'color_by' is by default set to be the phenotype values.
		"""
		if not ecotypes:
			ecotypes = self.phen_dict[pid]['ecotypes']
		eid = get_ecotype_id_info_dict()
		lats = []
		lons = []
		for e in ecotypes:
			r = eid[int(e)]
			try:
				latitude = float(r[2])
				longitude = float(r[3])

			except Exception, err_str:
				print "Latitude and Longitude, not found?:", err_str
				print 'Placing them in the Atlantic.'
				latitude = 40
				longitude = -20
			lats.append(latitude)
			lons.append(longitude)

		from mpl_toolkits.basemap import Basemap
		import numpy as np
		import matplotlib.pyplot as plt
		from pylab import cm
		if map_type == "global2":
			plt.figure(figsize=(14, 12))
			m = Basemap(width=21.e6, height=21.e6, projection='gnom', lat_0=76, lon_0=15)
			m.drawparallels(np.arange(20, 90, 20))
			m.drawmeridians(np.arange(-180, 180, 30))
		elif map_type == 'global':

			plt.figure(figsize=(16, 4))
			plt.axes([0.02, 0.02, 0.96, 0.96])
			m = Basemap(projection='cyl', llcrnrlat=10, urcrnrlat=80,
				    llcrnrlon= -130, urcrnrlon=150, lat_ts=20, resolution='c')
			m.drawparallels(np.arange(20, 90, 20))
			m.drawmeridians(np.arange(-180, 180, 30))
		else:
			raise Exception("map_type is invalid")

		#m.drawmapboundary(fill_color='aqua')
		m.drawcoastlines(zorder=0)
		m.fillcontinents(zorder=1)
		#m.fillcontinents(color='coral',lake_color='aqua')

		xs = []
		ys = []
		for lon, lat in zip(lons, lats):
			x, y = m(*np.meshgrid([lon], [lat]))
			xs.append(x)
			ys.append(y)

		if not color_by:
			color_vals = self.getPhenVals(p_i)
		else:
			color_vals = color_by
		if len(color_vals) != len(accessions):
			raise Exception("accessions and color_by_vals values don't match ! ")
		if not cmap:
			num_colors = len(set(color_vals))
			if num_colors <= 10:
				cmap = cm.get_cmap('jet', num_colors)
			else:
				cmap = cm.get_cmap('jet')
		lws = [0] * len(xs)
		plt.scatter(xs, ys, s=10, linewidths=lws, c=color_vals, cmap=cmap, alpha=0.7, zorder=2)
		plt.colorbar()
		if title:
			plt.title(title)
		if pdf_file:
			plt.savefig(pdf_file, format="pdf")
		if png_file:
			plt.savefig(png_file, format="png")
		if not pdf_file and not png_file:
			plt.show()


	def plot_marker_box_plot(self, pid, marker=None, marker_accessions=None, md=None, m_i=None,
				 pdf_file=None, png_file=None, title=None, marker_score=None, marker_missing_val='NA'):
		"""
		Plots a box plot for the given binary marker and phenotype. 
		"""
		raise NotImplementedError
		marker = list(marker)
		import pylab
		if md and m_i:
			r = self.get_marker_and_phen_vals(phen_i, md, m_i)
			marker = r['marker']
			phen_vals = r['phen_vals']
			missing_val = md.missingVal
		elif marker and marker_accessions:
#			new_marker = []
#			new_phen_vals = []
			phen_vals = self.getPhenVals(phen_i, noNAs=False)
			if len(marker_accessions) != len(phen_vals):
				raise Exception
#			for nt,a in zip(marker,marker_accessions):
#				if a in self.accessions:
#					a_i = self.accessions.index(a)
#					phen_val = phen_vals[a_i]
#					if phen_val!='NA' and nt!='NA':
#						new_marker.append(nt)
#						new_phen_vals.append(phen_val)
#			phen_vals = new_phen_vals
#			marker = new_marker		
			missing_val = marker_missing_val

		nt_counts = [(marker.count(nt), nt) for nt in set(marker) if nt != missing_val]
		nt_counts.sort()
		if len(nt_counts) > 2:
			import warnings
			warnings.warn("More than 2 alleles, box-plot might be wrong?")

		allele_phen_vals = []
		for c, nt in nt_counts:
			allele_phen_vals.append([v for i, v in enumerate(phen_vals) if marker[i] == nt])
			pylab.figure()
		pylab.boxplot(allele_phen_vals)
		minVal = min(phen_vals)
		maxVal = max(phen_vals)
		rangeVal = maxVal - minVal
		pylab.axis([0.2, 2.8, minVal - rangeVal * 0.3, maxVal + rangeVal * 0.3])
		pylab.text(0.4, minVal - 0.15 * rangeVal, "# of obs.: ", color='k')
		pylab.text(0.95, minVal - 0.15 * rangeVal, str(len(allele_phen_vals[0])), color='k')
		if len(allele_phen_vals) > 1:
			pylab.text(1.95, minVal - 0.15 * rangeVal, str(len(allele_phen_vals[1])), color='k')
		if marker_score:
			pylab.text(0.9, maxVal + 0.15 * rangeVal, "-log(p-value)/score: " + str(marker_score), color='k')
		if title:
			pylab.title(title)
		elif md and m_i:
			pylab.title(self.getPhenotypeName(phen_i) + ": chromosome " + str(md.chromosome) + ", position " + str(md.positions[m_i]))
		if pdf_file:
			pylab.savefig(pdf_file, format="pdf")
		if png_file:
			pylab.savefig(png_file, format="png")
		if not pdf_file and not png_file:
			pylab.show()




def parse_phenotype_file(file_name, delim=','):
	"""
	Parses a phenotype file of the new format, and returns a new phenotype_data object.
	"""
	phen_dict = {}
	last_pid = -1
	with open(file_name) as f:
		header = f.next()
		if len(header.split(delim)) != 5:
			raise Exception('Problems with delimiter', delim)
		for line in f:
			l = line.split(delim)
			pid = int(l[0])
			if pid != last_pid:
				if last_pid != -1:
					phen_dict[last_pid] = d
				phen_name = l[1]
				d = {'name':phen_name, 'ecotypes':[], 'values':[]}
			d['ecotypes'].append(l[2])
			d['values'].append(float(l[3]))
			last_pid = pid
	return phenotype_data(phen_dict=phen_dict, phen_ids=phen_dict.keys())







class PhenotypeData:
	"""
	A class that knows how to read simple tsv or csv phenotypefiles and facilitates their interactions with SnpsData objects.
	"""

	def __init__(self, accessions, phenotypeNames, phenotypeValues, accessionNames=None,
		use_phen_name_dict=False, with_db_ids=True):
		self.accessions = accessions
		self.phenotypeNames = phenotypeNames
		self.phenotypeNiceNames = phenotypeNames[:]
		self.phenotypeValues = phenotypeValues # list[accession_index][phenotype_index]
		self.accessionNames = accessionNames
		self.phenIds = []
		i = 0
		new_phenotype_names = []
		for phenName in self.phenotypeNames:
			id = phenName.split("_")
			if with_db_ids:
				try:
					self.phenIds.append(int(id[0]))
					new_phenotype_names.append("_".join(id[1:]))
				except Exception:
					i += 1
					self.phenIds.append(i)
					new_phenotype_names.append(phenName)
			else:
				i += 1
				self.phenIds.append(i)
				new_phenotype_names.append(phenName)

		self.phenotypeNames = new_phenotype_names

		if use_phen_name_dict:
			name_dict = self._getPhenotypeNameDict_()
		else:
			name_dict = {}
		index_map = self._getIndexMapping_()
		for pi in self.phenIds:
			if pi in name_dict:
				p_name = self.phenotypeNiceNames[index_map[pi]]
				if p_name == "11_Emoy*":
					p_name = "11_Emoy_star"
				#print p_name,"->", name_dict[pi]
				self.phenotypeNiceNames[index_map[pi]] = name_dict[pi]
			else:
				phenNameList = self.phenotypeNames[index_map[pi]].split("_")[1:]
				self.phenotypeNiceNames[index_map[pi]] = " ".join(phenNameList)



	def _getPhenotypeNameDict_(self):
		#filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/Phenotype_names_061809.csv"
		filename = "/home/cmbpanfs-01/bvilhjal/data/Phenotype_names_061809.csv"
		#filename = "/home/cmb-01/bvilhjal/Projects/data/Phenotype_names_061809.csv"
		name_dict = {}
		f = open(filename, "r")
		lines = f.readlines()
		for line in lines:
			ls = line.split(",")
			name_dict[int(ls[0])] = ls[1].rstrip()
		return name_dict



	def get_db_pid(self, pid, conn=None):
		"""
		Retrieves the DB pid, using the phenotype name.
		"""
		import dbutils
		if not conn:
			conn = dbutils.connect_to_default_lookup()
		cursor = conn.cursor()
		sql_statement = "SELECT id FROM stock_250k.phenotype_method WHERE short_name='%s'" % (self.getPhenotypeName(pid))
		print sql_statement
		numRows = int(cursor.execute(sql_statement))
		row = cursor.fetchone()
		if row:
			db_pid = int(row[0])
		else:
			print "No id found in DB for phenotype:%s" % (self.getPhenotypeName(pid))
			db_pid = None
		cursor.close()
		if not conn:
			conn.close()
		return db_pid



	def getTree(self):
		index_map = self._getIndexMapping_()
		print index_map
		pi_list = [0] * len(index_map)
		for i in index_map:
			pi_list[index_map[i]] = i
		import rpy, util
		import numpy as np
		import scipy as sp
		import scipy.cluster.hierarchy as hc
		import numpy as np
		phenotypeValues = map(list, zip(*self.phenotypeValues))
		phenvals = []

		for i in range(0, len(phenotypeValues)):
			vals = phenotypeValues[i]
			new_vals = []
			val_sum = 0
			val_count = 0
			for val in vals:
				if val != "NA":
					val_sum += float(val)
					val_count += 1.0

			if val_count > 0:
				val_mean = val_sum / val_count
			else:
				val_mean = 0
				print pi_map[i]
			for val in vals:
				if val == "NA":
					new_vals.append(0)
				else:
					new_vals.append((float(val) - val_mean))
			sd = util.calcSD(new_vals)
			for j in range(0, len(new_vals)):
				new_vals[j] = new_vals[j] / sd
			phenvals.append(new_vals)
		phenvals = np.array(phenvals)
		dist_mat = rpy.r.dist(phenvals)#,method=r("canberra"))
		print dist_mat
		Z = hc.average(dist_mat)
		#print "Z:",Z
		#hc.leaders(Z)
		import pylab
		pylab.figure(figsize=(24, 15))
		pylab.axes([0.03, 0.08, 0.96, 0.91])
		dend_dict = hc.dendrogram(Z, leaf_font_size=7, labels=pi_list)
		pylab.show()



	def getPC(self, pc_num=1):
		"""
		Returns principal component values of of the phenotypes
		"""
		index_map = self._getIndexMapping_()
		print index_map
		pi_list = [0] * len(index_map)
		for i in index_map:
			pi_list[index_map[i]] = i
		import rpy, util
		import numpy as np
		phenotypeValues = map(list, zip(*self.phenotypeValues))
		phenvals = []

		for i in range(0, len(phenotypeValues)):
			vals = phenotypeValues[i]
			new_vals = []
			val_sum = 0
			val_count = 0
			for val in vals:
				if val != "NA":
					val_sum += float(val)
					val_count += 1.0

			if val_count > 0:
				val_mean = val_sum / val_count
			else:
				val_mean = 0
				print pi_map[i]
			for val in vals:
				if val == "NA":
					new_vals.append(0)
				else:
					new_vals.append((float(val) - val_mean))
			sd = util.calcSD(new_vals)
			for j in range(0, len(new_vals)):
				new_vals[j] = new_vals[j] / sd
			phenvals.append(new_vals)
		phenvals = np.transpose(np.array(phenvals))
		print phenvals
		pc = rpy.r.princomp(phenvals)
		pc = zip(list(pc["scores"][pc_num - 1]), pi_list)
		for a in pc:
			print a
		pc.sort()
		print pc
		return pc


	def _getIndexMapping_(self):
		indexMapping = dict()
		for i in range(0, len(self.phenIds)):
			indexMapping[self.phenIds[i]] = i
		return indexMapping


	def getPhenIndex(self, phenId):
		return self.phenIds.index(phenId)

	def getVariance(self, phenId):
		import util
		vals = self.getPhenVals(phenId)
		variance = util.calcVar(vals)
		return variance

	def getMinValue(self, phenId):
		import util
		vals = self.getPhenVals(phenId)
		return min(vals)


	def getPhenVals(self, phenId, asString=False, noNAs=True):
		p_i = self.getPhenIndex(phenId)
		vals = []
		for phenVals in self.phenotypeValues:
			if not (noNAs and phenVals[p_i] == 'NA'):
				if asString or phenVals[p_i] == 'NA':
					vals.append(phenVals[p_i])
				else:
					vals.append(float(phenVals[p_i]))
		return vals

	def get_vals_accs(self, pid):
		"""
		Returns the values and accessions for a particular phenotype.
		(It doesn't return missing values)
		"""
		p_i = self.getPhenIndex(pid)
		vals = []
		accs = []
		for phenVals, acc in zip(self.phenotypeValues, self.accessions):
			if phenVals[p_i] != 'NA':
				vals.append(float(phenVals[p_i]))
				accs.append(acc)
		return (vals, accs)


#	def getAccessionsWithValues(self, phenId):
#		import warnings
#		warnings.warn("Use getNonNAEcotypes instead of getAccessionsWithValues!")
#		p_i = self.getPhenIndex(phenId)
#		accessions = []
#		for i in range(0, len(self.phenotypeValues)):
#			phenVals = self.phenotypeValues[i]
#			acc = self.accessions[i]
#			if phenVals[p_i] != 'NA':
#				accessions.append(acc)
#		return accessions


	def onlyBiologyCategory(self, phenotypeCategory, host=None, user=None, passwd=None):
		import MySQLdb

		#print "Connecting to db, host="+host
		if not user:
			import sys
			sys.stdout.write("Username: ")
			user = sys.stdin.readline().rstrip()
		if not passwd:
			import getpass
			passwd = getpass.getpass()
		try:
			conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="stock_250k")
		except MySQLdb.Error, e:
			print "Error %d: %s" % (e.args[0], e.args[1])
			sys.exit (1)
		cursor = conn.cursor ()
		#Retrieve the filenames
		#print "Fetching biological info on phenotypes."  

		numRows = int(cursor.execute("select distinct id, short_name, biology_category_id from stock_250k.phenotype_method order by id"))

		bioinfo = []
		currTairID = ""
		while(1):
			row = cursor.fetchone()
			if not row:
				break;
			bioinfo.append([int(row[0]), row[1], int(row[2])])

		cursor.close ()
		conn.close ()
		#print "Biol. info fetched"

		phenIds = []
		for phenName in self.phenotypeNames:
			id = phenName.split("_")
			phenIds.append(int(id[0]))

		indicesToKeep = []
		j = 0
		for i in range(0, len(phenIds)):
			while phenIds[i] > bioinfo[j][0]:
				j += 1
			if phenIds[i] == bioinfo[j][0]:
				if bioinfo[j][2] == phenotypeCategory:
					indicesToKeep.append(i)
				j += 1

		self.removePhenotypes(indicesToKeep)


	def removePhenotypes(self, indicesToKeep):
		"""
		Removes phenotypes from the data.
		"""
		numRemoved = len(self.phenotypeNames) - len(indicesToKeep)
		#print "Removing",numRemoved,"phenotypes in phenotype data, out of",len(self.phenotypeNames), "phenotypes."
		newPhenotypeNames = []
		newPhenotVals = []
		newPhenIds = []
		for j in range(0, len(self.accessions)):
			newPhenotVals.append([])
		for i in indicesToKeep:
			newPhenotypeNames.append(self.phenotypeNames[i])
			newPhenIds.append(self.phenIds[i])
			for j in range(0, len(self.accessions)):
				newPhenotVals[j].append(self.phenotypeValues[j][i])
		self.phenIds = newPhenIds
		self.phenotypeValues = newPhenotVals
		self.phenotypeNames = newPhenotypeNames

		#print "len(self.phenotypeNames):",len(self.phenotypeNames)
		#print "len(self.phenotypeValues[0]):",len(self.phenotypeValues[0])


	def removePhenotypeIDs(self, idsToKeep):
		"""
		Removes phenotypes from the data.
		"""
		idMap = self._getIndexMapping_()
		indicesToKeep = []
		for p_i in idsToKeep:
			indicesToKeep.append(idMap[p_i])
		self.removePhenotypes(indicesToKeep)


	def merge_repeated_accessions(self):
		"""
		Merges phenotypes for multiple repeated accessions..
		"""
		raise NotImplementedError #FINISH FUNCTION
		accs = list(set(self.accessions))
		acc_map = []
		for i, acc in enumerate(accs):
			if self.accessions.count(acc) > 1:
				merge_indices = []
				for j, acc2 in enumerate(self.accessions):
					if acc == acc2:
						merge_indices.append(j)
				print "merging: " + ', '.join([self.accessions[mi] for mi in merge_indices])

			else:
				merge_indices = [self.accessions.index(acc)]
			acc_map.append(merge_indices)


	def naOutliers(self, phenotypeID, iqrThreshold=10.0):
		"""
		Removes outliers from the data
		"""
		numNA = 0
		idMap = self._getIndexMapping_()
		#pdb.set_trace()
		vals = zip(*self.phenotypeValues[:])[idMap[phenotypeID]]
		vals = list(vals)
		indices = []
		values = []
		for i in range(0, len(vals)):
			if vals[i] != "NA":
				indices.append(i)
				values.append(float(vals[i]))
		import util
		quantiles = util.calcQuantiles(values[:])
		print "Quantiles:", quantiles
		median = quantiles[1]
		iqr = abs(quantiles[2] - quantiles[0])
		print iqr
		for i in range(0, len(values)):
			if abs(values[i] - median) > iqrThreshold * iqr:
				print "removed", values[i]
				vals[indices[i]] = "NA"
				numNA += 1
		for i in range(0, len(self.accessions)):
			#print vals[i],idMap,phenotypeID,i
			#pdb.set_trace()
			#print self.phenotypeValues
			self.phenotypeValues[i][idMap[phenotypeID]] = vals[i]
		print "NAed", numNA, "values."
		return numNA

	def getPhenotypeName(self, phenotypeIndex, rawName=False, niceName=False):
		indexMap = self._getIndexMapping_()
		if indexMap.has_key(phenotypeIndex):
			if niceName:
				phenName = self.phenotypeNiceNames[indexMap[phenotypeIndex]]
			else:
				phenName = self.phenotypeNames[indexMap[phenotypeIndex]]
				if not rawName:
					phenName = phenName.replace("<i>", "")
					phenName = phenName.replace("</i>", "")
					phenName = phenName.replace("/", "_div_")
					phenName = phenName.replace("*", "_star_")
					phenName = phenName.replace(" ", "_")
			return phenName
		else:
			print "Phenotype with ID", phenotypeIndex, "not found"
			return False

	def getPhenotypeNiceNamesDict(self):
		name_dict = {}
		for p_i in self.phenIds:
			phen_name = (" ".join(self.getPhenotypeName(p_i).split("_")[1:])).strip()
			nice_name = self.getPhenotypeName(p_i, niceName=True)
			print phen_name, ":" , nice_name
			name_dict[phen_name] = nice_name
		return name_dict


	def get_marker_and_phen_vals(self, phen_i, md, m_i):
		"""
		Returns the phenotype values matched together with a SNP for all non-missing values.
		"""


		phen_vals = self.getPhenVals(phen_i, False, False)
		m = md.snps[m_i]
		new_phen_vals = []
		new_accessions = []
		new_marker = []
		for i, a in enumerate(self.accessions):
			if a in md.accessions:
				m_ai = md.accessions.index(a)
				if phen_vals[i] != 'NA' and m[m_ai] != md.missingVal:
					new_phen_vals.append(phen_vals[i])
					new_marker.append(m[m_ai])
					new_accessions.append(a)


		return {"marker":new_marker, "phen_vals":new_phen_vals, "accessions":new_accessions}


	def plot_histogram(self, p_i, title=None , pdfFile=None, pngFile=None, withLabels=False):
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		plt.figure(figsize=(5, 4))
		plt.axes([0.14, 0.13, 0.81, 0.83])
		if withLabels:
			label = loadHistogramLabels()[p_i]
			plt.xlabel(label)
		phenValues = self.getPhenVals(p_i)

		minVal = min(phenValues)
		maxVal = max(phenValues)
		x_range = maxVal - minVal
		histRes = plt.hist(phenValues, bins=len(phenValues) / 4)
		y_max = max(histRes[0])
		plt.axis([minVal - 0.035 * x_range, maxVal + 0.035 * x_range, -0.035 * y_max, 1.1 * y_max])
		num_phen_vals = len(phenValues)
		plt.text(maxVal - 0.7 * x_range, 1.02 * y_max, "Number of values: " + str(num_phen_vals), size="x-small")
		print max(histRes[0])
		plt.ylabel("Frequency")
		if title:
			plt.title(title)
		if pdfFile:
			plt.savefig(pdfFile, format="pdf")
		if pngFile:
			plt.savefig(pngFile, format="png", dpi=300)
		elif not pdfFile:
			plt.show()
		plt.clf()


	def plot_marker_box_plot(self, phen_i, marker=None, marker_accessions=None, md=None, m_i=None,
				 pdf_file=None, png_file=None, title=None, marker_score=None, marker_missing_val='NA'):
		"""
		Plots a box plot for the given binary marker and phenotype. 
		"""
		marker = list(marker)
		import pylab
		if md and m_i:
			r = self.get_marker_and_phen_vals(phen_i, md, m_i)
			marker = r['marker']
			phen_vals = r['phen_vals']
			missing_val = md.missingVal
		elif marker and marker_accessions:
#			new_marker = []
#			new_phen_vals = []
			phen_vals = self.getPhenVals(phen_i, noNAs=False)
			if len(marker_accessions) != len(phen_vals):
				raise Exception
#			for nt,a in zip(marker,marker_accessions):
#				if a in self.accessions:
#					a_i = self.accessions.index(a)
#					phen_val = phen_vals[a_i]
#					if phen_val!='NA' and nt!='NA':
#						new_marker.append(nt)
#						new_phen_vals.append(phen_val)
#			phen_vals = new_phen_vals
#			marker = new_marker		
			missing_val = marker_missing_val

		nt_counts = [(marker.count(nt), nt) for nt in set(marker) if nt != missing_val]
		nt_counts.sort()
		if len(nt_counts) > 2:
			import warnings
			warnings.warn("More than 2 alleles, box-plot might be wrong?")

		allele_phen_vals = []
		for c, nt in nt_counts:
			allele_phen_vals.append([v for i, v in enumerate(phen_vals) if marker[i] == nt])
			pylab.figure()
		pylab.boxplot(allele_phen_vals)
		minVal = min(phen_vals)
		maxVal = max(phen_vals)
		rangeVal = maxVal - minVal
		pylab.axis([0.2, 2.8, minVal - rangeVal * 0.3, maxVal + rangeVal * 0.3])
		pylab.text(0.4, minVal - 0.15 * rangeVal, "# of obs.: ", color='k')
		pylab.text(0.95, minVal - 0.15 * rangeVal, str(len(allele_phen_vals[0])), color='k')
		if len(allele_phen_vals) > 1:
			pylab.text(1.95, minVal - 0.15 * rangeVal, str(len(allele_phen_vals[1])), color='k')
		if marker_score:
			pylab.text(0.9, maxVal + 0.15 * rangeVal, "-log(p-value)/score: " + str(marker_score), color='k')
		if title:
			pylab.title(title)
		elif md and m_i:
			pylab.title(self.getPhenotypeName(phen_i) + ": chromosome " + str(md.chromosome) + ", position " + str(md.positions[m_i]))
		if pdf_file:
			pylab.savefig(pdf_file, format="pdf")
		if png_file:
			pylab.savefig(png_file, format="png")
		if not pdf_file and not png_file:
			pylab.show()



	def plot_accession_map(self, p_i, pdf_file=None, png_file=None, color_by=None, map_type="global",
			       title=None, accessions=None, cmap=None):
		"""
		Plot accessions on a map.
		
		'color_by' is by default set to be the phenotype values.
		"""
		if not accessions:
			accessions = self.accessions
		eid = get_ecotype_id_info_dict()
		lats = []
		lons = []
		for e in accessions:
			r = eid[int(e)]
			try:
				latitude = float(r[2])
				longitude = float(r[3])

			except Exception, err_str:
				print "Latitude and Longitude, not found?:", err_str
				print 'Placing them in the Atlantic.'
				latitude = 40
				longitude = -20
			lats.append(latitude)
			lons.append(longitude)

		from mpl_toolkits.basemap import Basemap
		import numpy as np
		import matplotlib.pyplot as plt
		from pylab import cm
		if map_type == "global2":
			plt.figure(figsize=(14, 12))
			m = Basemap(width=21.e6, height=21.e6, projection='gnom', lat_0=76, lon_0=15)
			m.drawparallels(np.arange(20, 90, 20))
			m.drawmeridians(np.arange(-180, 180, 30))
		elif map_type == 'global':

			plt.figure(figsize=(16, 4))
			plt.axes([0.02, 0.02, 0.96, 0.96])
			m = Basemap(projection='cyl', llcrnrlat=10, urcrnrlat=80,
				    llcrnrlon= -130, urcrnrlon=150, lat_ts=20, resolution='c')
			m.drawparallels(np.arange(20, 90, 20))
			m.drawmeridians(np.arange(-180, 180, 30))
		else:
			raise Exception("map_type is invalid")

		#m.drawmapboundary(fill_color='aqua')
		m.drawcoastlines(zorder=0)
		m.fillcontinents(zorder=1)
		#m.fillcontinents(color='coral',lake_color='aqua')

		xs = []
		ys = []
		for lon, lat in zip(lons, lats):
			x, y = m(*np.meshgrid([lon], [lat]))
			xs.append(x)
			ys.append(y)

		if not color_by:
			color_vals = self.getPhenVals(p_i)
		else:
			color_vals = color_by
		if len(color_vals) != len(accessions):
			raise Exception("accessions and color_by_vals values don't match ! ")
		if not cmap:
			num_colors = len(set(color_vals))
			if num_colors <= 10:
				cmap = cm.get_cmap('jet', num_colors)
			else:
				cmap = cm.get_cmap('jet')
		lws = [0] * len(xs)
		plt.scatter(xs, ys, s=10, linewidths=lws, c=color_vals, cmap=cmap, alpha=0.7, zorder=2)
		plt.colorbar()
		if title:
			plt.title(title)
		if pdf_file:
			plt.savefig(pdf_file, format="pdf")
		if png_file:
			plt.savefig(png_file, format="png")
		if not pdf_file and not png_file:
			plt.show()


	def logTransform(self, pIndex):
		indexMap = self._getIndexMapping_() #Modified 9/26/08
		phenotypeIndex = indexMap[pIndex]

		if not self.isBinary(pIndex) and not self._lessOrEqualZero_(phenotypeIndex):
			import math
			for i in range(0, len(self.accessions)):
				if self.phenotypeValues[i][phenotypeIndex] != 'NA':
					self.phenotypeValues[i][phenotypeIndex] = str(math.log(float(self.phenotypeValues[i][phenotypeIndex])))
		else:
			print "Can't log-transform, since phenotype is binary OR values are out of logarithm range!"
			return False
		return True


	#def standard_transform(self, p_i):
	def log_transform(self, p_i):
		"""
		Adds a constant before applying log-transformation.
		"""
		self.addSDscaledConstant(p_i)
		self.logTransform(p_i)


	def sqrt_transform(self, p_i):
		"""
		Adds a constant before applying log-transformation.
		"""
		self.addSDscaledConstant(p_i)
		pi = self.getPhenIndex(p_i)

		if not self.isBinary(p_i) and not self._lessOrEqualZero_(pi):
			for i in range(len(self.accessions)):
				if self.phenotypeValues[i][pi] != 'NA':
					self.phenotypeValues[i][pi] = str(math.sqrt(float(self.phenotypeValues[i][pi])))
		else:
			print "Can't sqrt - transform, since phenotype is binary OR values are out of sqrt range ! "
			return False
		return True


	def transformToRanks(self, pIndex):
		"""
		Transformes the phenotypic values to ranks.
		"""
		indexMap = self._getIndexMapping_()
		phenotypeIndex = indexMap[pIndex]

		values = []
		for i in range(0, len(self.accessions)):
			if self.phenotypeValues[i][phenotypeIndex] != 'NA':
				values.append((self.phenotypeValues[i][phenotypeIndex], i))

		values.sort()
		ranks = []
		for i in range(0, len(values)):
			vals = values[i]
			ranks.append(vals[1])


		j = 1
		for r_i in ranks:
			self.phenotypeValues[r_i][phenotypeIndex] = float(j)
			j += 1




	def _lessOrEqualZero_(self, phenotypeIndex):
		lessOrEqualZero = False
		minVal = None
		for i in range(0, len(self.accessions)):
			if self.phenotypeValues[i][phenotypeIndex] != 'NA' and float(self.phenotypeValues[i][phenotypeIndex]) <= 0.0:
				lessOrEqualZero = True
				minVal = float(self.phenotypeValues[i][phenotypeIndex])
				break
			if self.phenotypeValues[i][phenotypeIndex] != 'NA' and (minVal == None or float(self.phenotypeValues[i][phenotypeIndex]) < minVal):
				minVal = float(self.phenotypeValues[i][phenotypeIndex])
		print "Minimum value =", minVal
		return lessOrEqualZero


	def isBinary(self, phenotypeIndex):
		indexMap = self._getIndexMapping_() #Modified 9/26/08
		phenotypeIndex = indexMap[phenotypeIndex]

		l = []
		for i in range(0, len(self.accessions)):
			val = self.phenotypeValues[i][phenotypeIndex]
			if val != 'NA':
				if not val in l:
					l.append(val)

		if 1 < len(l) < 3:
			return True
		elif 1 == len(l):
			raise Exception("Only one phenotype value")
		return False

	def countValues(self, phenotypeIndex):
		indexMap = self._getIndexMapping_() #Modified 9/26/08
		phenotypeIndex = indexMap[phenotypeIndex]

		valCount = 0
		for i in range(0, len(self.accessions)):
			val = self.phenotypeValues[i][phenotypeIndex]
			if val != 'NA':
				valCount += 1
		return valCount

	def countPhenotypes(self):
		return len(self.phenotypeNames)



	def orderAccessions(self, accessionMapping=None):
		"""
		Orders the accession alphabetically if no mapping is given.
		"""
		print "Ordering phenotype data accessions."
		newAccessions = ["" for a in self.accessions]
		newPhenotVals = [[] for a in self.accessions]

		if not accessionMapping:
			accessionMapping = []
			l = range(0, len(self.accessions))
			l1 = zip(self.accessions, l)
			l1.sort()
			j = 0
			for (acc, i) in l1:
				accessionMapping.append((i, j))
				j += 1

		for (i, j) in accessionMapping:
			#print j, len(newAccessions)
			newAccessions[j] = self.accessions[i]
			newPhenotVals[j] = self.phenotypeValues[i]
		self.accessions = newAccessions
		self.phenotypeValues = newPhenotVals

		if self.accessionNames:
			newAccessionNames = []
			for acc in self.accessions:
				newAccessionNames.append("")
			for (i, j) in accessionMapping:
				newAccessionNames[j] = self.accessionNames[i]
			self.accessionNames = newAccessionNames


	def addConstant(self, phenotypeID, constant):
		i = self.getPhenIndex(phenotypeID)
		for j in range(0, len(self.accessions)):
			if self.phenotypeValues[j][i] != "NA":
				self.phenotypeValues[j][i] = str(float(self.phenotypeValues[j][i]) + constant)


	def addSDscaledConstant(self, p_i, scale=0.1):
		import math
		addConstant = math.sqrt(self.getVariance(p_i)) * scale
		addConstant = addConstant - self.getMinValue(p_i)
		print "Adding a constant to phenotype", p_i, ":", addConstant
		self.addConstant(p_i, addConstant)

	def negateValues(self, phenotypeID):
		i = self.getPhenIndex(phenotypeID)
		for j in range(0, len(self.accessions)):
			self.phenotypeValues[j][i] = str(-float(self.phenotypeValues[j][i]))


	def removeAccessionsNotInSNPsData(self, snpsd):
		indicesToKeep = []
		for i in range(0, len(self.accessions)):
			acc = self.accessions[i]
			if acc in snpsd.accessions:
				indicesToKeep.append(i)
		self.removeAccessions(indicesToKeep)


	def removeAccessions(self, indicesToKeep):
		"""
		Removes accessions from the data. (indices based)
		"""
		numAccessionsRemoved = len(self.accessions) - len(indicesToKeep)
		print "Removing", numAccessionsRemoved, "accessions in phenotype data, out of", len(self.accessions), "accessions."
		newAccessions = []
		newPhenotVals = []
		for i in indicesToKeep:
			newAccessions.append(self.accessions[i])
			newPhenotVals.append(self.phenotypeValues[i])
		self.accessions = newAccessions
		self.phenotypeValues = newPhenotVals
#		print "len(self.accessions):", len(self.accessions)
#		print "len(self.phenotypeValues):", len(self.phenotypeValues)
		if self.accessionNames:
			newAccessionNames = []
			for i in indicesToKeep:
				newAccessionNames.append(self.accessionNames[i])
			self.accessionNames = newAccessionNames
			print "len(self.accessionNames):", len(self.accessionNames)


	def filterAccessions(self, accessionsToKeep):
		"""
		Filter accessions in the data.
		"""
		disregardedEcotypes = set(self.accessions)
 		disregardedEcotypes = disregardedEcotypes.difference(set(accessionsToKeep))

		indicesToKeep = []
		for i in range(0, len(self.accessions)):
		 	acc = self.accessions[i]
		 	if acc in accessionsToKeep:
		 		indicesToKeep.append(i)
		self.removeAccessions(indicesToKeep)

		return disregardedEcotypes

	def filter_accessions_w_missing_data(self, pid=None):
		if pid:
			i = self.getPhenIndex(pid)
			indices_to_keep = [ai for ai, a in enumerate(self.accessions) \
					if 'NA' != self.phenotypeValues[ai][i]]
			self.removeAccessions(indices_to_keep)
		else:
			indices_to_keep = [ai for ai, a in enumerate(self.accessions) \
					if not 'NA' in self.phenotypeValues[ai]]
			self.removeAccessions(indices_to_keep)




	def getNonNAEcotypes(self, phenotypeID):
		i = self.getPhenIndex(phenotypeID)
		ecotypes = []
		for j in range(0, len(self.accessions)):
			if self.phenotypeValues[j][i] != 'NA':
				ecotypes.append(self.accessions[j])
		return ecotypes


	def filter_na_ecotypes(self):
		"""
		Removes accessions with missing data in all traits.
		"""
		e_set = set()
		for pid in self.phenIds:
			e_set = e_set.union(set(self.getNonNAEcotypes(pid)))
		self.filterAccessions(list(e_set))


	def insert_into_DB(self, pids=None, phenotype_scoring='',
			   method_description='', growth_condition='', biology_category_id='',
			   citations='', data_description='', transformation_description=None,
			   method_id=None, data_type=None, comment=''):
		"""
		Inserts phenotypes into DB, but assumes that values are averages.
		(Assumes no transformation)
		"""

		if not pids:
			pids = self.phenIds
		import dbutils
		conn = dbutils.connect_to_default_insert("stock_250k")
		cursor = conn.cursor()
		for fpid, pid in enumerate(pids):
			try:
				phen_name = self.getPhenotypeName(pid)
				cur_method_id = self.get_db_pid(pid, conn=conn)
				phen_vals = self.getPhenVals(pid, noNAs=True)
				ecotypes = self.getNonNAEcotypes(pid)
				no_of_accessions = len(ecotypes)
				print cur_method_id
				if cur_method_id:
					print phen_name, 'is already in DB.'
					print 'Updating values only.'

					sql_statement = "UPDATE stock_250k.phenotype_method SET no_of_accessions=%d, biology_category_id=%d WHERE short_name='%s'" % (no_of_accessions, biology_category_id, phen_name)
					print sql_statement
					cursor.execute(sql_statement)
				else:
					if not data_type:
						if self.isBinary(pid):
							data_type = 'binary'
						else:
							data_type = 'quantitative'
					print "Inserting phenotype %s into DB." % phen_name
					if method_id:
						sql_statement = "INSERT INTO stock_250k.phenotype_method  (id, short_name, only_first_96, no_of_accessions," + \
						     " biology_category_id,phenotype_scoring, method_description, growth_condition, data_description, comment," + \
						     " data_type, transformation_description) VALUES (" + str(method_id) + ", '" + phen_name + "', false, " + \
						     str(no_of_accessions) + ", " + str(biology_category_id) + ", '" + phenotype_scoring + "', '" + method_description + \
						     "', '" + growth_condition + "', '" + data_description + "', '" + comment + "', '" + data_type + "', '" + \
						     str(transformation_description) + "')"
					else:
						sql_statement = "INSERT INTO stock_250k.phenotype_method  (short_name, only_first_96, no_of_accessions," + \
						     " biology_category_id,phenotype_scoring, method_description, growth_condition, data_description," + \
						     " comment, data_type, transformation_description) VALUES ('" + phen_name + "', false, " + \
						     str(no_of_accessions) + ", " + str(biology_category_id) + ", '" + phenotype_scoring + "', '" + method_description + \
						     "', '" + growth_condition + "', '" + data_description + "', '" + comment + "', '" + data_type + "', '" + \
						     str(transformation_description) + "')"

					print sql_statement
					numRows = int(cursor.execute(sql_statement))
					row = cursor.fetchone()
					if row:
						print row
					cur_method_id = self.get_db_pid(pid, conn=conn)

				print "Inserting values for phenotype method_id = %d" % cur_method_id
				for e_i, val in zip(ecotypes, phen_vals):
					sql_statement = "INSERT INTO stock_250k.phenotype_avg  (ecotype_id, value, ready_for_publication, " + \
							"method_id, transformed_value) VALUES ( " + str(e_i) + ", " + str(val) + ", 0, " + \
							str(cur_method_id) + ", " + str(val) + " )"
					try:
						print sql_statement
						numRows = int(cursor.execute(sql_statement))
						row = cursor.fetchone()
						if row:
							print row
					except Exception, err_str:
						print 'Failed at inserting data:', err_str
						print 'Trying to update.'
						sql_statement = "UPDATE stock_250k.phenotype_avg SET value=%d, transformed_value=%d \
								 WHERE ecotype_id=%d and method_id=%d" % (val, val, int(e_i), cur_method_id)
						print sql_statement
						numRows = int(cursor.execute(sql_statement))
						row = cursor.fetchone()
						if row:
							print row

				print "Committing"
				conn.commit()
			except Exception, err_str:
				print 'Failed at inserting phentoype nr.', fpid + 1, ':', err_str

		cursor.close ()
		conn.close ()




	def writeToFile(self, outputFile, phenotypes=None, delimiter=',', with_pid=False,
			with_accession_names=False):
		print "Writing out phenotype file:", outputFile
		outStr = "ecotype_id"
		if with_accession_names:
			if not self.accessionNames:
				ei_dict = get_ecotype_id_info_dict()
				self.accessionNames = []
				for ei in self.accessions:
					self.accessionNames.append(ei_dict[int(ei)][0])
			outStr += delimiter + "accession_name"
		if phenotypes:
			for i in phenotypes:
				if with_pid:
					name = str(self.phenIds[i]) + '_' + self.phenotypeNames[i]
				else:
					name = self.phenotypeNames[i]
				outStr += delimiter + name
			outStr += '\n'
			for i in range(0, len(self.accessions)):
				outStr += str(self.accessions[i])
				for j in phenotypes:
					outStr += delimiter + str(self.phenotypeValues[i][j])
				outStr += "\n"
		else:
			for i, name in enumerate(self.phenotypeNames):
				if with_pid:
					outStr += delimiter + str(self.phenIds[i]) + '_' + name
				else:
					outStr += delimiter + name
			outStr += '\n'
			for i in range(0, len(self.accessions)):
				outStr += str(self.accessions[i])
				if with_accession_names and self.accessionNames:
					outStr += delimiter + str(self.accessionNames[i])
				for j in range(0, len(self.phenotypeNames)):
					outStr += delimiter + str(self.phenotypeValues[i][j])
				outStr += "\n"
				#print outStr

		f = open(outputFile, "w")
		f.write(outStr)
		f.close()


def readPhenotypeFile(filename=None, file_object=None, delimiter=None, missing_value_decoder=None, accessionDecoder=None, with_db_ids=True):
	"""
	Reads a phenotype file and returns a phenotype object.
	"""
	if not missing_value_decoder:
		missing_value_decoder = {'NA':'NA'}
	if filename:
		f = open(filename, "r")
	elif file_object:
		f = file_object
	else:
		raise Exception('Error')

	lines = f.readlines()
	f.close()
	shift = 1
	if not delimiter: #Try different delimiters..
		for delimiter in ['\t', ',']:
			line = (lines[0].rstrip()).split(delimiter)
			if len(line) > 1:
				break
		if len(line) == 1:
			raise Exception("Either delimiter is wrong, or file format is bad.")

	if line[1] in ["accession_name", 'Ecotype'] :
		shift = 2
	phenotypeNames = line[shift:]
	accessions = []
	phenotypeValues = []
	for i in range(1, len(lines)):
		line = (lines[i].rstrip()).split(delimiter)
		if len(line) <= shift:
			print "Phenotypes are completely missing for accession: %s, hence skipping this one." % str(line)
			continue
		if accessionDecoder:
			accessions.append(accessionDecoder[line[0]])
		else:
			accessions.append(line[0])

		phen_vals = []
		for phen_val in line[shift:]:
			if phen_val in missing_value_decoder:
				phen_vals.append(missing_value_decoder[phen_val])
			else:
				phen_vals.append(float(phen_val))

		phenotypeValues.append(phen_vals)

	return PhenotypeData(accessions, phenotypeNames, phenotypeValues, with_db_ids=with_db_ids)


# A list of publishable phenotypes and their group.
_publishablePhenotypes_ = []



def get_all_phenotypes_from_db(file_name=None):
	"""
	Fetches all raw phenotype data from the DB, including replicates.
	
	Checks whether averages are consistent with replicates.
	"""
	import dbutils
	conn = dbutils.connect_to_default_lookup()
	cursor = conn.cursor ()

	#Get the phenotypic values
	sql_statement = "\
	SELECT DISTINCT pa.method_id, pa.ecotype_id, pa.value, pm.short_name\
	FROM stock_250k.phenotype_avg pa, stock_250k.phenotype_method pm \
	WHERE pa.method_id=pm.id ORDER BY pa.method_id, pa.ecotype_id"
	print sql_statement
	numRows = int(cursor.execute(sql_statement))

	phen_dict = {}
	num_values = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		pid = int(row[0])
		phen_name = row[3]
		if pid not in phen_dict:
			phen_dict[pid] = {'name':phen_name, 'ecotypes':[], 'values':[]}

		if row[2] != None:
			phen_dict[pid]['values'].append(float(row[2]))
			phen_dict[pid]['ecotypes'].append(int(row[1]))
			num_values += 1

	print 'Values for %d phenotypes were fetched from DB.' % len(phen_dict)
	print "%d values in total." % num_values

	avg_phend = phenotype_data(phen_dict)

	#Now fetch the replicates
#	sql_statement = "\
#	SELECT DISTINCT p.method_id, p.ecotype_id, p.value, p.replicate, pm.short_name\
#	FROM stock_250k.phenotype p, stock_250k.phenotype_method pm \
#	WHERE p.method_id=pm.id ORDER BY p.method_id, p.ecotype_id"
	sql_statement = "\
	SELECT DISTINCT p.method_id, ei.tg_ecotypeid, p.value, p.replicate, pm.short_name\
	FROM stock_250k.phenotype p, stock_250k.phenotype_method pm, stock.ecotypeid2tg_ecotypeid ei \
	WHERE ei.ecotypeid=p.ecotype_id and p.method_id=pm.id ORDER BY p.method_id, ei.tg_ecotypeid"
	print sql_statement
	numRows = int(cursor.execute(sql_statement))
	repl_phen_dict = {}
	num_values = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		pid = int(row[0])
		phen_name = row[4]
		if pid not in repl_phen_dict:
			repl_phen_dict[pid] = {'name':phen_name, 'ecotypes':[], 'values':[], 'replicate':[]}

		if row[2] != None:
			repl_phen_dict[pid]['values'].append(float(row[2]))
			repl_phen_dict[pid]['ecotypes'].append(int(row[1]))
			repl_phen_dict[pid]['replicate'].append(int(row[3]))
			num_values += 1

	print 'Values for %d phenotypes were fetched from DB.' % len(repl_phen_dict)
	print "%d values in total." % num_values

	cursor.close ()
	conn.close ()
	cm72_dict = get_250K_accession_to_ecotype_dict(dict_key='ecotype_id')
	#print cm72_dict
	repl_phend = phenotype_data(repl_phen_dict)
	no_problem_pids = []
	for pid in phen_dict:
		if pid in repl_phen_dict:
			d1 = avg_phend.get_avg_values(pid)
			d2 = repl_phend.get_avg_values(pid)
			ets1 = d1['ecotypes']
			ets2 = d2['ecotypes']
			vls1 = d1['values']
			vls2 = d2['values']
			l1 = zip(ets1, vls1)
			l2 = zip(ets2, vls2)
			l1.sort()
			l2.sort()
			et_set1 = set(ets1)
			et_set2 = set(ets2)
			phen_name = phen_dict[pid]['name']
			print phen_name, ', pid:', pid, ', et_set1==et_set2 :', et_set1 == et_set2
			common_ets = et_set1.intersection(et_set2)
			num_diffs = 0
			for et in common_ets:
				t = (sp.var(d1['values']) + sp.var(d2['values'])) / 100.0
				i1 = d1['ecotypes'].index(et)
				i2 = d2['ecotypes'].index(et)
				if (d1['values'][i1] - d2['values'][i2]) > t:
					print d1['values'][i1], d2['values'][i2]
					num_diffs += 1
			if not et_set1 == et_set2:
				print len(et_set1), len(et_set2)
				s = et_set1.symmetric_difference(et_set2)
				#print 'et_set1.difference(et_set2):', et_set1.difference(et_set2)
				#print 'et_set2.difference(et_set1):', et_set2.difference(et_set1)
				if len(s) < 8 and num_diffs == 0:
					print 'No major problems'
					no_problem_pids.append(pid)
					continue
#				for ei in s:
#					if str(ei) in cm72_dict:
#						print cm72_dict[str(ei)]
#					else:
#						print ei, 'is not in 250K data!'
			elif num_diffs == 0:
				print 'No problems?'
				no_problem_pids.append(pid)
	print no_problem_pids

	merged_phen_dict = {}
	for pid in phen_dict:
		if pid in no_problem_pids:
			d = {'name':phen_dict[pid]['name'], 'values': repl_phen_dict[pid]['values'],
				'ecotypes': repl_phen_dict[pid]['ecotypes']}

			et_set1 = set(phen_dict[pid]['ecotypes'])
			et_set2 = set(repl_phen_dict[pid]['ecotypes'])
			for et in et_set1 - et_set2:
				vi = phen_dict[pid]['ecotypes'].index(et)
				val = phen_dict[pid]['values'][vi]
				d['ecotypes'].append(et)
				d['values'].append(val)
			merged_phen_dict[pid] = d
		else:
			merged_phen_dict[pid] = phen_dict[pid]


	phend = phenotype_data(phen_dict=merged_phen_dict)
	if file_name:
		phend.write_to_file(file_name)
	return phend



def getPhenotypes(onlyBinary=False, onlyQuantitative=False, onlyCategorical=False, onlyReplicates=False,
		includeSD=False, rawPhenotypes=False, onlyPublishable=False):
	import dbutils
	conn = dbutils.connect_to_default_lookup()
	cursor = conn.cursor ()

	#Retrieve the ecotypes
	print "Fetching ecotypes"
	if onlyPublishable:
		sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename \
		from stock_250k.phenotype_avg pa, stock.ecotypeid2tg_ecotypeid ei \
		where ei.ecotypeid=pa.ecotype_id and pa.ready_for_publication=1 order by ei.tg_ecotypeid"
		numRows = int(cursor.execute(sql_statement))
	else:
		sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename \
		from stock_250k.phenotype_avg pa, stock.ecotypeid2tg_ecotypeid ei \
		where ei.ecotypeid=pa.ecotype_id order by ei.tg_ecotypeid"
		numRows = int(cursor.execute(sql_statement))


	ecotypes = []
	accessions = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		ecotypes.append(str(int(row[0])))
		accessions.append(row[1])
	print len(ecotypes), "ecotypes (accessions) were found."

	#Get the phenotypic values
	phenotypeValues = []
	phenotypeNames = []
	valueColumn = "pa.transformed_value"
	if rawPhenotypes:
		print "Fetching raw phenotypic values"
		valueColumn = "pa.value"
	else:
		print "Fetching phenotypic values"

	if onlyPublishable:
		sql_statement = "select distinct pa.method_id, ei.tg_ecotypeid, pa.value, \
				pa.transformed_value, pm.short_name from stock_250k.phenotype_avg pa, \
				stock.ecotypeid2tg_ecotypeid ei, stock_250k.phenotype_method pm \
				where ei.ecotypeid=pa.ecotype_id and pa.method_id=pm.id and pa.ready_for_publication=1 \
				order by pa.method_id, ei.tg_ecotypeid"
		numRows = int(cursor.execute(sql_statement))
	else:
		sql_statement = "select distinct pa.method_id, ei.tg_ecotypeid, pa.value, \
				pa.transformed_value, pm.short_name from stock_250k.phenotype_avg pa, \
				stock.ecotypeid2tg_ecotypeid ei, stock_250k.phenotype_method pm \
				where ei.ecotypeid=pa.ecotype_id and pa.method_id=pm.id order by pa.method_id, \
				ei.tg_ecotypeid"
		print sql_statement
		numRows = int(cursor.execute(sql_statement))

	pvalues = [[] for j in range(0, len(ecotypes))]
	trans_pvalues = [[] for j in range(0, len(ecotypes))]
	row = cursor.fetchone()
	currentMethod = int(row[0])
	phenName = str(int(row[0])) + "_" + row[4]
	#print phenName
	phenotypeNames.append(phenName)
	i = ecotypes.index(str(int(row[1])))
	if row[2] != None and row[2] != 'NA':
		pvalues[i] += [float(row[2])]
	if row[3] != None and row[3] != 'NA':
		trans_pvalues[i] += [float(row[3])]

	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		nextMethod = int(row[0])
		if currentMethod != nextMethod:  #If we've encountered the next phenotype_method
			phenName = str(int(row[0])) + "_" + row[4]
			#print phenName
			phenotypeNames.append(phenName)
			for j in range(0, len(pvalues)):
				if pvalues[j] == [] or pvalues[j] == None:
					pvalues[j] = "NA"
				else:
					pvalues[j] = sum(pvalues[j]) / float(len(pvalues[j]))
				if trans_pvalues[j] == [] or trans_pvalues[j] == None:
					trans_pvalues[j] = "NA"
				else:
					trans_pvalues[j] = sum(trans_pvalues[j]) / float(len(trans_pvalues[j]))

			if rawPhenotypes:
				phenotypeValues.append(pvalues)
			else:
				new_phen_vals = []
				for phen_val, trans_phen_val in zip(pvalues, trans_pvalues):
					if trans_phen_val == 'NA':
						new_phen_vals.append(phen_val)
					else:
						new_phen_vals.append(trans_phen_val)
				phenotypeValues.append(new_phen_vals)

			pvalues = [[] for j in range(0, len(ecotypes))]
			trans_pvalues = [[] for j in range(0, len(ecotypes))]
			currentMethod = nextMethod

		i = ecotypes.index(str(int(row[1])))
		if row[2] != None:
			pvalues[i] += [float(row[2])]
		if row[3] != None:
			trans_pvalues[i] += [float(row[3])]

	for j in range(0, len(pvalues)):
		if pvalues[j] == [] or pvalues[j] == None:
			pvalues[j] = "NA"
		else:
			pvalues[j] = sum(pvalues[j]) / float(len(pvalues[j]))
		if trans_pvalues[j] == [] or trans_pvalues[j] == None:
			trans_pvalues[j] = "NA"
		else:
			trans_pvalues[j] = sum(trans_pvalues[j]) / float(len(trans_pvalues[j]))
	if rawPhenotypes:
		phenotypeValues.append(pvalues)
	else:
		new_phen_vals = []
		for phen_val, trans_phen_val in zip(pvalues, trans_pvalues):
			if trans_phen_val == 'NA':
				new_phen_vals.append(phen_val)
			else:
				new_phen_vals.append(trans_phen_val)
		phenotypeValues.append(new_phen_vals)
	print len(phenotypeValues), "phenotypes were found."
	print len(phenotypeNames)

	cursor.close ()
	conn.close ()

	phenotypeValues = zip(*phenotypeValues)
	phenotypeValues = map(list, phenotypeValues)

	phenDat = PhenotypeData(ecotypes, phenotypeNames, phenotypeValues, accessionNames=accessions)
	return phenDat


def _getFirst96Ecotypes_(host="papaya.usc.edu", user="bvilhjal", passwd="*rri_bjarni@usc"):
	import MySQLdb
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	print "Fetching data"
	ecotypes = []
	numRows = int(cursor.execute("select distinct e2te.tg_ecotypeid, c2010.accession_id, e2te.nativename, e2te.stockparent from at.complete_2010_strains_in_stock c2010, stock.ecotypeid2tg_ecotypeid e2te where e2te.ecotypeid=c2010.ecotypeid and e2te.stockparent = c2010.stockparent and c2010.accession_id<97 order by c2010.accession_id"))
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		ecotypes.append(int(row[0]))
	cursor.close ()
	conn.close ()
	return ecotypes

def _getFirst192Ecotypes_(host="papaya.usc.edu", user="bvilhjal", passwd="*rri_bjarni@usc"):
	"""
	Result of a union of all the "192" accessions for the phenotypes to be used for the GWA 2009 paper.
	"""
	return map(int, ['100000', '5837', '6008', '6009', '6016', '6024', '6039', '6040', '6042', '6043', '6046', '6064', '6074', '6088', '6243', '6709', '6830', '6897', '6898', '6899', '6900', '6901', '6903', '6904', '6905', '6906', '6907', '6908', '6909', '6910', '6911', '6913', '6914', '6915', '6916', '6917', '6918', '6919', '6920', '6921', '6922', '6923', '6924', '6926', '6927', '6928', '6929', '6930', '6931', '6932', '6933', '6936', '6937', '6938', '6939', '6940', '6942', '6943', '6944', '6945', '6946', '6951', '6956', '6957', '6958', '6959', '6960', '6961', '6962', '6963', '6964', '6965', '6966', '6967', '6968', '6969', '6970', '6971', '6972', '6973', '6974', '6975', '6976', '6977', '6978', '6979', '6980', '6981', '6982', '6983', '6984', '6985', '6988', '7000', '7014', '7033', '7062', '7064', '7081', '7094', '7123', '7147', '7163', '7231', '7255', '7275', '7282', '7296', '7306', '7323', '7327', '7329', '7346', '7418', '7424', '7438', '7460', '7461', '7477', '7514', '7515', '7516', '7517', '7518', '7519', '7520', '7521', '7522', '7523', '7524', '7525', '7526', '8213', '8214', '8215', '8222', '8230', '8231', '8233', '8234', '8235', '8236', '8237', '8239', '8240', '8241', '8242', '8243', '8244', '8245', '8246', '8247', '8249', '8254', '8256', '8257', '8258', '8259', '8264', '8265', '8266', '8270', '8271', '8274', '8275', '8283', '8284', '8285', '8290', '8296', '8297', '8300', '8304', '8306', '8307', '8310', '8311', '8312', '8313', '8314', '8323', '8325', '8326', '8329', '8334', '8335', '8337', '8343', '8351', '8353', '8354', '8356', '8357', '8365', '8366', '8369', '8374', '8376', '8378', '8387', '8388', '8389', '8395', '8411', '8412', '8419', '8420', '8422', '8423', '8424', '8426', '8428', '8430', '9057', '9058'])

 	#accessions = set(['100000', '5837', '6008', '6009', '6016', '6024', '6039', '6040', '6042', '6043', '6046', '6064', '6074', '6088', '6243', '6709', '6830', '6897', '6898', '6899', '6900', '6901', '6903', '6904', '6905', '6906', '6907', '6908', '6909', '6910', '6911', '6913', '6914', '6915', '6916', '6917', '6918', '6919', '6920', '6921', '6922', '6923', '6924', '6926', '6927', '6928', '6929', '6930', '6931', '6932', '6933', '6936', '6937', '6938', '6939', '6940', '6942', '6943', '6944', '6945', '6946', '6951', '6956', '6957', '6958', '6959', '6960', '6961', '6962', '6963', '6964', '6965', '6966', '6967', '6968', '6969', '6970', '6971', '6972', '6973', '6974', '6975', '6976', '6977', '6978', '6979', '6980', '6981', '6982', '6983', '6984', '6985', '6988', '7000', '7014', '7033', '7062', '7064', '7081', '7094', '7123', '7147', '7163', '7231', '7255', '7275', '7282', '7296', '7306', '7323', '7327', '7329', '7346', '7418', '7424', '7438', '7460', '7461', '7477', '7514', '7515', '7516', '7517', '7518', '7519', '7520', '7521', '7522', '7523', '7524', '7525', '7526', '8213', '8214', '8215', '8222', '8230', '8231', '8233', '8234', '8235', '8236', '8237', '8239', '8240', '8241', '8242', '8243', '8244', '8245', '8246', '8247', '8249', '8254', '8256', '8257', '8258', '8259', '8264', '8265', '8266', '8270', '8271', '8274', '8275', '8283', '8284', '8285', '8290', '8296', '8297', '8300', '8304', '8306', '8307', '8310', '8311', '8312', '8313', '8314', '8323', '8325', '8326', '8329', '8334', '8335', '8337', '8343', '8351', '8353', '8354', '8356', '8357', '8365', '8366', '8369', '8374', '8376', '8378', '8387', '8388', '8389', '8395', '8411', '8412', '8419', '8420', '8422', '8423', '8424', '8426', '8428', '8430', '9057', '9058'])
	"""
	import MySQLdb
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	print "Fetching data"
	ecotypes = []
	numRows = int(cursor.execute("select distinct e2te.tg_ecotypeid, c2010.accession_id, e2te.nativename, e2te.stockparent from at.complete_2010_strains_in_stock c2010, stock.ecotypeid2tg_ecotypeid e2te, at.ecotype_192_vs_accession_192 ea192 where e2te.ecotypeid=c2010.ecotypeid and e2te.stockparent = c2010.stockparent and ea192.accession_id =c2010.accession_id order by c2010.accession_id"))
	#numRows = int(cursor.execute("select distinct e2te.tg_ecotypeid, ai2ei.accession_id, e2te.nativename, e2te.stockparent from at.accession2tg_ecotypeid ai2ei, stock.ecotypeid2tg_ecotypeid e2te where e2te.tg_ecotypeid = ai2ei.ecotype_id order by ai2ei.accession_id"))														   
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		ecotypes.append(row[0])
	cursor.close ()
	conn.close ()
	accessions = set(ecotypes)
	"""

#	#Get "192" ecotypes for which we have phenotype data...
#	phenotypes_192 = range(1,5)+range(8,80)+range(81,88)+range(89,161)+[168,169,171,172,174,175]+range(177,187)+range(191,195)
#	phenotypeFile = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phenotypes_all_raw_081009.tsv"
#	phed = readPhenotypeFile(phenotypeFile, delimiter = '\t')
# 	for p_i in phenotypes_192:
# 		accessions = accessions.union(set(phed.getAccessionsWithValues(p_i)))
#	print len(accessions)
#
#
#	#find 250K data overlap
#	import dataParsers
#	snps_accessions = dataParsers.parseCSVDataAccessions("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_data_t43_081009.csv")
#	print snps_accessions
#	accessions = accessions.intersection(set(snps_accessions))
#	print len(accessions)
#	e_dict = _getEcotypeIdToStockParentDict_()
#	f = open("/Users/bjarnivilhjalmsson/tmp/192_accessions_aug09.csv","w")
#	f.write("ecotype_id, accession_name, stock_parent\n")
#	for e_id in accessions:
#		#try:
#			(acc,sp) = e_dict[int(e_id)]
#			acc = unicode(acc,"latin-1")
#			#print acc
#			f.write(e_id+", "+acc+", "+sp+"\n")
##		except Exception, err_str:
##			print err_str,e_id
##			f.write(e_id+", NA, NA\n")
#		
#	f.close()
#	snpsds = dataParsers.parseCSVData("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_data_t43_081009.csv")
#	import snpsdata
#	snps_data = snpsdata.SNPsDataSet(snpsds,[1,2,3,4,5])
#	snps_data.filter_accessions(accessions)
#	print len(snps_data.accessions)
#	snps_data.writeToFile("/tmp/250K_t43_192.csv")
#	snpsds = dataParsers.parseCSVData("/tmp/250K_t43_192.csv")
#
#
#		
#	
#
#
#def _getEcotypeIdToStockParentDict_(host="papaya.usc.edu", user="bvilhjal", passwd="*rri_bjarni@usc", defaultValue=None):
#	import MySQLdb
#	import sys
#	print "Connecting to db, host="+host
#	if not user:
#		sys.stdout.write("Username: ")
#		user = sys.stdin.readline().rstrip()
#	if not passwd:
#		import getpass
#		passwd = getpass.getpass()
#	try:
#		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
#	except MySQLdb.Error, e:
#		print "Error %d: %s" % (e.args[0], e.args[1])
#		sys.exit (1)
#	cursor = conn.cursor ()
#	#Retrieve the filenames
#	print "Fetching data"
#	ecotypeDict = {}
#	numRows = int(cursor.execute("select distinct ei.tg_ecotypeid, ei.nativename, ei.stockparent from stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ai.paternal_ecotype_id=ei.tg_ecotypeid and ai.maternal_ecotype_id=ei.tg_ecotypeid "))
#	i = 0
#	while(1):
#		row = cursor.fetchone()
#		if not row:
#			break;
#		ecotypeDict[int(row[0])] = (row[1],row[2])
#		sp = "NA"
#		if row[2]:
#			sp = row[2]
#		#print str(int(row[0]))+","+str(row[1])+","+sp
#	cursor.close ()
#	conn.close ()
#	return ecotypeDict
#	
#	

def get_ecotype_id_info_dict(defaultValue=None):
	"""
	Returns a dict containing tuples of (nativename, stockparent, latitude, longitude, country)
	"""
	import dbutils
	conn = dbutils.connect_to_default_lookup(db='stock')
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"
	ecotypeDict = {}

	sql_statment = """
select distinct ei.tg_ecotypeid, ei.nativename, ei.stockparent, e.latitude, e.longitude, c.abbr
from stock.ecotype e, stock.ecotypeid2tg_ecotypeid ei, stock.site s, stock.address a, stock.country c
where e.id = ei.tg_ecotypeid and e.siteid = s.id and s.addressid = a.id and a.countryid = c.id
"""
#	sql_statment= """
#select distinct ei.tg_ecotypeid, ei.nativename, ei.stockparent, e.latitude, e.longitude, c.abbr
#from stock.ecotype e, stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai, stock.site s, stock.address a, stock.country c
#where ai.paternal_ecotype_id=ei.tg_ecotypeid and e.id=ei.tg_ecotypeid and ai.maternal_ecotype_id=ei.tg_ecotypeid and e.siteid=s.id and s.addressid=a.id and a.countryid=c.id
#"""
	numRows = int(cursor.execute(sql_statment))
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		latitude = None
		longitude = None
		country = None
		if row[3]:
			latitude = float(row[3])
		if row[4]:
			longitude = float(row[4])
		if row[5]:
			country = row[5]
		ecotypeDict[int(row[0])] = (row[1], row[2], latitude, longitude, country)
		sp = "NA"
		if row[2]:
			sp = row[2]
		#print str(int(row[0]))+","+str(row[1])+","+sp
	cursor.close()
	conn.close()
	return ecotypeDict



#def _get_stock_parent_info_dict_(host="papaya.usc.edu", user="bvilhjal", passwd="*rri_bjarni@usc", defaultValue=None):
#	import MySQLdb
#	import sys
#	print "Connecting to db, host="+host
#	if not user:
#		import sys
#		sys.stdout.write("Username: ")
#		user = sys.stdin.readline().rstrip()
#	if not passwd:
#		import getpass
#		passwd = getpass.getpass()
#	try:
#		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
#	except MySQLdb.Error, e:
#		print "Error %d: %s" % (e.args[0], e.args[1])
#		sys.exit (1)
#	cursor = conn.cursor ()
#	#Retrieve the filenames
#	print "Fetching data"
#	ecotypeDict = {}
#
#	sql_statment= """
#select distinct ei.tg_ecotypeid, ei.nativename, ei.stockparent, e.latitude, e.longitude, c.abbr
#from stock.ecotype e, stock.ecotypeid2tg_ecotypeid ei, stock.site s, stock.address a, stock.country c
#where e.id=ei.tg_ecotypeid and e.siteid=s.id and s.addressid=a.id and a.countryid=c.id
#"""
#
#	numRows = int(cursor.execute(sql_statment))
#	i = 0
#	while(1):
#		row = cursor.fetchone()
#		if not row:
#			break;
#		latitude = None
#		longitude = None
#		country = None
#		if row[3]:
#			latitude = float(row[3])
#		if row[4]:
#			longitude = float(row[4])
#		if row[5]:
#			country = row[5]
#		if row[2]:
#			ecotypeDict[str(row[2])] = (row[0],row[1],latitude,longitude,country)
#		#print row #str(int(row[0]))+","+str(row[1])+","+sp
#	cursor.close()
#	conn.close()
#	return ecotypeDict
#	
#
#
#
#
#def _getEcotype2TgEcotypeDict_(host="papaya.usc.edu", user="bvilhjal", passwd="*rri_bjarni@usc", defaultValue=None):
#	import MySQLdb
#	#print "Connecting to db, host="+host
#	if not user:
#		import sys
#		sys.stdout.write("Username: ")
#		user = sys.stdin.readline().rstrip()
#	if not passwd:
#		import getpass
#		passwd = getpass.getpass()
#	try:
#		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "stock")
#	except MySQLdb.Error, e:
#		print "Error %d: %s" % (e.args[0], e.args[1])
#		sys.exit (1)
#	cursor = conn.cursor ()
#	#Retrieve the filenames
#	#print "Fetching data"
#	eDict = {}
#	numRows = int(cursor.execute("select distinct ecotypeid, tg_ecotypeid from stock.ecotypeid2tg_ecotypeid "))
#	while(1):
#		row = cursor.fetchone()
#		if row:
#			eDict[int(row[0])] = int(row[1])
#		else:
#			break
#	cursor.close ()
#	conn.close ()
#	print "Ecotype to tg_eoctype dictionary was retrieved from the DB."
#	return eDict
#
#
#

def get_250K_accession_to_ecotype_dict(call_method=72, dict_key='nativename'):
	"""
	"""
	delim = '\t'
	fn = env['data_dir'] + 'call_method_72_ecotype_info.tsv'
	f = open(fn)
	header = f.next().split(delim)
	i = header.index(dict_key)
	info_table = []
	ret_dict = {}
	for line in f:
		l = line.split(delim)
		info_table.append(l)
		ret_dict[l[i].strip().lower()] = l
	return ret_dict



def get_accession_to_ecotype_id_dict(accessions, stockParents=None, defaultValue=None, only_250K_accessions=False, lower_cases=True):
	import warnings
	warnings.warn("This function is possibly outdated, please update SQL statement before use!")
	import dbutils
	conn = dbutils.connect_to_default_lookup(db='stock')
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data in DB."
	accDict = {}
	not_found_count = 0
	for i in range(0, len(accessions)):
		acc = accessions[i]
		#acc = unicode(acc,"latin-1")
		#print acc
		if stockParents != None and stockParents[i] != "NA":
			sp = stockParents[i]
			if only_250K_accessions:
				sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename from stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ai.paternal_ecotype_id=ei.tg_ecotypeid and ai.maternal_ecotype_id=ei.tg_ecotypeid and ei.nativename like '" + str(acc.lower()) + "' and ei.stockparent like '" + str(sp) + "'"
			else:
				sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename from stock.ecotypeid2tg_ecotypeid ei where ei.nativename like '" + str(acc.lower()) + "' and ei.stockparent like '" + str(sp) + "'"
		else:
			if only_250K_accessions:
				sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename from stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ai.paternal_ecotype_id=ei.tg_ecotypeid and ai.maternal_ecotype_id=ei.tg_ecotypeid and ei.nativename like '" + str(acc.lower()) + "'"
			else:
				sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename from stock.ecotypeid2tg_ecotypeid ei where ei.nativename like '" + str(acc.lower()) + "'"
		#sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename from stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ei.nativename like '"+str(acc)+"'"
		#print sql_statement
		numRows = int(cursor.execute(sql_statement))
		i = 0
		accession_found = False
		while(1):
			row = cursor.fetchone()
			#print row
			if not row:
				if i == 0:
					not_found_count += 1
					#print "Accession", acc, "wasn't found."
				break
			else:
				if i == 0:
					ecotype = int(row[0])
				else:
					ecotype = min(int(row[0]), ecotype)
				i += 1
		if i:
			if not accDict.has_key(acc.lower()):
				accDict[acc.lower()] = ecotype
				if not lower_cases:
					accDict[acc.upper()] = ecotype
					accDict[acc] = ecotype
			elif  accDict[acc.lower()] > ecotype:
				accDict[acc.lower()] = ecotype
				if not lower_cases:
					accDict[acc.upper()] = ecotype
					accDict[acc] = ecotype
				#print acc.lower(),":", row[1], int(row[0])
	cursor.close ()
	conn.close ()
	if not_found_count:
		print not_found_count, "accessions weren't found."
	else:
		print "All accessions were found in DB."
	return accDict



def _simpleInsertPhenotypesIntoDb_(host="papaya.usc.edu", user="bvilhjal", passwd="bamboo123"):
	ion_number = 19
	method_description = "Ratio: Delta_wet/Delta_dry. (Thomas Juenger)"

	p_i = 3  #FIXME
	m_i = 262  #FIXME
	data_description = ""
	raw_phen_dir = "/Network/Data/250k/tmp-bvilhjal/thomas_juenger/"
	filename = raw_phen_dir + "Delta_Magnus_Set_DEC08.csv"
	f = open(filename, "r")
	lines = f.readlines()
	#lines = lines[0].split("\r")
	print lines
	phenValues = []
	accessions = []
	stockParents = []
	line = lines[0].strip()
	values = line.split(",")
	phenName = "Delta_plasticity"
	for i in range(1, len(lines)):
		line = lines[i].strip()
		values = line.split(",")
		#print values
		if values[p_i] != '' :
			accName = values[0].strip()
			accessions.append(accName)
			phenValues.append(float(values[p_i].strip()))


	print accessions
	print phenValues


	#Insert data into db.
	import MySQLdb
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Inserting data"

	print phenName

	sql_statement = "INSERT INTO stock_250k.phenotype_method  (id, short_name, only_first_96, biology_category_id, method_description, data_description, data_type) VALUES (" + str(m_i) + ", '" + phenName + "', true, 3, '" + method_description + "', '" + data_description + "', 'quantitative')"
	numRows = int(cursor.execute(sql_statement))
	row = cursor.fetchone()
	if row:
		print row

	for i in range(0, len(accessions)):
		val = phenValues[i]
		e_i = accessions[i]
		sql_statement = "INSERT INTO stock_250k.phenotype_avg (ecotype_id, value, ready_for_publication, method_id, transformed_value) VALUES ( " + str(e_i) + ", " + str(val) + ", 0, " + str(m_i) + ", " + str(val) + " )"
		#print sql_statement
		numRows = int(cursor.execute(sql_statement))
		row = cursor.fetchone()
		if row:
			print row
	conn.commit()
	cursor.close ()
	conn.close ()



def _insert_pecinka_phen_into_db_():
	"""
	Insert Ales Pecinka's phenotypes into the papaya DB. 
	"""
	#phen_file = env.home_dir+"/Projects/Data/phenotypes/UV_survial_phen.csv"
	phen_file = '/Users/bjarnivilhjalmsson/Projects/Ales_Pecinka/phen_pecinka_170310.tsv'
	phend = readPhenotypeFile(phen_file)
	phend.removePhenotypeIDs([1, 2, 3, 4])
	for i in [1, 2, 3, 4]:
		print phend.getPhenotypeName(i)
	phenotype_scoring = ""
	method_description = "UV effects on roots"
	growth_condition = ""
	biology_category_id = 19
	citations = ""
	data_description = ""
	transformation_description = "None"
	#transformation_description = "Log(SD/10+x-minVal)"
	method_id = None
	data_type = ""
	comment = "Data from Ales Pecinka (GMI)"
	assert len(phend.accessions) == len(phend.phenotypeValues), \
			       "length of phenotype values doesn't equal the number of accessions"
	phend.insert_into_DB(phenotype_scoring=phenotype_scoring, method_description=method_description,
			     growth_condition=growth_condition, biology_category_id=biology_category_id,
			     citations=citations, data_description=data_description,
			     transformation_description=transformation_description, method_id=method_id,
			     data_type=data_type, comment=comment)


def _insert_wilzcek_phen_into_db_():
	"""
	Insert Amity Wilzcek's phenotypes into the papaya DB. 
	"""
	#phen_file = env.home_dir+"/Projects/Data/phenotypes/phen_wilzcek_042710.tsv"
	phen_file = env['home_dir'] + "Projects/Amity_Wilczek/Wilczek_PTU_eID_new.csv"
	phend = readPhenotypeFile(phen_file, missing_value_decoder={'':'NA'}, with_db_ids=False)
	s = set(phend.accessions)
	#s.remove('')
	phend.filterAccessions(list(s))
	phenotype_scoring = ""
	method_description = "Various flowering times and date to bolting (PTU)"
	growth_condition = "Field growth conditions"
	biology_category_id = 15
	citations = ""
	data_description = ""
	transformation_description = "None"
	#transformation_description = "Log(SD/10+x-minVal)"
	method_id = None
	data_type = ""
	comment = "Data from Amity Wilczek"
	print phend.accessions
	print phend.phenotypeValues
	assert len(phend.accessions) == len(phend.phenotypeValues), "ARRG"
	print phend.phenIds
	phend.insert_into_DB(phenotype_scoring=phenotype_scoring, method_description=method_description,
			     growth_condition=growth_condition, biology_category_id=biology_category_id,
			     citations=citations, data_description=data_description,
			     transformation_description=transformation_description, method_id=method_id,
			     data_type=data_type, comment=comment)



def _insert_peijin_phen_into_db_():
	"""
	Insert Peijin Li's phenotypes into the papaya DB. 
	"""
	phen_file = env.home_dir + "Projects/Data/phenotypes/FLC_phenotypes_011710.tsv"
	phend = readPhenotypeFile(phen_file)
	phenotype_scoring = ""
	method_description = "FLC expression"
	growth_condition = ""
	biology_category_id = 1
	citations = ""
	data_description = ""
	transformation_description = "None"
	#transformation_description = "Log(SD/10+x-minVal)"
	method_id = None
	data_type = ""
	comment = "Data from Peijin Li"
	phend.insert_into_DB(phenotype_scoring=phenotype_scoring, method_description=method_description,
			     growth_condition=growth_condition, biology_category_id=biology_category_id,
			     citations=citations, data_description=data_description,
			     transformation_description=transformation_description, method_id=method_id,
			     data_type=data_type, comment=comment)



def _insert_bergelsson_phen_into_db_():
	"""
	Insert Joy Bergelsson's phenotypes into the DB.
	"""

	phen_file = env.home_dir + "Projects/Joy_Bergelsson/phen_bergelsson_051710.tsv"
	phend = readPhenotypeFile(phen_file)
	phenotype_scoring = ""
	method_description = "Glucosinolate levels"
	growth_condition = ""
	biology_category_id = 21
	citations = ""
	data_description = ""
	transformation_description = "None"
	#transformation_description = "Log(SD/10+x-minVal)"
	method_id = None
	data_type = ""
	comment = "Data from Joy Bergelsson"
	print len(phend.accessions), len(set(phend.accessions))


	phend.insert_into_DB(phenotype_scoring=phenotype_scoring, method_description=method_description,
			     growth_condition=growth_condition, biology_category_id=biology_category_id,
			     citations=citations, data_description=data_description,
			     transformation_description=transformation_description, method_id=method_id,
			     data_type=data_type, comment=comment)


#def _insert_wilczek_phen_into_db_2_():
#	"""
#	"""
#	
#	phen_file = "/Users/bjarni.vilhjalmsson/Projects/Data/phenotypes/wilczek_seaset.txt"
#	phend = readPhenotypeFile(phen_file)
#	phenotype_scoring = ""
#	method_description = ""
#	growth_condition = ""
#	biology_category_id = 22
#	citations = ""
#	data_description = ""
#	transformation_description = "None"
#	#transformation_description = "Log(SD/10+x-minVal)"
#	method_id = None
#	data_type = ""
#	comment = "Data from Amity Wilczek"			
#	print len(phend.accessions),len(set(phend.accessions))
#	
#	
#	phend.insert_into_DB(phenotype_scoring=phenotype_scoring, method_description=method_description, 
#			     growth_condition=growth_condition, biology_category_id=biology_category_id, 
#			     citations=citations, data_description=data_description, 
#			     transformation_description=transformation_description, method_id=method_id, 
#			     data_type=data_type, comment=comment)
#	
#	return phend
#
#	
#
#
#def _createRatioPhenotype_(pi1,pi2,methodID,short_name,method_description,data_description="",data_type="quantitative",created_by="bvilhjal",comment="",biology_category_id=1,only_first_96=0,readyForPublication = 0,user = "bvilhjal",host="papaya.usc.edu",passwd="bamboo123"):
#	phed = getPhenotypes(user=user, passwd=passwd,rawPhenotypes=True) 
#	pi1 = phed.getPhenIndex(pi1)
#	pi2 = phed.getPhenIndex(pi2)
#	print len(phed.accessions)
#	
#	newVals = []
#	for i in range(0,len(phed.accessions)):
#		newVal = 'NA'
#		if phed.phenotypeValues[i][pi1] != 'NA' and phed.phenotypeValues[i][pi2]!='NA':
#			v1 = float(phed.phenotypeValues[i][pi1])
#			v2 = float(phed.phenotypeValues[i][pi2])
#			#print v1, v2
#			if v2 != 0.0:
#				newVal = v1/v2
#		newVals.append(newVal)
#
#	import MySQLdb
#	print "Connecting to db, host="+host
#	if not user:
#		import sys
#		sys.stdout.write("Username: ")
#		user = sys.stdin.readline().rstrip()
#	if not passwd:
#		import getpass
#		passwd = getpass.getpass()
#	try:
#		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
#	except MySQLdb.Error, e:
#		print "Error %d: %s" % (e.args[0], e.args[1])
#		sys.exit (1)
#	cursor = conn.cursor ()
#
#	sqlStat = """INSERT INTO stock_250k.phenotype_method 
#    (id,short_name,only_first_96,biology_category_id,method_description, data_description, comment, created_by, data_type) 
#    VALUES """
#	sqlStat += "("+str(methodID)+",'"+short_name+"',"+str(only_first_96)+","+str(biology_category_id)+",'"+method_description+"','"+data_description+"','"+comment+"','"+created_by+"','"+data_type+"')"	
#	#print sqlStat	
#	#numRows = int(cursor.execute(sqlStat))
#	#print "Inserted data into stock_250k.phenotype_method:",numRows
#	#row = cursor.fetchone()
#	#if row:
#	#	print row
#	
#	for i in range(0,len(phed.accessions)):
#		
#		val = newVals[i]
#		if val !='NA':
#			e_i = phed.accessions[i]
#			sqlStatement = "INSERT INTO stock_250k.phenotype_avg (ecotype_id, value, ready_for_publication, method_id) VALUES ( "+str(e_i)+", "+str(val)+", "+str(readyForPublication)+", "+str(methodID)+")"
#			#print sqlStatement
#			numRows = int(cursor.execute(sqlStatement))
#			row = cursor.fetchone()
#			if row:
#				print row
#	print "Done inserting data into table!"
#	conn.commit()
#	cursor.close ()
#	conn.close ()
#
#
#	
#
#		
#def _runTest_():
#	filename = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_120308.tsv"
#	phed = readPhenotypeFile(filename)
#	print	phed._getIndexMapping_() #Modified 9/26/08
#	p_i = 39
#	print "phentoype:",p_i
#	print "Number of non-zero values:",phed.countValues(p_i)
#	phenIndex = phed.getPhenIndex(p_i)
#	vals = [] 
#	accessions = []
#	for i in range(0,len(phed.phenotypeValues)):
#		phenVals = phed.phenotypeValues[i]
#		if phenVals[phenIndex]!='NA':
#			vals.append(phenVals[phenIndex])
#			accessions.append(phed.accessions[i])
#	
#	print "len(phenVals):",len(vals)
#	print "len(accessions):",len(accessions)
#	print "len(set(accessions)):",len(set(accessions))
#	print vals
#	print accessions
#	
#
#	
#	
#def _runTest_():
#	filename = "/Network/Data/250k/dataFreeze_011209/250K_f13_012609.csv"
#	import dataParsers,snpsdata
#	snpsds = dataParsers.parseCSVData(filename, format=1, deliminator=",")#,debug=True)
#	snpsd = snpsdata.SNPsDataSet(snpsds,[1,2,3,4,5])
#	eDict = _getEcotypeIdToStockParentDict_()
#	accessions = _getFirst192Ecotypes_()
#	accessions = map(int,accessions)
#	a_set = set(accessions) 
#	accessions = map(int,snpsd.accessions)
#	a_set = a_set.intersection(set(accessions))
#	accessions = list(a_set)
#	accessions.sort()
#	print "ecotype_id, native_name, stock_parent"
#	i = 0
#	import csv
#	filename = "/Users/bjarni/tmp/192_accessions_031009.csv"
#	f = open(filename,"w")
#	w = csv.writer(f)
#	w.writerow(["ecotype_id","accession_name","stock_parent_id"])
#	#l = zip(log_true_pvals,sd_log_pvals)
#	for et in accessions:
#		w.writerow([et,eDict[et][0],eDict[et][1]])
#	f.close()
#
#	for et in accessions:
#		et = int(et)
#		print str(et)+", "+str(eDict[et][0])+", "+str(eDict[et][1])
#	print len(accessions)
#	
#def _runTest2_():
#	filename = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_120308.tsv"
#	phed = readPhenotypeFile(filename)	
#	phed.getTree()
#	
#def _runtest3_():
#	d = _getEcotypeIdInfoDict_()
#	tg_d = _getEcotype2TgEcotypeDict_()
#	import pickle
#	f = open("/Users/bjarnivilhjalmsson/tmp/eco_dict.pickle","w")
#	for ei in [6899,6909,6916,8215,6046,7515,7514,6968,6970,8268,8269,5719,7059,8281,8299,8324,8346,7336,7374,8400]:
#		print "%i,%i,%s"%(ei, tg_d[ei], d[tg_d[ei]])
#	pickle.dump(d,f)
#	f.close()
#	
#def _test_phen_map_plot_():
#	filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phenotypes_all_raw_081009.tsv"
#	phed = readPhenotypeFile(filename)	
#	accessions = phed.getAccessionsWithValues(1)
#	phed.plot_accession_map(1, "/Users/bjarnivilhjalmsson/tmp/test.pdf",accessions=accessions)
#	
#def _getTransformedPhenoype_():
#	filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phenotypes_all_raw_081009.tsv"
#	phed = readPhenotypeFile(filename)	
#	phed.logTransform(43)
#	phed.writeToFile("/tmp/test.tsv", [43], "\t")
#
#def get_AW_common_dataset():
#	filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_wilzcek_wo_OF_NS06_060210.csv"
#	phed = readPhenotypeFile(filename)  
#	phed.filter_accessions_w_missing_data()
#	phed.writeToFile('/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_wilzcek_wo_OF_NS06_060210.tsv', delimiter='\t')
#	 
#	     
#
#def get_hypocotyl_lenghts():
#	filename = "/home/GMI/bjarni.vilhjalmsson/Projects/data/phen_raw_071910.tsv"
#	phed = readPhenotypeFile(filename)	
#	phed.removePhenotypes([182])
#	phed.filter_accessions_w_missing_data()
#	phed.writeToFile('/tmp/hypocotyl.tsv',with_accession_names=True)


def combine_resistance_nc14():
	filename = "/Users/bjarni.vilhjalmsson/Projects/Data/phenotypes/phen_raw_092910.tsv"
	phed = readPhenotypeFile(filename)

	print phed.getPhenotypeName(1316)
	phed.filter_accessions_w_missing_data(1316)
	eids1 = phed.getNonNAEcotypes(1316)
	phen_vals1 = phed.getPhenVals(1316)

	phed = readPhenotypeFile(filename)
	print phed.getPhenotypeName(1328)
	phed.filter_accessions_w_missing_data(1328)
	eids2 = phed.getNonNAEcotypes(1328)
	phen_vals2 = phed.getPhenVals(1328)
	intersection_eids = set(eids1).intersection(set(eids2))
	for eid in intersection_eids:
		i1 = eids1.index(eid)
		i2 = eids2.index(eid)
		print phen_vals1[i1], phen_vals2[i2]
	new_ecotype_ids = []
	new_phen_vals = []
	for eid in set(eids1).difference(intersection_eids):
		i = eids1.index(eid)
		new_ecotype_ids.append(eid)
		new_phen_vals.append(phen_vals1[i])
	for eid in set(eids2).difference(intersection_eids):
		i = eids2.index(eid)
		new_ecotype_ids.append(eid)
		new_phen_vals.append(phen_vals2[i])
	for eid in intersection_eids:
		i = eids1.index(eid)
		new_ecotype_ids.append(eid)
		new_phen_vals.append(phen_vals1[i])
	print new_ecotype_ids
	phen_vals = [[pv] for pv in new_phen_vals]
	phed = PhenotypeData(new_ecotype_ids, ['Resistance_Nc14_combined_binary'], phen_vals, with_db_ids=False)
	phed.insert_into_DB([1], method_description='', growth_condition='Field', biology_category_id=2,
			citations='Jonathan Jones', data_description='Albug laibachii Nc14 combined',
			data_type='binary')



def combine_telomere_lengths():
	fn1 = env['phen_dir'] + 'telomere_lengths_192.csv'
	fn2 = env['phen_dir'] + 'telomere_lengths_swedish.csv'
	phed1 = readPhenotypeFile(fn1)
	phed2 = readPhenotypeFile(fn2)
	accs = list(set(phed1.accessions).union(set(phed2.accessions)))
	i_accs = list(set(phed1.accessions).intersection(set(phed2.accessions)))
	print 'i_accs:', i_accs
	pvs1 = phed1.getPhenVals(1)
	pvs2 = phed2.getPhenVals(1)
	pvs = []
	for acc in accs:
		if acc in i_accs:
			i1 = phed1.accessions.index(acc)
			i2 = phed2.accessions.index(acc)
			print acc, ':', pvs1[i1], pvs2[i2]
			pvs.append([(pvs1[i1] + pvs2[i2]) / 2.0])
		else:
			if acc in phed1.accessions:
				i = phed1.accessions.index(acc)
				pvs.append([pvs1[i]])
			elif acc in phed2.accessions:
				i = phed2.accessions.index(acc)
				pvs.append([pvs2[i]])
	print zip(accs, pvs)
	phed = PhenotypeData(accs, ['telomere_length'], pvs)
	phed.writeToFile(env['phen_dir'] + 'telomere_lengths_all.csv', delimiter=',')



def _get_fri_flc_():
	file_name = '/Users/bjarni.vilhjalmsson/test.phen_csv'
	filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_raw_100810.csv"
	phed = readPhenotypeFile(filename)
	phed.removePhenotypes([0, 4, 5, 43])
	phed.filter_na_ecotypes()
	phed.writeToFile(file_name)


if __name__ == '__main__':
	#_createRatioPhenotype_(184,183,222,"Trich_avg_JA_div_Trich_avg_C","Ratio: Trich_avg_JA/Trich_avg_C",biology_category_id=7)
	#pass
	#_getEcotypeIdToStockParentDict_()
	#eDict = _getEcotype2TgEcotypeDict_()
	#_insertPhenotypesIntoDb_()
	#_simpleInsertPhenotypesIntoDb_()
	#_getTransformedPhenoype_()
	#_getFirst192Ecotypes_()
	#_test_phen_map_plot_()
	#_insert_bergelsson_phen_into_db_()
	#get_AW_common_dataset()
	#get_hypocotyl_lenghts()
	#combine_telomere_lengths()
	get_all_phenotypes_from_db('/tmp/phen_raw_112210.csv')
	phed = parse_phenotype_file('/tmp/phen_raw_112210.csv')
	phed.convert_to_averages()
	phed.write_to_file('/tmp/phen_avg_112210.csv')
	print "Done!"

