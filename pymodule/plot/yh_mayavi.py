#!/usr/bin/env python

from enthought.mayavi.tools.helper_functions import Pipeline, document_pipeline, bar_mode_dict
from enthought.mayavi.tools import tools
from enthought.traits.api import Array, Callable, CFloat, HasTraits, \
	List, Trait, Any, Instance, TraitError

from enthought.mayavi.tools.sources import vector_scatter, vector_field, scalar_scatter, \
			scalar_field, line_source, array2d_source, grid_source, \
			triangular_mesh_source, vertical_vectors_source

from enthought.mayavi.tools.modules import VectorsFactory, StreamlineFactory, GlyphFactory, \
	IsoSurfaceFactory, SurfaceFactory, ContourSurfaceFactory, \
	ImageActorFactory, glyph_mode_dict

class CustomBarChart(Pipeline):
	"""
	2012.2.21
		custom version of BarChart. It has two more keyword arguments, x_scale, y_scale,
			which is in charge of the lateral_scale in the X and Y direction.
		
	Plots vertical glyphs (like bars) scaled vertical, to do
	histogram-like plots.

	This functions accepts a wide variety of inputs, with positions given
	in 2-D or in 3-D.

	**Function signatures**::

		barchart(s, ...)
		barchart(x, y, s, ...)
		barchart(x, y, f, ...)
		barchart(x, y, z, s, ...)
		barchart(x, y, z, f, ...)

	If only one positional argument is passed, it can be a 1-D, 2-D, or 3-D
	array giving the length of the vectors. The positions of the data
	points are deducted from the indices of array, and an
	uniformly-spaced data set is created.

	If 3 positional arguments (x, y, s) are passed the last one must be
	an array s, or a callable, f, that returns an array. x and y give the
	2D coordinates of positions corresponding to the s values. 
	
	If 4 positional arguments (x, y, z, s) are passed, the 3 first are
	arrays giving the 3D coordinates of the data points, and the last one
	is an array s, or a callable, f, that returns an array giving the
	data value.
	"""

	_source_function = Callable(vertical_vectors_source)

	_pipeline = [VectorsFactory, ]

	mode = Trait('cube', bar_mode_dict,
					desc='The glyph used to represent the bars.')

	lateral_scale = CFloat(0.9, desc='The lateral scale of the glyph, '
				'in units of the distance between nearest points')

	def __call__(self, *args, **kwargs):
		""" Override the call to be able to scale automaticaly the axis.
		"""
		g = Pipeline.__call__(self, *args, **kwargs)
		gs = g.glyph.glyph_source
		# Use a cube source for glyphs.
		if not 'mode' in kwargs:
			gs.glyph_source = gs.glyph_list[-1]
		# Position the glyph tail on the point.
		gs.glyph_position = 'tail'
		gs.glyph_source.center = (0.0, 0.0, 0.5)
		g.glyph.glyph.orient = False
		if not 'color' in kwargs:
			g.glyph.color_mode = 'color_by_scalar'
		if not 'scale_mode' in kwargs:
			g.glyph.scale_mode = 'scale_by_vector_components'
		g.glyph.glyph.clamping = False
		x, y, z = g.mlab_source.x, g.mlab_source.y, g.mlab_source.z
		scale_factor = g.glyph.glyph.scale_factor* \
					tools._min_axis_distance(x, y, z)
		x_scale = kwargs.pop('x_scale', self.lateral_scale)
		y_scale = kwargs.pop('y_scale', self.lateral_scale)
		try:
			g.glyph.glyph_source.glyph_source.y_length = \
					y_scale/(scale_factor)
			g.glyph.glyph_source.glyph_source.x_length = \
					x_scale/(scale_factor)
		except TraitError:
			" Not all types of glyphs have controlable y_length and x_length"

		return g
	
	def get_all_traits(self):
		"""
		2012.2.21
			this function determines the set of possible keys in kwargs.
			add two keywords, x_scale, y_scale.
			
		Returns all the traits of class, and the classes in the pipeline.
		"""
		#call the parental one first.
		traits = Pipeline.get_all_traits(self)
		
		traits['x_scale'] = Trait(0.9, desc='The scale of the glyph on the x-axis, '
				'in units of the distance between nearest points')
		traits['y_scale'] = Trait(0.9, desc='The scale of the glyph on the y-axis, '
				'in units of the distance between nearest points')
		return traits

customBarchart = document_pipeline(CustomBarChart())

def testVerticalSurface():
	import pylab, numpy
	pylab.clf()
	
	no_of_phenotypes = 5
	x, y = numpy.mgrid[0:2*no_of_phenotypes:1, 0:10:1]	#added a gap of 1 column between two phenotypes. one phenotype occupies two rows & two columns.
	
	
	#remove the gap in x & y
	needed_index_ls = [0,5]
	#for i in [0, 5]:
	#	for j in xrange(no_of_phenotypes):
	#		needed_index_ls.append(no_of_phenotypes*i+j)
	#for i in range(0, no_of_phenotypes):
	#	needed_index_ls.append(2*i)
		#needed_index_ls.append(3*i+1)
		#y[3*i+1][1]=2
	x = x[:, needed_index_ls]
	y = y[:, needed_index_ls]
	
	enrichment_matrix = numpy.ones(x.shape, numpy.float)
	enrichment_matrix[:,:] =10
	enrichment_matrix[0,0]=3
	
	from enthought.mayavi import mlab
	mlab.clf()
	#from pymodule.yh_mayavi import customBarchart
	bar = customBarchart(x, y , enrichment_matrix, x_scale=0.9, y_scale=4.5, opacity=1, mode='cube', color=(0,1,0), scale_factor=1.0)
	"""
	#mlab.ylabel("KW")
	#mlab.xlabel("Emma")
	#mlab.zlabel("Enrichment Ratio")
	from pymodule.DrawMatrix import get_font 
	font = get_font()
	
	for i in range(len(xlabel_ls)):
		label = xlabel_ls[i]
		char_width, char_height = font.getsize(label)	#W is the the biggest(widest)
		
		mlab.text(2*i, 0, label, z=0, width=char_width/1500.)	#min(0.0075*len(label), 0.04))
	
	"""
	s = numpy.zeros(x.shape, numpy.int)
	#s[0,1]=0.5
	surf = mlab.surf(x, y, s, opacity=0.6, extent=[-1, 2*no_of_phenotypes, -1, 10, 0.0,0.0])
	mlab.show()
	
if __name__ == '__main__':
	testVerticalSurface()