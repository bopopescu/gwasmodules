"""
An interface to the hapmap DB.
"""

import sys, csv
from elixir import *
from sqlalchemy import create_engine

metadata.bind = "mysql://root:blazroca@localhost/hapmap3"
#metadata.bind.echo = True


pop_abbrv_map = {
'ASW': ("African ancestry in Southwest USA",'USA'),
'CEU': ("Utah residents with Northern and Western European ancestry from the CEPH collection",'USA'),
'CHD': ("Chinese in Metropolitan Denver, Colorado",'USA'),
'GIH': ("Gujarati Indians in Houston, Texas",'USA'),
'LWK': ("Luhya in Webuye, Kenya",'Kenya'),
'MEX': ("Mexican ancestry in Los Angeles, California",'USA'),
'CHB': ("Han Chinese in Beijing, China",'China'),
'JPT': ("Japanese in Tokyo, Japan",'Japan'),
'MKK': ("Maasai in Kinyawa, Kenya",'Kenya'),
'TSI': ("Toscani in Italia",'Italia'),
'YRI': ("Yoruba in Ibadan, Nigeria",'Nigeria'),
}

data_dir="/Users/bjarnivilhjalmsson/Projects/Data/hapmap3/hapmap3_r2_b36_fwd.qc.poly/"

class Family(Entity):
	using_options(tablename='family') 
    	
    	id = Field(Unicode(30), primary_key=True)
	population = ManyToOne('Population')
	individuals = OneToMany('Individual')
    	description = Field(UnicodeText)

	def __repr__(self):
		return '<Family "%s" >' % (self.id)

class Population(Entity):
	using_options(tablename='population') 
	abbreviation=Field(Unicode(30), primary_key=True)
	sample_country = Field(UnicodeText)
	individuals = OneToMany('Individual')
	families = OneToMany('Family')
    	description = Field(UnicodeText)
    	map_file = Field(UnicodeText)
    	ped_file = Field(UnicodeText)
	def __repr__(self):
		return '<Population "%s", "%s">' % (self.abbreviation, self.description)

class Individual(Entity):
	using_options(tablename='individual') 
	id = Field(Unicode(30), primary_key=True)
	paternal_id = Field(Unicode(30))
    	maternal_id = Field(Unicode(30))
    	sex = Field(Unicode(30))
    	family = ManyToOne('Family')
    	population = ManyToOne('Population')
     	description = Field(UnicodeText)

	def __repr__(self):
		return '<Individual "%s", "%s", "%s", "%s">' % (self.id, self.sex, self.family, self.population)

class Marker(Entity):
	using_options(tablename='marker')
	id = Field(Unicode(30),primary_key=True)
	chromosome = Field(Unicode(10))
	position = Field(Integer)
	genetic_distance = Field(Float) #Morgans
	population_count = Field(Integer)

	def __repr__(self):
		return '<Marker "%s", "%s", "%s", "%s">' % (self.chromosome, self.position, self.genetic_distance, self.snp_id)


def load_family_data():
	r = csv.reader(open(data_dir+"relationships_w_pops_051208.txt"), delimiter='\t', quotechar='|')
	r.next()
	abbrvs = []
	populations = []
	for abbrv in pop_abbrv_map:
		(desc,country) = pop_abbrv_map[abbrv]
		m_file = "hapmap3_r3_b36_fwd."+abbrv+".qc.poly.map"
		p_file = "hapmap3_r3_b36_fwd."+abbrv+".qc.poly.ped"
		pop = Population(abbreviation=unicode(abbrv),description=unicode(desc),sample_country=unicode(country),ped_file=unicode(p_file),map_file=unicode(m_file))
		populations.append(pop)
		pop_abbrv_map[abbrv] = pop 
	session.commit()			
	for row in r:
		(fid,iid,pid,mid,sex,phem,pop_abr) = row
		if pid=='0':
			pid = None
		if mid=='0':
			mid = None
		pop = pop_abbrv_map[pop_abr]
		#try: 
		f = Family.get_by(id=unicode(fid))
		#except Exception, err_str:
		#	print err_str
			
		if not f:
			f = Family(id=unicode(fid),population=pop)
		Individual(id = unicode(iid),parental_id=unicode(pid),maternal_id=unicode(mid),
					sex=unicode(sex),population=pop,family=f)
		
		session.commit()
	
	
def load_marker_data(chromosome):
	
	#Now Marker Data...
	populations = Population.query.filter_by().all()
	for pop in populations:
		print "Looking at markers in population:",pop
		r = csv.reader(open(data_dir+pop.map_file), delimiter='\t', quotechar='|')
		for (chrs,m_id,gen_dist,pos) in r:
			if chrs==chromosome:
				m = Marker.get_by(id=unicode(m_id))
				if not m:
					gen_dist = float(gen_dist)
					if gen_dist==0:
						gen_dist=None
					m = Marker(id=unicode(m_id),chromosome=unicode(chrs),
						position=int(pos),genetic_distance=gen_dist,population_count=1)
				else:
					m.population_count += 1
				
		session.commit()		
	

#
def 


if __name__=='__main__':
	argv = sys.argv
	if 'load_data_to_db' in argv:
		setup_all()
		create_all()
		load_marker_data('3')
		#r = load_family_data()
		#INSERT INTO DB!!
