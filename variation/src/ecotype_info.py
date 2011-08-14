"""
A class which contains various ecotype info related help functions.
"""
import sqlalchemy as sqla
engine = sqla.create_engine('mysql://bvilhjal:*rri_bjarni@usc@papaya.usc.edu/stock')  #What engine?
md = sqla.MetaData()
md.bind = engine

def get_ecotype_table():
        return sqla.Table('ecotype', md, autoload=True)

def get_ecotype_table():
        return sqla.Table('ecotype', md, autoload=True)
#def

def get_ecotype_info_dict():
        from sqlalchemy.orm import sessionmaker
        et = get_ecotype_table()
        Session = sessionmaker(bind=et.bind)
        s = Session()
        s.query(et.columns['ID'],et.columns['NAME'],et.columns['NATIVENAME'],et.columns['STOCKPARENT'],
                et.columns['LATITUDE'],et.columns['LONGITUDE']).all()[:20]
        e_dict = {}
        