#!/usr/bin/env python

from variation.src.misc import *
import Stock_250kDB

db_250k = Stock_250kDB.Stock_250kDB(drivername='mysql', username='yh',
                                   password='yh324', hostname='gmi-ara-devel-be', database='stock_250k')
db_250k.setup(create_tables=False)


input_fname1 = '/Network/Data/250k/db/results/type_1/4023_results.tsv'  # LM_cofactor_G. orontii conidiophore_32
input_fname2 = '/Network/Data/250k/db/results/type_1/4116_results.tsv'  # LM_G. orontii conidiophore_32
gwr = GWA.subtractTwoGWAS(input_fname1, input_fname2)
output_fname_prefix = '/tmp/LM_cofactor-LM_G. orontii conidiophore_32'

GWA.drawGWANicer(db_250k, gwr, output_fname_prefix, min_value=None, ylim_type=2)

