"""
2010-4-14
	a data structure by Dazhe to hold all types of variation data: SNP or indels (CNV)
"""
import tools # For bufcount
import time

class vardata:
    """
    This class is designed to encompass all types of variation data: SNP or indels (CNV)
    accessions (x): accession list
    data (x*y): [chr, pos]+[list of chars]
    2010.03.01 overhauled readfromfile with inclusion of ability to compound genotypes (for plink tped) and display progress
    """
    def __init__(self, filename="", copyaccfrom = None):
        self.accessions = []
        self.posvec = []
        self.data = []
        self.dict = {'A':'A','T':'T','G':'G','C':'C','-':'-','NA':'?','AG':'?','CT':'?','AT':'?','AC':'?','CG':'?','GT':'?'}
        self.revdict = {'A':'A','T':'T','G':'G','C':'C','-':'-','?':'NA'}
        self.acclen = 0
        if filename!="":
            self.readfromfile(filename)
        elif copyaccfrom != None:
            self.copyacc(copyaccfrom)
    
    def copyacc(self, vardata2):
        'Copy the accessions from vardata2'
        self.accessions = list(vardata2.accessions)
        self.acclen = len(self.accessions)
    
    def readfromfile(self, filename, format=1, displayprogress=1, sep=",", compound_genotype=0, header=1):
        'Initialize from a file, 1 for bjarni format'
        fi = open(filename)
        if displayprogress==1:
            print "Reading variation data from "+filename
            starting_time = time.time()
        if header == 1:
            line = fi.readline()
            u = line.strip().split(sep)
            self.accessions = u[2:]
            self.acclen = len(self.accessions)
        self.read_format = format
        self.read_sep = sep
        self.read_compound_genotype = compound_genotype
        if displayprogress == 1:
            total_num_lines = tools.bufcount(filename)
            print "%s used for line counting"%(time.time()-starting_time)
            report_multiple = max(total_num_lines/20, 1)
            num_lines_read = 0
            for line in fi:
                num_lines_read += 1
                if num_lines_read % report_multiple == 0:
                    print "%s%% (%s lines) read, %ss total used"%((num_lines_read*20/total_num_lines)*5, num_lines_read, time.time()-starting_time)
                self.data.append(self.__parse_data_line__(line))
        else:
            for line in fi:
                self.data.append(self.__parse_data_line__(line))      

    def __parse_data_line__(self, line):
        'Internal slave function to parse the data line for readfromfile'
        u = line.strip().split(self.read_sep)
        if self.read_format == 1:
            returned_data=[int(u[0]),int(u[1]),[self.dict[i] for i in u[2:]]]
        elif self.read_format == 3: #plink tped
            if self.read_compound_genotype == 0:
                returned_data=[int(u[0]),int(u[3]),[u[1]]+u[4:]]
            else: # this is for the case where plink tped file have multiple column for each allele in multiploidal dataset
                returned_data=[int(u[0]),int(u[3]),[u[1]]+self.__compound_genotype__(u[4:])]
        else:
            returned_data=[int(u[0]),int(u[1]),[i for i in u[2:]]]   
        return returned_data       

    def __compound_genotype__(self, genotype_list):
        'Internal, compounds a genotype in the case of plink tped format'
        new_genotypelist = []
        if len(genotype_list)%2 != 0:
            raise "Incorrect number of snp data in compound genotype list"
        for i in xrange(len(genotype_list)/2):
            new_genotypelist.append("".join(sorted([genotype_list[i],genotype_list[i+1]])))
        return new_genotypelist
    
    def searchpos(self,chr,pos):
        'Internal: finds the index position of the given chromosome number and base position, or the one closest to it downstream'
        if len(self.data)==0:
            return 0
        l = 0; r = len(self.data)
        while (r-l)>1:
            i=(l+r)/2
            if chr==self.data[i][0] and pos==self.data[i][1]:
                l=i
                r=i
                break
            elif chr>self.data[i][0]:
                l=i
            elif chr<self.data[i][0]:
                r=i
            elif pos>self.data[i][1]:
                l=i
            elif chr<self.data[i][1]:
                r=i
        if chr<=self.data[l][0] and pos<=self.data[l][1]:
            return l
        else:
            return r
        
    def convert_to_01(self):
        'Convert snp data into 0 and 1'
        pass
    
    def insert(self, entry):
        'Insert a line into the data; The accessions must match'
        if len(entry[2])!=len(self.accessions):
            print 'Incorrect number of accessions in data being added'
        else:
            i = self.searchpos(entry[0],entry[1])
            self.data.insert(i,entry)
            
    def expand_acc(self, moreaccs):
        'Expand the accession list to include new accessions from the moreaccs list'
        new_acc_no = len(set(moreaccs)-set(self.accessions))
        if new_acc_no == 0:
            print "Nothing new in the accessions"
        else:
            self.accessions += list(set(moreaccs)-set(self.accessions))
            for i in self.data:
                i[2]+=["?"]*new_acc_no # Treat as missing data
            self.acclen += new_acc_no
    
    def index_acc(self, accs):
        'Index the entries in accs in self.accessions, missing ones would not be indexed and receive -1!'
        indexlist = []
        for i in accs:
            try:
                indexlist.append(self.accessions.index(i))
            except ValueError:
                indexlist.append(-1)
        return indexlist
    
    def merge(self, vardata2, indicator = 0):
        'merge 2 vardata into a new one; the accessions are merged as well!'
        self.expand_acc(vardata2.accessions)
        indl = self.index_acc(vardata2.accessions)
        for i in vardata2.data:
            newdat = ["?"]*self.acclen
            for ind in xrange(len(i[2])):
                newdat[indl[ind]] = i[2][ind]
            if indicator == 0:
                self.insert([i[0],i[1],newdat])
            else:
                self.insert([i[0],i[1],newdat,1])
    
    def iterate(self):
        'iterate through the variation data'
        for i in self.data:
            yield i
    
    def trimaccession(self, newaccessions):
        'trim the data and accessions, leaving only those in newaccessions'     
        pass
    
    def output(self, filename = "", format = 'Bjarni'):
        'output the vardata in a format indicated'
        if filename == "":
            filename = "VarData.csv"
        outfile = open(filename, "w")
        if format == 'Bjarni':
            outfile.write("Chromosome, Position, "+",".join(self.accessions)+"\n")
            outfile.write("\n".join([",".join([str(i[0]),str(i[1])]+[self.revdict[str(x)] for x in i[2]]) for i in self.data]))
        elif format == 'DataOnly':
            outfile.write(",".join(self.accessions)+"\n")
            outfile.write("\n".join([",".join([str(x) for x in i[2]]) for i in self.data]))
        elif format == 'NoConversion':
            outfile.write("Chromosome,Position,"+",".join(self.accessions)+"\n")
            outfile.write("\n".join([",".join([str(i[0]),str(i[1])]+[str(x) for x in i[2]]) for i in self.data]))
            
    def filter_maf(self, threshold):
        'Removes all snp whose maf is below threshold; threshold can be given either as absolute number or fraction'
        if threshold < 1:
            threshold = int(self.acclen*threshold)
        for entry in self.data:
            if count_ma(entry[2])<threshold:
                self.data.remove(entry)
                
    def filter_mis(self, threshold):
        'Removes all var whose missing rate is higher or equal than threshold; threshold can be given either as absolute number or fraction'
        if threshold<1:
            threshold = int(self.acclen*threshold)
        newvardata = vardata(copyaccfrom = self)
        for entry in self.data:
            if entry[2].count("?")<threshold:
                newvardata.insert(entry)
        self.data = newvardata.data
            
def count_ma(snp):
    'Count the minor allele of the given variation, data only!'
    alphabet = ["A","T","G","C","-","0","1"]
    counts = [snp.count(i) for i in alphabet]
    counts.sort(reverse=True)
    return counts[1]