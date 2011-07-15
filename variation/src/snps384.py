import dataParsers as dp
import sequences

def generate_384_snps_illumina_file():
        #sd = dp.parse_snp_data("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_data_t43_081009.csv")
        sd = dp.parse_snp_data("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_192_043009.csv")
        
        locus_names = []
        locus_sequences = []
        allele_1s = []
        allele_2s = []
        chromosomes = []
        positions = [] 
        col_alleles = []
        for ci,chromosome in enumerate([1,2,3,4,5]):
                col = sequences.get_col_sequence(chromosome)
                seq_len = len(col.seq)
                snpsd = sd.snpsDataList[ci]
                for i, pos in enumerate(snpsd.positions):
                        start_pos = max(0,pos-61)
                        end_pos = min(seq_len,pos+60)
                        snp = list(set(snpsd.snps[i]))
                        if len(snp)==2 and pos < len(col.seq) and '-' not in snp:
                                col_allele = col.seq[pos-1]
                                if col_allele==str(snp[0]):                                        
                                        other_allele = str(snp[1])
                                else:
                                        other_allele = str(snp[0])
                                        if col_allele!=str(snp[1]):
                                                raise Exception
                                allele_1s.append(col.seq[pos-3:pos+2])
                                allele_2s.append(col.seq[pos-3:pos-1]+other_allele+col.seq[pos:pos+2])
                                snp_str = '['+str(snp[0])+'/'+str(snp[1])+']'
                                local_seq = col.seq[start_pos:pos-1]+snp_str+col.seq[pos:end_pos]
                                locus_names.append('c'+str(chromosome)+'_p'+str(pos))
                                positions.append(pos)
                                locus_sequences.append(local_seq)
                                chromosomes.append(chromosome)
                                if not col.seq[pos-1] in snp:
                                        print col.seq[pos-1], snp,chromosome,pos
                                col_alleles.append(col.seq[pos-1])
        
        import csv
        w = csv.writer(open("/Users/bjarnivilhjalmsson/tmp/test.csv",'w'))
#        w.writerow(['Locus_Name','Target_Type','Sequence','Chromosome','Coordinate','Genome_Build_Version',
#                    'Source','Source_Version','Sequence_Orientation','Plus_Minus'])
#        w.writerow(['Chromosome','Coordinate','Allele_1','Allele_2','Genome_Build_Version','Sequence_Orientation','Plus_Minus'])
        w.writerow(['Chromosome','Coordinate','Allele_1','Allele_2','Genome_Build_Version','Sequence_Orientation','Plus_Minus'])
#        for (ln,ls,c,p) in zip(locus_names,locus_sequences,chromosomes,positions):
#                w.writerow([ln,'SNP',ls,c,p,'TAIR8','TAIR','8','Forward','Plus'])
        for (ln,a1,a2,c,p) in zip(locus_names,allele_1s,allele_2s,chromosomes,positions):
                w.writerow([c,p,a1,a2,'TAIR8','Forward'])

        print locus_sequences[100:110]
        print col_alleles[100:110]


def write_simple_toomaijan_file(filename="/Users/bjarnivilhjalmsson/tmp/test.csv" , window=25):
        sd = dp.parse_snp_data("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv")
        
        locus_names = []
        locus_sequences = []
        allele_1s = []
        allele_2s = []
        chromosomes = []
        positions = [] 
        col_alleles = []
        for ci,chromosome in enumerate([1,2,3,4,5]):
                col = sequences.get_col_sequence(chromosome)
                seq_len = len(col.seq)
                snpsd = sd.snpsDataList[ci]
                for i, pos in enumerate(snpsd.positions):
                        start_pos = max(0,pos-window-1)
                        end_pos = min(seq_len,pos+window)
                        snp = list(set(snpsd.snps[i]))
                        if len(snp)==2 and pos < len(col.seq) and '-' not in snp:
                                col_allele = col.seq[pos-1]
                                if col_allele==str(snp[0]):                                        
                                        other_allele = str(snp[1])
                                else:
                                        other_allele = str(snp[0])
                                        if col_allele!=str(snp[1]):
                                                raise Exception
                                allele_1s.append(col.seq[pos-3:pos+2])
                                allele_2s.append(col.seq[pos-3:pos-1]+other_allele+col.seq[pos:pos+2])
                                snp_str = '['+col_allele+'/'+other_allele+']'
                                local_seq = col.seq[start_pos:pos-1]+snp_str+col.seq[pos:end_pos]
                                locus_names.append('c'+str(chromosome)+'_p'+str(pos))
                                positions.append(pos)
                                locus_sequences.append(local_seq)
                                chromosomes.append(chromosome)
                                if not col.seq[pos-1] in snp:
                                        print col.seq[pos-1], snp,chromosome,pos
                                col_alleles.append(col.seq[pos-1])
        
        import csv
        w = csv.writer(open(filename,'w'))
        w.writerow(['Chromosome','Coordinate','Sequence','Genome_Build_Version','Sequence_Orientation'])
        for (ls,c,p) in zip(locus_sequences,chromosomes,positions):
                w.writerow([c,p,ls,'TAIR8','Forward'])

        print locus_sequences[100:110]
        print col_alleles[100:110]
        

def load_illumina_results():
        file_name = "/Users/bjarnivilhjalmsson/Projects/384_SNPs/GGGTScoreResults.csv"
        f = open(file_name,'r')
        import csv
        r = csv.reader(f)
        line = r.next()
        while line[0]!='[DATA]':
                line = r.next()
        line = r.next()
        score_i = line.index('Final_Score')
        chr_pos_list = []
        quality_scores = []
        for row in r:
                chr_pos_list.append((int(row[3]),int(row[4])))
                try:
                        quality_scores.append(float(row[score_i]))
                except Exception, err_str:
                        print row
                        print err_str
                        quality_scores.append(0.0)
  
        return chr_pos_list,quality_scores
        

def remove_overlapping_snps():
        from bisect import bisect
        #loading perlegen data.
        import dataParsers as dp
        sd = dp.parse_snp_data("/Users/bjarnivilhjalmsson/Projects/Data/perlegen/perlegen_011609.csv")
        perl_chr_pos = sd.getChrPosList()        
        chr_pos_list,quality_scores = load_illumina_results()
        in_perlegen = []
        nearby_snp_counts = [] 
        for i, (chr,pos) in enumerate(chr_pos_list):
                j = bisect(perl_chr_pos,(chr,pos))
                if perl_chr_pos[j-1]!=(chr,pos):
                        in_perlegen.append(False)
                else:
                        in_perlegen.append(True)
                k = j-2
                (n_chr,n_pos) = perl_chr_pos[k]
                n_count = 0                
                while pos-n_pos<61 and n_chr==chr:
                        n_count +=1
                        k-=1
                        (n_chr,n_pos) = perl_chr_pos[k]
                        
                k = j
                (n_chr,n_pos) = perl_chr_pos[k]
                while n_pos-pos<61 and n_chr==chr:
                        n_count +=1
                        k+=1
                        (n_chr,n_pos) = perl_chr_pos[k]
                nearby_snp_counts.append(n_count)
                if i%(len(chr_pos_list)/10)==0:
                        print '%d%% done.'%(((i+1.0)/len(chr_pos_list))*100)
        qc = zip(nearby_snp_counts,chr_pos_list,in_perlegen,quality_scores)
        qc.sort()
        k = bisect(qc,(1,(0,0),False,0))
        good_snp_chr_pos = []
        good_q_scores = []
        for n_count,chr_pos,in_perlegen,q_score in qc[:k]:
                if n_count==0 and in_perlegen:
                        good_snp_chr_pos.append(chr_pos)
                        good_q_scores.append(q_score)
        print len(good_snp_chr_pos)
        return good_snp_chr_pos,good_q_scores

def removing_imputed_snps():
        from bisect import bisect
        chr_pos_list,quality_scores = remove_overlapping_snps()
        sd = dp.parse_snp_data("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t52.csv")
        sd.filter_na_snps()
        sd_chr_pos = sd.getChrPosList()
        new_qs = []
        new_chr_pos = []
        for i, (chr,pos) in enumerate(chr_pos_list):
                j = bisect(sd_chr_pos,(chr,pos))
                if sd_chr_pos[j-1]!=(chr,pos):
                        if quality_scores[i]>0.8:
                                new_chr_pos.append((chr,pos))
                                new_qs.append(quality_scores[i])
        print len(new_chr_pos),len(new_qs)
        return new_chr_pos,new_qs
        
def output_validated_snps(filename):
        import csv
        f = open(filename,"w")
        chr_pos,qs = removing_imputed_snps()
        f.write("chromosome,position,score\n")
        for (chr,pos),score in zip(chr_pos,qs):
                f.write("%d,%d,%f\n"%(chr,pos,score))
        f.close()         
        
        

if __name__=='__main__':
        #generate_384_snps_illumina_file()
        write_simple_toomaijan_file()
        #output_validated_snps("/tmp/test.csv")
        