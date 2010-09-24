def tair8_to_tair9_positions(positions,chromosome):
        """
        Converts TAIR 8 positions on the given chromosome to the corresponding TAIR 9 positions. 
        """
        import env
        #Read conversion manual.
        f = open(env.home_dir+"Projects/Data/tair/TAIR9_assembly_updates_relative_to_TAIR8_TAIR9_assemblies.csv","r")
        f.readline()
        r = csv.reader(f)
        offsets = []
        offset_positions = [] 
        curr_offset = 0
        for row in r:
                if row[0] == 'chr'+str(chromosome):
                        pos = int(row[2])
                        offset_positions.append(pos)
                        type = row[3]
                        
                        if type=="insertion":
                                curr_offset +=len(row[4])
                                offsets.append(curr_offset)
                        if type=="deletion":
                                if row[4][0]!='N':
                                        curr_offset -=len(row[4])
                                else:
                                        curr_offset -= 1
                                offsets.append(curr_offset)
        
        f.close()
        curr_offset = 0
        offset_i = 0
        new_positions = []
        for pos in positions:
                if pos >= offset_positions[offset_i]:
                        curr_offset = offsets[offset_i]
                        offset_i += 1
                new_positions.append(pos+curr_offset)
        
        #Verfiy conversion!
        return new_positions
                        