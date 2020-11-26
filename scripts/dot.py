def seq_inversor(seq):
    seq_inv = ""
    for nuc in seq:
        if nuc == "A":
            seq_inv = seq_inv + "T"
        
        if nuc == "T":
            seq_inv = seq_inv + "A"
        
        if nuc == "C":
            seq_inv = seq_inv + "G"
        
        if nuc == "G":
            seq_inv = seq_inv + "C"
    
    return(seq_inv[::-1])

def compare(seq,seq_test,count,Miss_pair):
    coordy = []
    for len_seq in range(len(seq)):
        miss = 0
        
        try: # Comparing each nucleotide 
            for len_seq_test in range(len(seq_test)):
                if seq[len_seq + len_seq_test] != seq_test[len_seq_test]:
                    miss +=1

                if miss/int(len(seq_test)) > Miss_pair:
                    break
                    

        except IndexError:            
            break

        if miss/int(len(seq_test)) <= Miss_pair:
            coordy.append(len_seq)
            count += 1
    
    return(coordy)

def dot_(seq,seqy,Miss_pair):
    # List Seq to be tested
    list_seq = [i*20 for i in range(1,(len(seq)+20)//20)]

    # Lists
    lenghts, lenghts_inv, list_lenghts, list_lenghts_inv, list_coordy = {}, {}, [], [], []

    for coordx in list_seq:
        coordx = coordx-20
        seq_test = seq[coordx:coordx+20]
        seq_test_inv = seq_inversor(seq_test)
        
        list_coordy = compare(seqy,seq_test,0,Miss_pair)
        list_coordy_inv = compare(seqy,seq_test_inv,0,Miss_pair)
        
        # Seq normal
        for coordy in list_coordy:
            
            b = coordy-coordx
            x = -b
    
            if x >= 0:
                list_lenghts.append(x)

                if x in lenghts.keys():
                    if coordx < lenghts[x][1]:
                        lenghts[x] = [coodx,coordy]

                else:
                    lenghts[x] = [coordx,coordy]
        
        # Seq invert
        for coordy in list_coordy_inv:
            
            b = coordy+coordx
            x = b
            
            if x >= 0:
                list_lenghts_inv.append(x)

                if x in lenghts_inv.keys():
                    if coordx < lenghts_inv[x][1]:
                        lenghts_inv[x] = [coordx,coordy]

                else:
                    lenghts_inv[x] = [coordx,coordy]

    

    return(lenghts,lenghts_inv)









