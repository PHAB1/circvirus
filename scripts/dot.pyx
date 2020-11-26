
cpdef seq_inversor(str seq):
    cdef str seq_inv = ""
    cdef str nuc
    for nuc in seq:
        if nuc == "A":
            seq_inv += "T"
        
        if nuc == "T":
            seq_inv += "A"
        
        if nuc == "C":
            seq_inv += "G"
        
        if nuc == "G":
            seq_inv += "C"
    
    return(seq_inv[::-1])

cpdef compare(str seq,str seq_test,int count,double Miss_pair):
    cdef list coordy = [] 
    cdef int len_seq
    cdef int miss
    cdef int len_seq_test 
    for len_seq in range(len(seq)):
        miss = 0 

        try:
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

cpdef dot_(str seq,str seqy,double Miss_pair):
    # List Seq to be tested
    cdef int seq_len = len(seq)-len(seq)%20-20
    
    # Lists
    cdef dict lenghts = {}
    cdef dict lenghts_inv = {} 
    cdef list list_lenghts = [] 
    cdef list list_lenghts_inv = []
    cdef list list_coordy = []
    
    cdef str seq_test
    cdef str seq_test_inv
    
    cdef int coordx
    cdef int coordy
    cdef int b
    cdef int x
    
    for coordx in range(0,seq_len,20):
        seq_test = seq[coordx:coordx+20]
        seq_test_inv = seq_inversor(seq_test)
        
        list_coordy = compare(seqy,seq_test,0,Miss_pair)
        list_coordy_inv = compare(seqy,seq_test_inv,0,Miss_pair)
         
        # Seq normal
        for coordy in list_coordy:
            
            b = coordy-coordx
            x = -b
    
            if x > -1:
                list_lenghts.append(x)

                if x in lenghts.keys():
                    if coordx < lenghts[x][1]:
                        lenghts[x] = [coordx,coordy]

                else:
                    lenghts[x] = [coordx,coordy]
        

        # Seq invert
        for coordy in list_coordy_inv:
            
            b = coordy+coordx
            x = b
            

            if x > -1:
                list_lenghts_inv.append(x)

                if x in lenghts_inv.keys():
                    if coordx < lenghts_inv[x][1]:
                        lenghts_inv[x] = [coordx,coordy]

                else:
                    lenghts_inv[x] = [coordx,coordy]

    
    return(lenghts,lenghts_inv,list_lenghts,list_lenghts_inv)


