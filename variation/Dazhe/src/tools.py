# Binary search

def bufcount(filename):
    ' Buffer count as offered by Mykola Kharechko, arguably the fastest way to count lines in 2.6'
    f = open(filename)                  
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines

def bsearch(ls, element ,type="="):
    """
    Return the index of the element =,<=,<,>= or > than the element
    if such index does not exist return -1
    
    works for a SORTED list
    """
    li=-1
    ri=len(ls)    
    while (ri-li)>1:
        i=(li+ri)/2
        if element==ls[i]:
            li=i-1
            ri=i+1
            break
        elif element>ls[i]:
            li=i
        elif element<ls[i]:
            ri=i
    # set all to -1 if they do not exist
    if ri-li==1:
        i=-1
    if ri==len(ls):
        ri=-1
    if type=="=":
        return i
    elif type==">":
        return ri
    elif type=="<":
        return li
    elif type==">=":
        if i==-1:
            return ri
        else:
            return i
    elif type=="<=":
        if i==-1:
            return li
        else:
            return i
    else:
        return -1


