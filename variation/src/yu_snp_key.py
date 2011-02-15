"""
Various keys to converting SNP to and from Yu's format.
"""

# 2009-6-5 add lower-case letter to the dictionary
nt_2_number = {'|': -2,   #2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
        '-': -1,        #deletion
        'N': 0,
        'NA': 0,
        'n': 0,
        'na': 0,
        'X': 0,
        'x': 0,
        None: 0,
        'A': 1,
        'a': 1,
        'C': 2,
        'c': 2,
        'G': 3,
        'g': 3,
        'T': 4,
        't': 4,
        'AC': 5,
        'CA': 5,
        'M': 5,
        'm': 5,
        'AG': 6,
        'GA': 6,
        'R': 6,
        'r': 6,
        'AT': 7,
        'TA': 7,
        'W': 7,
        'w': 7,
        'CG': 8,
        'GC': 8,
        'S': 8,
        's': 8,
        'CT': 9,
        'TC': 9,
        'Y': 9,
        'y': 9,
        'GT': 10,
        'TG': 10,
        'K': 10,
        'k': 10
        }


number_2_UIPAC_nt = {-2: '|',       #2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
        -1: '-',
        0: 'NA',  #Missing data value.
        1:'A',
        2:'C',
        3:'G',
        4:'T',
        5:'M',
        6:'R',
        7:'W',
        8:'S',
        9:'Y',
        10:'K'
        }

