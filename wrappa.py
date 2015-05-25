def get_dehydrons(file_name):
    """
    Given a wrapping file returns a list of a pairs of residue numbers that represent
    dehydrons.
    """

    a_prev = None
    b_prev = None
    res = []
    with open(file_name) as w:
        for line in w:
            if line.find("HB_") == 0:
                a = line[48:52].strip()
                b = line[71:75].strip()
                chain_a = line[46:47].strip()
                chain_b = line[69:70].strip()
                if a and b and a != a_prev and b != b_prev:
                    a_prev = a
                    b_prev = b
                    res.append(((a, chain_a), (b, chain_b)))
    return res
