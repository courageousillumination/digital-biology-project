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
                a = line[49:52].strip()
                b = line[72:75].strip()
                if a.isdigit() and b.isdigit()  and a != a_prev and b != b_prev:
                    a_prev = a
                    b_prev = b
                    res.append((int(a), int(b)))
    return res
