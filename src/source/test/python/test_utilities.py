from nupack import *
import numpy as np

################################################################################

def test_strand():
    assert str(RawStrand('A4')) == 'AAAA'
    assert str(RawStrand('AAAA')) == 'AAAA'
    assert RawStrand('A4') == RawStrand('AAAA')
    a = RawStrand('AAAA')
    assert a == RawStrand(a)
    a2 = TargetStrand([Domain('AAAA', name='blah')], name='blah2')

def test_domain():
    d = Domain('A4', name='blah')
    e = Domain('A4', name='blah2')
    assert d == d
    assert d != e

def test_complex():
    assert RawComplex('A+A') == RawComplex(['A', 'A'])
    assert RawComplex([RawStrand('A'), RawStrand('T')]) == RawComplex('A+T')

################################################################################

def test_pfunc():
    pf, fe = pfunc(['AAAA', 'CCCC'], model=Model())
    assert pf == 0
    assert fe == float('inf')

    pf, fe = pfunc(['AAAA'], model=Model())
    assert pf == 1
    assert float(pf.ln()) == 0
    assert fe == 0

################################################################################

def test_pairs():
    P = pairs(('AAAA', 'TTTT'), model=Model()).to_array()
    assert np.array_equal(P, P.T)
    assert np.all(abs(P.sum(0) - 1) < 1e-8)
    assert np.all(P >= 0)
    assert np.all(P <= 1)

################################################################################

def test_mfe():
    structures = mfe(['A4', 'T4'], model=Model())
    print(structures)

################################################################################

def test_energy():
    e = energy(['AAAAT'], '(...)', model=Model())
    print(e)

################################################################################

def test_prob():
    p = structure_probability(['AAAAT'], '(...)', model=Model())
    print(p)

################################################################################

def test_count():
    n = ensemble_size(['AAAAT'], model=Model())
    print(n, type(n))

################################################################################

def test_subopt():
    s = subopt(['AAAATTTT'], energy_gap=1.5, model=Model())
    print(s)

################################################################################

def test_sample():
    s = sample(['AAAATTTT'], num_sample=100, model=Model())
    print(s)

################################################################################

def test_loop_energy():
    dGloop1 = loop_energy(['AA', 'TT'],
        model=Model(material='RNA', ensemble='stacking')) # stack energy

################################################################################

def test_structure_energy():
    dGstruc1 = structure_energy(['AAAA', 'TTTT'], structure='((((+))))',
        model=Model(material='DNA', celsius = 25, ensemble='stacking'))

################################################################################

def test_design():
    designed_sequences = des('(((+)))', model=Model())
    print(designed_sequences)
    # --> ['CCC', 'GGG']

# ################################################################################

def test_defect():
    my_defect = defect('(((+)))', ['CCC', 'GGG'], model=Model())
    print(my_defect)
    # --> ['CCC', 'GGG']

def test_seq_distance():
    assert seq_distance('AA+TT', 'AA+CC') == 2
    assert seq_distance('TT', 'CC') == 2
    assert seq_distance(['TT'], ['CC']) == 2


def test_struc_distance():
    assert struc_distance('.....', '(...)') == 2
    assert struc_distance('.7', '((.3))') == 4
