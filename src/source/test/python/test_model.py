import nupack as nu

def test_loop_energy():
    em = nu.Model(ensemble='some-nupack3', material='rna95-nupack3')
    assert abs(em.beta - 1.6225317071110157) < 1e-8
    assert em.hairpin_energy('AAAAT') == 4.1

def test_stack_energy():
    mod = nu.Model()
    pairs = 'AT CG GC TA GT TG'.split()
    pairs = 'AT CG GC TA'.split()
    s = ''
    for i in pairs:
        for j in pairs:
            seqs = [i[0] + j[0], j[1] + i[1]]
            s += '%d ' % int(100*mod.loop_energy(seqs))
        s += '\n'
    return s

def test_structure_energy():
    model = nu.Model(ensemble='some-nupack3', material='rna95-nupack3')
    assert model.structure_energy('CCCCTTTGGGG', '((((...))))') == -4.6
    assert model.structure_energy('GGAAACC', '.(...).') == 3.3
    assert abs(nu.Model(ensemble='stacking', material='rna95-nupack3').structure_energy('GGAAACC', '.(...).') - 2.7239996321033035) < 1e-8
