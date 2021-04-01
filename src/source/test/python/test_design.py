from nupack import *
import pickle

################################################################################

def test_pickle():
    def test(x):
        assert pickle.loads(pickle.dumps(x)) == x

    s = RawStrand('AAA')
    test(s)

    c = RawComplex('AA+TT')
    test(c)

    a = Domain('AAAA', name='Domain a')
    test(a)

    A = TargetStrand([a], name='RawStrand A')
    test(A)

    c1 = TargetComplex([A], structure='.4', name='RawComplex c1')
    test(c1)

    # t4 = TargetTube({c1: 2e-4}, max_size=2, name='t4')
    # test(t4)

################################################################################

def test_inputs():
    a = Domain('AAAA', name='Domain a')
    b = Domain('A4', name='Domain b') # equivalent sequence specification
    c = Domain('NNNNNNNNNN', name='Domain c')
    d = Domain('N10', name='Domain d') # equivalent sequence specification
    e = Domain('RRSSAAACCA', name='Domain e')
    f = Domain('R2S2A3C2A', name='Domain f') # equivalent sequence specification
    g = Domain('N10', name='Domain g')

        # Domains should not be specified inline
    A = TargetStrand([a, b, g], name='RawStrand A')
    B = TargetStrand([d, ~e], name='RawStrand B')
    C = TargetStrand([e, a, f], name='RawStrand C')
    D = TargetStrand([d, d, d], name='RawStrand D')


    c1 = TargetComplex([A], structure='.18', name='RawComplex c1')
    c2 = TargetComplex([A, B, B, C], structure='.18+.20+.20+.24', name='RawComplex c2')
    c3 = TargetComplex([A, A], structure='.18+.18', name='RawComplex c3')
    c4 = TargetComplex([A, A], structure='.18+.18', name='RawComplex c4')
    assert c3 == c4
    assert c3.name == 'RawComplex c3'
    assert c4.name == 'RawComplex c4'

    assert isinstance(c1[0], TargetStrand)
    try:
        c1[4]
        raise AssertionError
    except IndexError:
        pass

    # dot-parens-plus notation
    C1 = TargetComplex([A, B, C], structure='........((((((((((+))))))))))((((((((((+))))))))))..............', name='C1')

    # DU+ notation
    C2 = TargetComplex([D, D], structure='D30 +', name='C2')
    C3 = TargetComplex([B, B, B], structure='D10(D10 + D10 +)', name='C3')
    C4 = TargetComplex([B, A, B], structure='D8(U12 +) D10(+) U10', name='C4')

    # run-length encoded dot-parens-plus notation
    C5 = TargetComplex([B, C], structure='.%d+.%d' % (B.nt(), C.nt()), name='C5')

    # define target test tubes

    t1 = TargetTube({c1: 1e-8, c2: 1e-8}, SetSpec(3, [c3], [c1]), name='Tube t1')
    t2 = TargetTube({c1: 1e-8, c2: 1e-8}, SetSpec(3, [c3], [c1]), name='Tube t2')

    # specify tubes by their names and on-target complexes with on-target concentrations
    t1 = TargetTube({C1: 1e-6}, SetSpec(include=[C4, C5]), name='t1')

    # specify unnamed off-targets each denoted by a TargetStrand ordering
    t2 = TargetTube({C2: 1e-6}, SetSpec(include=[[D, D, D], [D, D, D, D]]), name='t2')

    # Mix the two kinds of specifications
    t3 = TargetTube({C1: 0.000001, C2: 1e-3}, SetSpec(include=[C4, [D, D, D]]), name='t3')

    # all complexes of up to 2 TargetStrands that are not on-targets in tube `T4'
    t4 = TargetTube({C1: 2e-4, C3: 3e-5}, SetSpec(max_size=2), name='t4')

    # specify off-targets as the sum of sets
    t5 = TargetTube({C4: 4e-6, C5: 5e-7}, SetSpec(max_size=2, include=[C3, [B, B, B, B]]), name='t5')

    # specify off-targets as the difference of sets
    t6 = TargetTube({C5: 6e-8}, SetSpec(max_size=3, exclude=[C3, [B, B]]), name='t6')

    # define hard constraints
    toeholds = ['CTAC', 'TAGT']
    gfp = 'auggugagcaagggcgaggagcuguucaccgggguggugcccauccuggucgagcuggacggcgacguaaacggccacaaguucagcguguccggcgagggcgagggcgaugccaccuacggcaagcugacccugaaguucaucugcaccaccggcaagcugcccgugcccuggcccacccucgugaccacccugaccuacggcgugcagugcuucagccgcuaccccgaccacaugaagcagcacgacuucuucaaguccgccaugcccgaaggcuacguccaggagcgcaccaucuucuucaaggacgacggcaacuacaag'

    hard = [
        # Match([c], [b, ~e]),
        # Match([a, b], [d, d, e]),
        # Complementarity(allow_wobble=True), # global flag (?)
        # Complementarity([a, b], [c, d, e], allow_wobble=True), # local flag (?)
        Similarity([b], 'S4', limits=[0.45, 0.55]), # GC content
        Library([a], catalog = [toeholds]),
        Window([a, ~b], sources = [gfp]),
        Pattern(['A5', 'C5', 'G5', 'U5'], scope=A),
        Pattern(['A5', 'C5', 'G5', 'U5'], scope=[b]),
        Pattern(['A4', 'C4', 'G4', 'U4', 'M6', 'K6', 'W6', 'S6', 'R6', 'Y6']),
        Diversity(word=4, types=2),
        Diversity(word=6, types=3),
        Diversity(word=10, types=4, scope=B)
    ]


    e = Domain('S2', name='e')
    f = Domain('S2', name='f')

    #add another constraint to the constrain set
    hard += [Complementarity([e],[f], wobble_mutations=True)]
    hard.append(Complementarity([e],[f], wobble_mutations=True)) # same thing

    a = Domain('N'*10, name='a')
    b = Domain('N'*4, name='b')
    c = Domain('H'*6, name='c')
    d = Domain('N'*6, name='d')
    e = Domain('S'*2, name='e')

    match1 = Match([c], [b, ~e])
    match2 = Match([a, b], [d, d, e])

    comp = Complementarity([a, b], [c, d, e])

    a = Domain('N'*10, name='a')
    sim1 = Similarity([a], 'S5K5', limits=[0.25, 0.75])

    # "composition constraint" special case: enforce 45-55% GC content
    b = Domain('N20', name='b')
    sim2 = Similarity([b], 'S20', limits=[0.45, 0.55])

    a = Domain('N10', name='a')
    b = Domain('N10', name='b')
    c = Domain('N10', name='c')
    e = Domain('N10', name='e')

    gfp = 'AUGGUGAGCAAGGGCGAGGAGCUGUUCACCGGGGUGGUGCCCAUCCUGGUCGAGCUGGACGGCGACGUAAACGGCCACAAGUUCAGCGUGUCCGGCGAGGGCGAGGGCGAUGCCACCUACGGCAAGCUGACCCUGAAGUUCAUCUGCACCACCGGCAAGCUGCCCGUGCCCUGGCCCACCCUCGUGACCACCCUGACCUACGGCGUGCAGUGCUUCAGCCGCUACCCCGACCACAUGAAGCAGCACGACUUCUUCAAGUCCGCCAUGCCCGAAGGCUACGUCCAGGAGCGCACCAUCUUCUUCAAGGACGACGGCAACUACAAG'

    rfp = 'CCUGCAGGACGGCGAGUUCAUCUACAAGGUGAAGCUGCGCGGCACCAACUUCCCCUCCGACGGCCCCGUAAUGCAGAAGAAGACCAUGGGCUGGGAGGCCUCCUCCGAGCGGAUGUACCCCGAGGACGGCGCCCUGAAGGGCGAGAUCAAGCAGAGGCUGAAGCUGAAGGACGGCGGCCACUACGACGCUGAGGUCAAGACCACCUACAAGGCCAAGAAGCCCGUGCAGCUGCCCGGCGCCUACAACGUCAACAUCAAGUUGGACAUCACCUCCCACAACGAGGACUACACCAUCGUGGAACAGUACGAACGCGCCGAGGGCCGCCACUCCACCGGCGGCAUGGACGAGCUGUACAAGUAA'

    # constrain window to be drawn from source
    window1 = Window([a, ~b], [gfp])
    # OR constrain window to be drawn from more than once source
    window2 = Window([~c, e], [gfp, rfp])

    a = Domain('N6', name='a')
    b = Domain('N12', name='b')

    # define a library of sequences
    toeholds = ['CAGUGG', 'AGCUCG', 'CAGGGC']

    # define a library of codons for each amino acid
    aaI = ['AUU', 'AUC', 'AUA']
    aaL = ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG']
    aaV = ['GUU', 'GUC', 'GUA', 'GUG']
    aaF = ['UUU', 'UUC']
    aaM = ['AUG']
    aaC = ['UGU', 'UGC']
    aaA = ['GCU', 'GCC', 'GCA', 'GCG']
    aaG = ['GGU', 'GGC', 'GGA', 'GGG']
    aaP = ['CCU', 'CCC', 'CCA', 'CCG']
    aaT = ['ACU', 'ACC', 'ACA', 'ACG']
    aaS = ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC']
    aaY = ['UAU', 'UAC']
    aaW = ['UGG']
    aaQ = ['CAA', 'CAG']
    aaN = ['AAU', 'AAC']
    aaH = ['CAU', 'CAC']
    aaE = ['GAA', 'GAG']
    aaD = ['GAU', 'GAC']
    aaK = ['AAA', 'AAG']
    aaR = ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
    aaSTOP = ['UAA', 'UAG', 'UGA']

    # domain a is drawn from the `toeholds' library
    lib1 = Library([a], [toeholds])

    # domain b is drawn from a concatenation of library sequences representing codons
    lib2 = Library([b], [aaI, aaM, aaC, aaG])

    a = Domain('N12', name='a')
    b = Domain('N12', name='b')
    A = TargetStrand([a, ~a], name='A')
    B = TargetStrand([b, ~b], name='B')

    # pattern prevention for a domain
    pat1 = Pattern(['A4', 'U4'], scope=[a])

    # pattern prevention for a TargetStrand
    pat2 = Pattern(['A4', 'U4'], scope=B)

    # preventing the same patterns for TargetStrand `A' and domain `b'
    pat3 = Pattern(['A5', 'C5', 'G5', 'U5'], scope=A)

    # global pattern prevention
    pat4 = Pattern(['A4', 'C4', 'G4', 'U4', 'M6', 'K6', 'W6', 'S6', 'R6', 'Y6'])

    div1 = Diversity(word=4, types=2)
    div2 = Diversity(word=6, types=3)

    div3 = Diversity(word=10, types=4, scope=[a])

    # define soft for soft constraints
    soft = [
        Pattern(['A4', 'U4'], scope=[a]),
        Pattern(['A5', 'C5', 'G5', 'U5'], scope=A), # default weight 1
        Pattern(['A4', 'C4', 'G4', 'U4', 'M6', 'K6', 'W6', 'S6', 'R6', 'Y6'], weight=0.5),
        Similarity([b], 'S12', limits=[0.45, 0.55], weight=0.25),
        EnergyMatch([a, b]), # min energy diff to median
        EnergyMatch([a, b], energy_ref=-17, weight=0.5) # energy diff to reference
    ]

    pat = Pattern(patterns=['A4', 'C4', 'G4', 'U4', 'M6',
        'K6', 'W6', 'S6', 'R6', 'Y6'], weight=0.5)

    a = Domain('N10', name='a')
    b = Domain('N20', name='b')

    # explicitly specify weight
    sim = Similarity([b], 'S20', limits=[0.45, 0.55], weight=0.25)

    a = Domain('N12', name='a')
    b = Domain('N12', name='b')
    A = TargetStrand([a, ~a], name='A')
    B = TargetStrand([b, ~b], name='B')

    C = TargetComplex([A], "(10.4)10", name='C')
    D = TargetComplex([A, A], "D24 +", name='D')

    ssm1 = SSM(scope=[C, D], word=4, weight=0.15)

    # the same complexes with larger windows weighted higher
    ssm2 = SSM(scope=[C, D], word=5, weight=0.25)
    ssm3 = SSM(scope=[C, D], word=6, weight=0.45)

    # equalize to median value
    diff1 = EnergyMatch([a, b])

    # equalize to reference value, with explicit weight
    diff2 = EnergyMatch([a, b], energy_ref=-17, weight=0.5)

    weights = Weights([t1, t2])

    weights[a] *= 2
    weights[:,:,:,t2] = 2
    weights[:, :, c1, t1] = 5,
    weights[b, A, :, :] = 0.75
    weights[a, D, c3, t2] = 0.5
    weights[d, :, :, t2] = 3

    print(weights[:, :, :, t1])
    print(weights[:, :, c2, t1])
    print(weights[:, A, c2, t1])
    print(weights[d, :, c2, t1])

    options = DesignOptions(
        seed=0,     # random number generation seed
        f_stop=0.5,  # stop condition
        f_passive=0.01,
        H_split=2,
        N_split=12,
        f_split=0.99,
        f_stringent=0.99,
        dG_clamp=-20,
        M_bad=300, # number of bad
        M_reseed=50,
        M_reopt=3,
        f_redecomp=0.03,
        f_refocus=0.03,
        f_sparse=1e-05,
        slowdown=0,
        log='',
        decomposition_log=None,
        thermo_log=None,
        time_analysis=False
    )


    # run the job
    try:
        des = Design(tubes=[t1, t2],
            hard_constraints=[], soft_constraints=[],
            defect_weights=None, options=options, model=Model()).run(trials=1)
    except RuntimeError:
        pass

    try:
        des = complex_design(complexes=[c1, c2],
            hard_constraints=[], soft_constraints=[],
            defect_weights=None, options=options, model=Model()).run(trials=1)
    except RuntimeError:
        pass

    a = Domain('N20', name='a')
    A = TargetStrand([a], name='A')
    B = TargetStrand([~a], name='B')

    C = TargetComplex([A, B], structure='(20+)20', name='C')

    tube = TargetTube({C: 1e-6}, SetSpec(2), name='tube1')

    des = Design([tube], model=Model(), options=options)
    result = des.run(trials=1)[0]

    check_print(result)

    print(result.to_analysis(tube))
    print(result.to_analysis(C))
    print(result.to_analysis(A))
    print(result.to_analysis(a))

    print(result.defects.ensemble_defect)
    print(result.defects.tubes)
    print(result.defects.complexes)
    print(result.defects.tube_complexes)

    #########

    designed = result.to_analysis(tube)
    # Compute the MFEs of the designed complexes that were in t1
    tube_results = tube_analysis(tubes=[designed], compute=['mfe'], model=Model())
    # Compute complex concentrations with a different set of strand concentrations
    conc_results = complex_concentrations(designed, result.analysis_result.complexes)

    conc_results = complex_concentrations(designed, result.analysis_result.complexes, {
        result.to_analysis(A): 1e-8, result.to_analysis(B): 1e-9
    })

    #######

    import pathlib
    out = pathlib.Path('design-result.o')
    result.save(out)
    result = DesignResult.load(out)
    out.unlink()

    newer_result = des.run(trials=1, restart=[result])[0]

    from nupack.design import TimeInterval, WriteToFileCheckpoint

    result = des.run(trials=1, checkpoint_condition=[TimeInterval(1)], checkpoint_handler=[WriteToFileCheckpoint("design-checkpoint")])

################################################################################

def test_example8():
    # set physical properties
    model = Model(material='RNA', celsius=23)

    # define domains
    a = Domain('ACCUCCAAGCACAACUGUGGCCCCAUA', name='a')
    b = Domain('GGGGCCGGAUUACAACUUUCCCUGUGAAC', name='b')
    c = Domain('AUCACAGACAGUUAACCACUUGAGG', name='c')
    d = Domain('AUCAAGUGGGCUUGGAGC', name='d')

    # define strands from domains
    left = TargetStrand([a], name='left')
    top = TargetStrand([b], name='top')
    right = TargetStrand([c], name='right')
    bottom = TargetStrand([d], name='bottom')

    # define complex compsed of strands in a given order AND
    # Define target structure for complex
    stickfigure = TargetComplex((left, top, right, bottom), structure="U2D8(U2D6(D6(U3+)D3U9D6(U2+U1))U2D8(U2+U1))U1", name='stickfigure')

    # define test tube
    figuretube = TargetTube({stickfigure: 1e-6}, SetSpec(3), name='figuretube')

    # evaluate
    design = Design([figuretube], model)
    result = design.evaluate()
    check_print(result)

################################################################################

def check_print(x):
    assert str(x)
    assert isinstance(x._repr_html_(), str)

def test_example7():
    # set physical parameters
    model = Model(material='rna', celsius=37)

    # define domains
    a = Domain('N6', name='a')
    c = Domain('N8', name='c')
    b = Domain('N4', name='b')
    w = Domain('N2', name='w')
    y = Domain('N4', name='y')
    x = Domain('N12', name='x')
    z = Domain('N3', name='z')
    s = Domain('N5', name='s')

    # define strands from domains
    Cout_s   = TargetStrand([w, x, y, s], name='Cout_s')
    A_s      = TargetStrand([~c, ~b, ~a, ~z, ~y], name='A_s')
    A_toe_s  = TargetStrand([~c], name='A_toe_s')
    C_s      = TargetStrand([w, x, y, s, ~a, ~z, ~y, ~x, ~w], name='C_s')
    C_loop_s = TargetStrand([s, ~a, ~z], name='C_loop_s')
    B_s      = TargetStrand([x, y, z, a, b], name='B_s')
    Xs_s     = TargetStrand([a, b, c], name='Xs_s')

    # define complexes composed of one or more strands in a given order AND
    # define target structures for each complex
    C      = TargetComplex([C_s],       'D2 D12 D4( U5 U6 U3 )', name='C')
    B      = TargetComplex([B_s],       'U12 U4 U3 U6 U4', name='B')
    C_loop = TargetComplex([C_loop_s],  'U14', name='C_loop')
    A_B    = TargetComplex([A_s, B_s],  'U8 D4 D6 D3 D4(+ U12)', name='A_B')
    X      = TargetComplex([Xs_s],      'U18', name='X')
    X_A    = TargetComplex([Xs_s, A_s], 'D6 D4 D8(+) U3 U4', name='X_A')
    C_out  = TargetComplex([Cout_s],    'U23', name='C_out')
    B_C    = TargetComplex([B_s, C_s],  'D12 D4 D3 D6 (U4 + U2 U12 U4 U5) U2', name='B_C')
    A_toe  = TargetComplex([A_toe_s],   'U8', name='A_toe')

    # on-target tubes
    Step_0 = TargetTube({C: 1e-08, X: 1e-08, A_B: 1e-08}, SetSpec(2, include=[[A_s], [B_s]], exclude=[X_A]), name='Step_0')

    Step_1 = TargetTube({X_A: 1e-08, B: 1e-08}, SetSpec(2, include=[X, A_B]), name='Step_1')

    Step_2 = TargetTube({B_C: 1e-08}, SetSpec(2, include=[B, C]), name='Step_2')

    # global orthogonality tube
    Crosstalk = TargetTube({
        A_B: 1e-08,
        C: 1e-08,
        X: 1e-08,
        B: 1e-08,
        C_out: 1e-08,
        C_loop: 1e-08,
        A_toe: 1e-08,
    }, SetSpec(2, exclude=[X_A, B_C, [Xs_s, A_toe_s], [B_s, C_loop_s]]), name='Crosstalk')

    # GC content constraints
    hard = [
        Similarity(d, 'S' * d.nt(), limits=(0.45, 0.55))
            for d in [Cout_s, A_s, C_s, C_loop_s, B_s, Xs_s]
    ]

    # sources lines
    tpm3 = 'GAACACTATTAGCTATTTGTAGTACTCTAAAGAGGACTGCAGAACGCATCGCAGTAGTGGTGAAAAGCCGTGCGTGCGCGTGAAACATCTGATCCTCACGTTACTTCCACTCGCTCTGCGTTTGACTTGTTGGCGGGGCGTTGGTGCCTTGGACTTTTTTTTCCTCCTTCTCTTCTTCGCGGCTCGGTCCACTACGCTGCTCGAGAGGAATCTGCTTTATTCGACCACACTACTCCTAAAGTAACACATTAAAATGGCCGGATCAAACAGCATCGATGCAGTTAAGAGAAAAATCAAAGTTTTACAACAGCAAGCAGATGAGGCAGAAGAAAGAGCCGAGATTTTGCAGAGACAGGTCGAGGAGGAGAAGCGTGCCAGGGAGCAGGCTGAGGCAGAGGTGGCTTCTCTGAACAGGCGTATCCAGCTGGTTGAGGAGGAGTTGGATCGTGCTCAGGAGAGACTGGCCACAGCCCTGCAAAAGCTGGAGGAAGCCGAGAAGGCCGCAGATGAGAGCGAGAGAGGGATGAAGGTGATTGAGAACAGGGCTCTGAAGGATGAGGAGAAGATGGAGCTGCAGGAGATCCAGCTTAAGGAGGCCAA'
    hard += [Window([a, b, c], [tpm3])]

    desm = 'CATTTACACAGCGTACAAACCCAACAGGCCCAGTCATGAGCACGAAATATTCAGCCTCCGCCGAGTCGGCGTCCTCTTACCGCCGCACCTTTGGCTCAGGTTTGGGCTCCTCTATTTTCGCCGGCCACGGTTCCTCAGGTTCCTCTGGCTCCTCAAGACTGACCTCCAGAGTTTACGAGGTGACCAAGAGCTCCGCTTCTCCCCATTTTTCCAGCCACCGTGCGTCCGGCTCTTTCGGAGGTGGCTCGGTGGTCCGTTCCTACGCTGGCCTTGGTGAGAAGCTGGATTTCAATCTGGCTGATGCCATAAACCAGGACTTCCTCAACACGCGTACTAATGAGAAGGCCGAGCTCCAGCACCTCAATGACCGCTTCGCCAGCTACATCGAGAAGGTGCGCTTCCTCGAGCAGCAGAACTCTGCCCTGACGGTGGAGATTGAGCGTCTGCGGGGTCGCGAGCCCACCCGTATTGCAGAGCTGTACGAGGAGGAGATGAGAGAGCTGCGCGGACAGGTGGAGGCACTGACCAATCAGAGATCCCGTGTGGAGATCGAGAGGGACAACCTAGTCGATGACCTACAGAAACTAAAGCTCAGACTTC'
    hard += [Window([w, x, y, z], [desm])]

    # Prevented patterns
    hard += [Pattern(['A4','C4','G4','U4'])]

    options = DesignOptions(seed=93, f_stop=0.1)
    tubes = [Step_0, Step_0, Step_0, Crosstalk]
    weights = Weights(tubes)

    check_print(weights)
    weights[:, :, :, Step_0] = 2
    weights[:, :, :, 'Step_0'] = 1

    mydes = Design(tubes, model=Model(), hard_constraints=hard, options=options, defect_weights=weights)
    result = mydes.run(trials=1)[0]

    check_print(result)
    check_print(result.defects.tubes)
    check_print(result.defects.tube_complexes)
    check_print(result.defects.complexes)
    check_print(result.to_analysis)

    opt = mydes.launch(trials=1)
    results = opt.wait()

################################################################################

def test_old_design():
    x = design.Specification(model=Model())
    n = x.domains.__len__()
    assert n == 0
    assert len(x.tubes) == 0
    x.add_domain('a', 'N'*12)
    x.add_domain('b', 'N'*24)
    x.add_domain('c', 'N'*12)

    x.add_strand('i1', ('b*', 'a*'))
    x.add_strand('h2', ('b*', 'a*', 'b', 'c'))
    x.add_strand('i2', ('c*', 'b*'))
    x.add_strand('h1', ('a', 'b', 'c*', 'b*'))


    model = Model(material='DNA', ensemble='some-nupack3', celsius=25)
    x.model = model

    x.add_complex('A', ('h1'), '.12(24.12)24')
    x.add_complex('B', ('h2'), '(24.12)24.12')
    x.add_complex('X', ('i1'), '.36')
    x.add_complex('X__A', ('h1', 'i1'), '(36.36+)36')
    x.add_complex('X__A__B', ('h1', 'h2', 'i1'), '(72+.36)36+)36')
    x.add_complex('X__A__A__B', ('h1', 'h1', 'h2', 'i1'), '(72+(36.36+)72+)36')
    x.add_complex('A_out', ('i2'), '.36')

    def on_targets(l):
        return {k: 1e-8 for k in l}

    def split(l):
        return [x.split() for x in l]

    x.add_tube('Step_0', on_targets(['X', 'A', 'B']))
    assert len(x.tubes) == 1
    exclude = split(['X__A__B', 'X__A', 'X__A__A__B', 'h1 h1 h2', 'h2 h2 h1', 'h1 h1 h2 h2'])
    explicit = split(['h1 h2'])
    x.add_off_targets('Step_0', max_size=2, explicit=explicit, exclude=exclude)

    x.add_tube('Step_1', on_targets(['X__A']))
    explicit = split(['h1', 'i1'])
    x.add_off_targets('Step_1', max_size=2, explicit=explicit)

    x.add_tube('Step_2', on_targets(['X__A__B']))
    explicit = split(['h2', 'X__A'])
    exclude = split(['h1 h1 h2', 'h2 h2 h1', 'h1 h1 h2 h2', 'X__A__A__B'])
    x.add_off_targets('Step_2', max_size=2, explicit=explicit, exclude=exclude)

    x.add_tube('Step_3', on_targets(['X__A__A__B']))
    explicit = split(['h1', 'X__A__B'])
    exclude = split(['h1 h1 h2', 'h2 h2 h1', 'h1 h1 h2 h2'])
    x.add_off_targets('Step_3', max_size=2, explicit=explicit, exclude=exclude)

    x.add_tube('Crosstalk', on_targets(['A', 'B', 'X', 'A_out']))
    exclude = split(['i1 h1', 'i2 h2', 'h1 h2 i1', 'h1 h1 h2 i1', 'h1 h1 h2', 'h2 h2 h1', 'h1 h1 h2 h2', 'i2 h1 h2', 'i2 h1 h2 h2'])
    x.add_off_targets('Crosstalk', max_size=2, exclude=exclude)

    # x.add_source("tmp3", "GAACACTATTAGCTATTTGTAGTACTCTAAAGAGGACTGCAGAACGCATCGCAGTAGTGGTGAAAAGCCGTGCGTGCGCGTGAAACATCTGATCCTCACGTTACTTCCACTCGCTCTGCGTTTGACTTGTTGGCGGGGCGTTGGTGCCTTGGACTTTTTTTTCCTCCTTCTCTTCTTCGCGGCTCGGTCCACTACGCTGCTCGAGAGGAATCTGCTTTATTCGACCACACTACTCCTAAAGTAACACATTAAAATGGCCGGATCAAACAGCATCGATGCAGTTAAGAGAAAAATCAAAGTTTTACAACAGCAAGCAGATGAGGCAGAAGAAAGAGCCGAGATTTTGCAGAGACAGGTCGAGGAGGAGAAGCGTGCCAGGGAGCAGGCTGAGGCAGAGGTGGCTTCTCTGAACAGGCGTATCCAGCTGGTTGAGGAGGAGTTGGATCGTGCTCAGGAGAGACTGGCCACAGCCCTGCAAAAGCTGGAGGAAGCCGAGAAGGCCGCAGATGAGAGCGAGAGAGGGATGAAGGTGATTGAGAACAGGGCTCTGAAGGATGAGGAGAAGATGGAGCTGCAGGAGATCCAGCTTAAGGAGGCCAAGCACATTGCTGAGGAGGCTGACCGCAAATATGAAGAGGTGGCTCGTAAGCTGGTGATCGTTGAGGGAGAGTTGGAGCGTACAGAGGAGAGAGCAGAGCTTGCAGAGAGCCATGTCAAGCAGATGGAGGAGGAGCTGAGAGCTCTTGACCAGACACTGAAGACTCTTCAGGCCTCAGAGGAGAAGTATTCCCAGAAGGAGGACAAGTATGAGGAAGAAATCAAGATCCTCACTGATAAGCTGAAGGAGGCTGAGACCCGTGCAGAGTTTGCTGAGAGGTCTGTGGCCAAACTGGAGAAAACCATTGATGATTTGGAAGAGAAACTGAGAGATGCTAAAGAGGAGAACATCAAGATCCATGCTACTTTGGACCAGACCCTGAGCGAGCTCAATAGTTTCTAAAGAAGACCTGGAGCAGAAAAAAGGCCTTTTCTTCCCTTCTTGACTCCCTCATCTCATTTTGGTTTCTTTGTCTCTGCACATCTGATTCTCCCCCTTTTTTTTTCTTCTCTTCTTCTGCTGGAGGATAAGCTCACCAAGCCAACCAGCAAAAATGTGGTGCCTCTCAATTTTTCCAAACTACTATTCCAAGTGATTTGAGAAATGATCTACTACGATACTCCTCAAGAGTCAAATGTTGACCTCGGGGAGCCTTTTTTGGTATTGCTCCATGATCAGAGCTTTACGAGCTAGTGTTTTTTCTGCATATCAGCCCAAACTCTCAATGATAATTTTACTGGAGGCTGATTTTTGTAAAATTTTGTGCCATAAAAGCCTTGTTGGCTTGTCTCTTGCTTGGCTTTAGATCATTCTCAAGCCATTTTTTTCCTGCTGTTGCTCTGACACAGGTTGTTTTTGCTGGTCTTGTTGGTGCCTGATCCACTGCTATCCTTTTCACACCTCTTTTTTTTTTTTCTTCATCCTGCACAAGTTTCTGCTGCCTGTTAGTCGGCATCACCGGTTTTGGGACCAAAACCACATCATGTGGTCTGTAACAGTATGCACAACCATGCCGTGAGGACCAAATTTGTTTTATTATTGTTATTATTATTAAAAGCCTTTGCTTCCATTCGGAGTTTGTTTTTTTGAGTAATATATGTATTCATTGTTTGGGTCGAATCCCCTTGCTTTTTTAACACAAATGTTTTGCAAACCACTATTTGAAATGGTGCACTGTTATGGGCTTATGGTGAGCAGATGAGGCCAAGTCATGGTTTCTTCATTATAATTTTCTTTTCATTTGCTTTAAAGAGCCATATTCTACCCAGGGAAGAAAGGTTGAAGTTGTTTTGTTTTTTTACCGTGAGTTCAAAGCAGTGGCACTGCCAGATTTAAAAGGTTCAAAAGCCGTGCAGATCTAAAATATGTATTATGAACACAGTAATGGGAGCGAATTGTAACACTTAATAGTATACAAATTTAAGAAACAGGGGTGAACACATAGTTTTAACTGGAAAAAGCCCACAATGATGTGTAATCACTTTGTTACTGTCTGTATCTTGTGTAATGATACCTAAATTCTTTTTTTAAATAAAAACCATGATTTTTACTGTCACTGAAAAAAAAAAAAAAAAAAA")
    # x.add_source("desm", "CATTTACACAGCGTACAAACCCAACAGGCCCAGTCATGAGCACGAAATATTCAGCCTCCGCCGAGTCGGCGTCCTCTTACCGCCGCACCTTTGGCTCAGGTTTGGGCTCCTCTATTTTCGCCGGCCACGGTTCCTCAGGTTCCTCTGGCTCCTCAAGACTGACCTCCAGAGTTTACGAGGTGACCAAGAGCTCCGCTTCTCCCCATTTTTCCAGCCACCGTGCGTCCGGCTCTTTCGGAGGTGGCTCGGTGGTCCGTTCCTACGCTGGCCTTGGTGAGAAGCTGGATTTCAATCTGGCTGATGCCATAAACCAGGACTTCCTCAACACGCGTACTAATGAGAAGGCCGAGCTCCAGCACCTCAATGACCGCTTCGCCAGCTACATCGAGAAGGTGCGCTTCCTCGAGCAGCAGAACTCTGCCCTGACGGTGGAGATTGAGCGTCTGCGGGGTCGCGAGCCCACCCGTATTGCAGAGCTGTACGAGGAGGAGATGAGAGAGCTGCGCGGACAGGTGGAGGCACTGACCAATCAGAGATCCCGTGTGGAGATCGAGAGGGACAACCTAGTCGATGACCTACAGAAACTAAAGCTCAGACTTCAAGAGGAGATCCACCAGAAAGAGGAAGCTGAAAACAACCTTTCTGCTTTCAGAGCTGATGTCGATGCTGCCACTCTGGCCAGGCTGGACCTGGAAAGACGTATCGAGGGTCTTCACGAAGAGATTGCATTCCTCAGGAAGATTCATGAGGAGGAGATCCGTGAGCTGCAGAACCAGATGCAGGAGAGTCAGGTGCAGATCCAAATGGACATGTCCAAACCAGACCTGACTGCGGCCCTCAGAGACATTCGCCTGCAGTACGAGGCTATCGCTGCCAAGAATATCAGCGAGGCCGAGGACTGGTATAAGTCTAAGGTTTCAGATTTGAACCAGGCAGTGAACAAGAATAACGAGGCTCTCAGAGAAGCCAAGCAGGAGACCATGCAGTTCCGTCACCAGCTCCAGTCCTACACCTGCGAGATTGACTCTCTCAAGGGCACCAATGAGTCTCTGAGGAGGCAAATGAGTGAGATGGAGGAGCGGCTGGGACGTGAGGCCGGTGGTTATCAGGACACTATCGCCCGTCTCGAGGCTGAGATCGCAAAAATGAAAGACGAGATGGCCCGCCACCTCCGCGAGTACCAGGATCTGCTGAATGTGAAGATGGCTCTGGATGTGGAGATCGCCACCTACAGGAAGCTTTTGGAAGGAGAGGAGAGCAGGATCTCGCTGCCCGTGCAGTCCTTTTCATCCCTGAGTTTCAGAGAGAGCAGTCCAGAGCAGCACCACCACCAGCAGCAGCAACCACAACGCTCATCTGAAGTCCACTCCAAGAAAACAGTCCTGATCAAGACCATCGAGACCCGCGATGGCGAGGTCGTCAGCGAGTCCACACAGCACCAGCAGGACGTCATGTAAAGCTTGAGAAACAGATCGAGTTTCACAGAATGCCTTGCATTTTCACTGATGGCCTCAGGCTTTTTTAAGCACACACCCAGTATTGCCGTGACCCATTACCGCATGTGGATGACGCATGGAGACAAAAGGAAAGTGAGCTGAAAAACCAGAGGGAGGAAAAGTGGAATGGTGTGATGCTGAGCGTTCAGAAAGTGGCCAGATGAGCTCAGAGTTTCTGATTTAATGAATGTATGTGTGCGTGTGTGTGTGGTTGGGTCATATCTGAGACACTGTTCCACAGCAACAAAAACAATAAAATTCACTGTATTTTCTCCTAAAAAAAAAAAAAAAAAAAAAAAAAA")
    # x.add_window_constraint(('b*', 'a*'), "tmp3")
    # x.add_window_constraint(['c*', 'b*'], "tmp3")
    # x.add_pattern_constraints([x * 4 for x in list("ACTG")] + [x * 6 for x in list("WSMKYR")])
    # x.add_similarity_constraint("h1", "S"*72, (0.45, 0.55))

    # x.add_similarity_constraint("i1", "H" * 36, (0.8, 1))
    # x.add_similarity_constraint("i2", "H" * 36, (0.8, 1))
    # x.add_similarity_constraint("h2", "H" * 72, (0.8, 1))
    # x.add_similarity_constraint("h1", "H" * 72, (0.8, 1))

    x.parameters.f_stop = 0.05
    x.log = "design_test.log"

    y = x()
    print(y)

    # print(x.to_json())
    x = design.Specification(model=Model(), json=x.to_json())
    # print(x.to_json())
    print(len(x.complexes))
    for tube in x.tubes:
        print(tube.name, len(tube.targets))

    y = x()
    print(y)
    # print(y.to_json())
    # y = design.RawResult(json=y.to_json())
    # print(y.to_json())
    # print([+i.normalized_defect/sum(+j.normalized_defect_contribution for j in i.complexes) for i in y.tubes])
    # defs = [i.normalized_defect for i in y.tubes]
    # print(sum(defs)/len(defs))
    # y = x.designer()
    # y.run().print_results()

def test_evaluate_wobble():
    my_model = Model(material='rna')

    seq1 = 'GGGGGAAACCCCC'
    seq1_mfe = mfe(seq1, model=my_model)
    st1 = str(seq1_mfe[0].structure) # (((((...)))))
    ensemble_defect = defect(strands=seq1, structure=st1, model=my_model)
    print(ensemble_defect)

    seq2 = 'GGGGGAAACCUCC'
    seq2_mfe = mfe(seq2, model=my_model)
    st2 = str(seq2_mfe[0].structure) # (((((...)))))
    ensemble_defect = defect(strands=seq2, structure=st2, model=my_model)
    print(ensemble_defect)
