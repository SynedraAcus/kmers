#  Tests for the LazyDict

from kmers.kmers import LazyDict


def test_lazy_dict_recompute_on_use():
    # Key is an index in source
    # value is twice the value in it
    source = [1, 2, 3]
    ld = LazyDict(recompute=lambda x: source[x]*2,
                  keys_source=lambda: range(len(source)))
    expect = {x: 2*source[x] for x in range(len(source))}
    # Use value to recompute
    a = ld[1]
    assert ld == expect
    source.append(4)
    source.pop(0)
    ld.values_fresh = False
    # Recompute even if no value was directly used
    expect = {0: 4, 1: 6, 2: 8}
    assert ld == expect


def test_recompute_on_indirect_use():
    # Recompute if __eq__, keys() or values() is called
    source = [1, 2, 3]
    ld = LazyDict(recompute=lambda x: source[x]*2,
                  keys_source=lambda: range(len(source)))
    assert sorted(list(ld.keys())) == [0, 1, 2]
    source.append(1)
    ld.values_fresh = False
    assert sorted(list(ld.values())) == [2, 2, 4, 6]
    source.append(1)
    expect = {0: 2, 1: 4, 2: 6, 3: 2, 4: 2}
    ld.values_fresh = False
    assert ld == expect

#  Tests for Composition and distances
from kmers.kmers import Composition, euclidean, ffp_distance
from Bio.SeqRecord import SeqRecord


def test_relative():
    record = SeqRecord(seq='ACACAT')
    compo = Composition(k=3, seq=record)
    assert compo.relative_distribution == {'ACA': 0.5, 'CAC': 0.25, 'CAT': 0.25}


def test_log():
    record = SeqRecord(seq='AAAAA')
    compo = Composition(k=3, seq=record)
    # No more complex tests b/c rounding logarithms is a pain in the ass
    assert compo.log_distribution == {'AAA': 0}


def test_similars():
    # Zero distance between sequences of similar k-mer composition
    a = SeqRecord(seq='ACAC')
    b = SeqRecord(seq='CACA')
    assert euclidean(Composition(seq=a, k=3), Composition(seq=b, k=3)) == 0
    assert ffp_distance(Composition(seq=a, k=3), Composition(seq=b, k=3)) == 0

def test_distances():
    a = Composition(3, SeqRecord(seq='ACTGACTG'))
    b = Composition(3, SeqRecord(seq='CTGACTGA'))
    #  Let's forget small rounding errors for a while
    assert 0.235 < euclidean(a, b) < 0.236
    assert 0.040 < ffp_distance(a,b) < 0.041
