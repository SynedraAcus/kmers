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
    source = [1,2,3]
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
