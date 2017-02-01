import math
import collections.abc
from decimal import Decimal, getcontext
from Bio.SeqRecord import SeqRecord


class LazyDict(collections.MutableMapping):
    """
    A dictionary that can calculate its values lazily.
    This class behaves like a regular dictionary with the only exception:
    its values can be recomputed when they are asked for. This is controlled by
    three parameters:
    `self.recompute` is a function that is used whenever values need to be found
    It should be a one-argument function that accepts a key and returns a value.
    It is the duty of user to make sure all object(s) referenced in this
    function do still exist when self.recompute is called. Most likely use case
    is a lambda that refers to some external object(s).
    `self.keys_source` is a function that is used to generate key list. It should
    be a zero-argument function that returns an iterable with all keys LazyDict
    is supposed to have right now. As above, it's user's duty to make sure all
    referenced objects exist whenever it's called. Most likely use case is
    something like lambda: source_dict.keys()
    `self.values_fresh` is a boolean attribute that is True if values currently
    in dict are useful and to False if they are not. Initially set to False and
    changed to True on every recalculation of values, this attribute needs to be
    manually changed to False by caller whenever there is a reason to think that
    values are getting out of date.

    Items are refreshed whenever any LazyDict element is accessed alone (via
    `lazydict['key']`) or a dictionary as a whole is used (via lazydict.keys(),
    lazydict.values() or lazydict.items()). It happens iff by that moment
    lazydict.values_fresh is set to False.
    """
    def __init__(self, keys_source=None, recompute=None, *args, **kwargs):
        super(LazyDict, self).__init__(*args, **kwargs)
        if not hasattr(recompute, '__call__'):
            raise ValueError('Only callables accepted as recompute')
        self.recompute = recompute
        if not hasattr(keys_source, '__call__'):
            raise ValueError('Only callables accepted as keys_source')
        self.keys_source = keys_source
        self.values_fresh = False
        self._dict = {}

    def __getitem__(self, item):
        if self.values_fresh:
            return self._dict[item]
        else:
            self._update_values()
            return self._dict[item]

    def __setitem__(self, key, value):
        self._dict[key] = value

    def __delitem__(self, key):
        self._dict.__delitem__(key)

    def __iter__(self):
        return self._dict.__iter__()

    def __len__(self):
        return len(self._dict)

    def keys(self):
        if not self.values_fresh:
            self._update_values()
        return self._dict.keys()

    def values(self):
        if not self.values_fresh:
            self._update_values()
        return self._dict.values()

    def items(self):
        if not self.values_fresh:
            self._update_values()

    def __eq__(self, other):
        if not self.values_fresh:
            self._update_values()
        if not isinstance(other, (dict, LazyDict)):
            raise ValueError(
                'LazyDict can be compared only with other LazyDict or dict')
        if sorted(list(self.keys())) == sorted(list(other.keys())):
            for x in self.keys():
                if other[x] != self[x]:
                    # Different value
                    return False
        else:
            # Different keys
            return False
        return True

    def _update_values(self):
        """Build new values. Not meant to be called explicitly"""
        k = self.keys_source()
        self._dict = {x: self.recompute(x) for x in k}
        self.values_fresh = True
        # for x in self._dict.keys():
        #     self._dict[x] = self.recompute(x)


class Composition(collections.abc.MutableMapping):
    '''
    A class for aminoacid k-mer composition of a sequence or a sequence set.
    May be initialized with a single sequence, sequence iterable, file-like object
    or nothing (in latter case use Composition.process(SeqRecord)).
    '''
    def __init__(self, k, seq=None, fh=None):
        """
        Create a new composition
        :param k: int. Length of k-mer. This parameter is obligatory and cannot
        be changed for an exisiting model.
        :param seq: SeqRecord. If not None, this sequence will be used to build
        a model.
        :param fh: a filehandle or file-like IO stream. If not None, model will
        be read from this file.
        :return:
        """
        self.k = k
        self.kmer_count = 0
        self.sequence_count = 0
        self.abs_distribution = {}
        self.log_distribution = LazyDict(keys_source=lambda: self.keys(),
              recompute=lambda x: math.log10(self.abs_distribution[x]) - math.log10(self.kmer_count))
        self.relative_distribution = LazyDict(keys_source=lambda: self.keys(),
              recompute=lambda x: self.abs_distribution[x]/self.kmer_count)
        if seq:
            self.process(seq)
        if fh:
            self.read(fh)

    def __getitem__(self, item):
        return self.log_distribution[item]

    def __setitem__(self, key, value):
        raise NotImplementedError('Manually adding values not supported')
        pass

    def keys(self):
        return self.abs_distribution.keys()

    def __delitem__(self, key):
        raise NotImplementedError('Element deletion not supported')
        pass

    def __iter__(self):
        return iter(self.log_distribution)

    def __len__(self):
        pass

    def process_single_sequence(self, sequence):
        """Add statistics of a single sequence to this Composition"""
        if not isinstance(sequence, SeqRecord):
            raise ValueError('Only SeqRecord objects can be added to Composition')
        s = str(sequence.seq)
        for j in range(len(s) - self.k+1):
            self.kmer_count += 1
            try:
                self.abs_distribution[s[j:j + self.k]] += 1
            except KeyError:
                self.abs_distribution.update({s[j:j + self.k]: 1})
        self.sequence_count += 1
        self.log_distribution.values_fresh = False
        self.relative_distribution.values_fresh = False

    def process(self, item):
        '''
        Add a sequence statistic to this composition object.
        :param item: SeqRecord or an iterable of SeqRecords
        :param update: Boolean
        :return:
        '''
        if isinstance(item, SeqRecord):
            self.process_single_sequence(item)
        elif isinstance(item, list) or isinstance(item, tuple):
            for single_item in item:
                self.process_single_sequence(single_item)

    def distribution(self, sequence):
        '''
        Count distribution of k-mers in a given sequence and return it as a dictionary of
        {'AAA': count('AAA'), 'AAB': count('AAB'), etc.}
        :param seq: Bio.SeqRecord
        :return:
        '''
        s = str(sequence.seq)
        d = {}
        count = 0
        for j in range(len(s) - self.k):
            count += 1
            try:
                d[s[j:j+3]] += 1
            except KeyError:
                d.update({s[j:j+3]: 1})
        return d

    def prob(self, sequence):
        """
        Takes seq and calculates probability of it being generated by this Composition
        :param seq: Bio.SeqRecord
        :return:
        """
        pseudocount = math.log10(1/self.kmer_count)
        p = 0
        s = str(sequence.seq)
        for j in range(len(s) - self.k):
            try:
                pr = self.log_distribution[s[j:j+self.k]]
            except KeyError:
                pr = pseudocount
            p += pr
        return p

    def write(self, fh):
        '''
        Write an absolute k-mer distribution to a given filehandle.
        '''
        pseudocount = math.log10(1 / self.kmer_count)
        fh.write('Pseudo\t{0}\n'.format(self.pseudocount))
        fh.write('SCount\t{0}\n'.format(self.sequence_count))
        for j in self.keys():
            fh.write('{0}\t{1}\n'.format(j, self[j]))

    def read(self, fh):
        """
        Read a model from a filehandle.
        Throws an error if a filehandle contains a model declared for a different k
        :param fh:
        :return:
        """
        getcontext().prec = 100
        model = {}
        checked_len = False
        for line in fh:
            kmer, fl = line.split(sep='\t')
            if kmer == 'Pseudo':
                self.pseudocount = float(fl)
                continue
            elif kmer == 'SCount':
                self.sequence_count = int(fl)
            elif not checked_len:
                checked_len = True
                if len(kmer) != self.k:
                    raise ValueError(
                        'Model and file are using a different kmer length!')
            model[kmer] = fl


def euclidean(a, b):
    '''
    Count euclidean distances between two kmer distribution dicts
    :param a: dict
    :param b: dict
    :return: int
    '''
    getcontext().prec = 500
    n = Decimal(0)
    for j in set(a.keys()).intersection(set(b.keys())):
        aj = a.relative_distribution[j]
        bj = b.relative_distribution[j]
        n += Decimal((aj - bj)**2)
    return n.sqrt()


def ffp_distance(a, b):
    '''
    Calculate Jensen-Shannon divergence.
    Sims, Se-Ran Jun, Wu, Kim 2008 "Alignment-free genome comparison with
    feature frequency profiles and optimal resolutions" have shown it to be a
    viable distance estimator (although it's still not a metric), hence the name
    :param a: Composition
    :param b: Composition
    :return: float
    '''
    if not (a.k == b.k):
        raise ValueError('Comparison possible only for models with the same k-mer length')
    avg_model = {}
    # Not sure if these midpoint calculations are strictly correct, but at the
    # very least they seem to produce reasonable result
    for j in set(a.keys()).intersection(set(b.keys())):
        avg_model[j] = (a.relative_distribution[j]+b.relative_distribution[j])/2
    return kullback_leibler(a.relative_distribution, avg_model)/2 +\
           kullback_leibler(b.relative_distribution, avg_model)/2


def kullback_leibler(d, e):
    '''
    Calculate kullback-Leibler divergence on a pair of compositions.
    :param d: dict
    :param e: dict
    :return:
    '''
    kl = Decimal(0.0)
    for j in set(d.keys()).intersection(set(e.keys())):
        kl_local = Decimal(d[j]) * (Decimal(math.log2(d[j])) - Decimal(math.log2(e[j])))
        kl += kl_local
    return float(kl)
