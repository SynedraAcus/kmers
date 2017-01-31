import math
import collections.abc
from decimal import *
from Bio.SeqRecord import SeqRecord


class LazyDict(dict):
    """
    A dictionary that can calculate its values lazily
    This class behaves like a regular dictionary with the only exception:
    its values can be recomputed when they are asked for. This is controlled by
    two parameters:
    `self.recompute` is a function that is used whenever values need to be found
    It should be a one-argument function that accepts a key and returns a value.
    It is the duty of user to make sure all object(s) referenced in this
    function do still exist when self.recompute is called. Most likely use case
    is a lambda that refers to some external object(s).
    `self.values_fresh` is a boolean attribute that is True if values currently
    in dict are useful and to False if they are not. Initially set to False and
    changed to True on every recalculation of values, this attribute needs to be
    manually changed to False by caller whenever there is a reason to think that
    values are getting out of date.

    """
    def __init__(self, recompute=None, *args, **kwargs):
        super(self, LazyDict).__init__(*args, **kwargs)
        if not hasattr(recompute, '__call__'):
            raise ValueError('Only callables accepted as recompute')
        self.recompute = recompute
        self.values_fresh = False
        self._dict = {}

    def __getitem__(self, item):
        if self.values_fresh:
            return self._dict[item]
        else:
            self.update_values()
            return self._dict[item]

    def __setitem__(self, key, value):
        self._dict[key] = value

    def _update_values(self):
        """Build new values.
        Old ones are discarded. Not meant to be called explicitly"""
        for x in self._dict.keys():
            self._dict[x] = self.recompute(x)



class Composition(collections.abc.MutableMapping):
    '''
    A class for aminoacid k-mer composition of a sequence or a sequence set.
    May be initialized with a single sequence, sequence iterable, file-like object
    or nothing (in latter case use Composition.process(SeqRecord)).
    '''
    def __init__(self, k, seq=None, fh=None):
        """
        :param seq: SeqRecord
        :param k: int
        :return:
        """
        self.k = k
        self.kmer_count = 0
        self.sequence_count = 0
        self.abs_distribution = {}
        self.log_distribution = {}
        if seq:
            self.process(seq)
        if fh:
            self.read(fh)
        self.pseudocount = None

    def __getitem__(self, item):
        return self.log_distribution[item]

    def __setitem__(self, key, value):
        raise NotImplementedError('Manually adding values not supported')
        pass

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
        for j in range(len(s) - self.k):
            self.kmer_count += 1
            try:
                self.abs_distribution[s[j:j + self.k]] += 1
            except KeyError:
                self.abs_distribution.update({s[j:j + self.k]: 1})
        self.sequence_count += 1

    def process(self, item, update=True):
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
                # Even if we accidentrally start iterating over SeqRecord-like
                # object (say, Seq), it returns `str` for every letter, thus
                # raising exception in process_single_sequence()
                self.process_single_sequence(single_item)
        if update:
            self.update_relative()

    def update_relative(self):
        '''
        Update self.relative_distribution to be in accord with abs distribution
        :return:
        '''
        getcontext().prec = 100
        for j in self.abs_distribution.keys():
            self.log_distribution[j] = math.log10(self.abs_distribution[j]) - math.log10(self.kmer_count)
        #  Define a 'pseudocount' value: for k-mers that were not present in a traning set the log
        #  probability should be the same as for k-mers that were present exactly once
        self.pseudocount = 0 - math.log10(self.kmer_count)

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
        p = 0
        s = str(sequence.seq)
        for j in range(len(s) - self.k):
            try:
                pr = self.log_distribution[s[j:j+self.k]]
            except KeyError:
                pr = self.pseudocount
            p += pr
        return p

    def write(self, fh):
        '''
        Write a relative k-mer usage to a given filehandle
        '''
        fh.write('Pseudo\t{0}\n'.format(self.pseudocount))
        fh.write('SCount\t{0}\n'.format(self.sequence_count))
        for j in self.keys():
            fh.write('{0}\t{1}\n'.format(j, self[j]))

    def read(self, fh):
        """
        Read a model from a filehandle. Throw an error if a filehandle contains a model declared for a different k
        :param fh:
        :return:
        """
        getcontext().prec = 100
        for line in fh:
            kmer, fl = line.split(sep='\t')
            if kmer == 'Pseudo':
                self.pseudocount = float(fl)
                continue
            if kmer == 'SCount':
                self.sequence_count = int(fl)
            self.log_distribution.update({kmer: float(fl)})


def euclidean(a, b):
    '''
    Count euclidean distances between two kmer distribution dicts
    :param a: dict
    :param b: dict
    :return: int
    '''
    getcontext().prec = 500
    n = Decimal(0)
    for j in set(a.keys()).union(set(b.keys())):
        aj = a.get(j, 0)
        bj = b.get(j, 0)
        n += (aj - bj)**2
    return n.sqrt()


def ffp_distance(a, b):
    '''
    Calculate "Feature frequency profile" distance, as per Sims, Se-Ran Jun, Wu, Kim 2008 and protein variants thereof
    Important note: k-mer length is fixed, unlike in described papers, where they derive it for every dataset
    :param a: Composition
    :param b: Composition
    :return: float
    '''
    if not (a.k == b.k):
        raise ValueError('Comparison possible only for models with the same k-mer length')
    avg_model = Composition(a.k)
    for j in set(a.keys()).intersection(set(b.keys())):
        avg_model.relative_distribution[j] = (a.relative_distribution[j]+b.relative_distribution[j])/2
    return kullback_leibler(a, avg_model)/2 + kullback_leibler(b, avg_model)/2


def kullback_leibler(d, e):
    '''
    Calculate kullback-Leibler divergence on a pair of compositions.
    :param d:
    :param e:
    :return:
    '''
    kl = Decimal(0.0)
    for j in set(d.keys()).union(set(e.keys())):
        try:
            kl_local = d.relative_distribution[j]*Decimal(math.log2(d.relative_distribution[j]) - math.log2(e.relative_distribution[j]))
            kl += kl_local
        except KeyError:
            pass
    return float(kl)
