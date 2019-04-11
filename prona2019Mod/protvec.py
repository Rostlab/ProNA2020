import prona2019Mod.utils as utils
import itertools as it
from six import iteritems, string_types, PY2, next
import numpy as np
import sys 

def _is_single(obj):
    """
    Check whether `obj` is a single document or an entire corpus.
    Returns (is_single, new) 2-tuple, where `new` yields the same
    sequence as `obj`.

    `obj` is a single document if it is an iterable of strings.  It
    is a corpus if it is an iterable of documents.
    """
    obj_iter = iter(obj)
    temp_iter = obj_iter
    try:
        peek = next(obj_iter)
        obj_iter = it.chain([peek], obj_iter)
    except StopIteration:
        # An empty object is a single document
        return True, obj
    if isinstance(peek, string_types):
        # It's a document, return the iterator
        return True, obj_iter
    if temp_iter == obj:
        # Checking for iterator to the object
        return False, obj_iter
    else:
        # If the first item isn't a string, assume obj is a corpus
        return False, obj
'''
def _apply(corpus, chunksize=None, **kwargs):
    """Apply the transformation to a whole corpus and get the result as another corpus.

    Parameters
        ----------
    corpus : iterable of list of (int, number)
        Corpus in BoW format.
    chunksize : int, optional
        If provided - more effective processing (by group of documents) will performed.
    kwargs
        Arbitrary keyword arguments.

    Returns
    -------
    :class:`~gensim.interfaces.TransformedCorpus`
        Transformed corpus.

    """
    return TransformedCorpus(self, corpus, chunksize, **kwargs)
'''

def score_item(worda, wordb, components, scorer, phrasegrams):
    """score is retained from original dataset
    """
    try:
        return phrasegrams[tuple(components)][1]
    except KeyError:
        return -1

def analyze_sentence(sentence, threshold, common_terms, scorer,phrasegrams):
    """Analyze a sentence

    `sentence` a token list representing the sentence to be analyzed.

    `threshold` the minimum score for a bigram to be taken into account

    `common_terms` the list of common terms, they have a special treatment

    `scorer` the scorer function, as given to Phrases
    """
    s = [utils.any2utf8(w) for w in sentence]
    last_uncommon = None
    in_between = []
    # adding None is a trick that helps getting an automatic happy ending
    # has it won't be a common_word, nor score
    for word in s + [None]:
        is_common = word in common_terms
        if not is_common and last_uncommon:
            chain = [last_uncommon] + in_between + [word]
            # test between last_uncommon
            score = score_item(
                worda=last_uncommon,
                wordb=word,
                components=chain,
                scorer=scorer,
                phrasegrams=phrasegrams
            )
            if score > threshold:
                yield (chain, score)
                last_uncommon = None
                in_between = []
            else:
                # release words individually
                for w in it.chain([last_uncommon], in_between):
                    yield (w, None)
                in_between = []
                last_uncommon = word
        elif not is_common:
            last_uncommon = word
        else:  # common term
            if last_uncommon:
                # wait for uncommon resolution
                in_between.append(word)
            else: 
                yield (word, None)




def get_phrase(sentence,phrase_model):
        is_single, sentence = _is_single(sentence)
        if not is_single:
            # if the input is an entire corpus (rather than a single sentence),
            # return an iterable stream.
            sys.exit("It is not a protein sequence")

        delimiter = phrase_model['delimiter']
        bigrams = analyze_sentence(
            sentence,
            threshold=phrase_model['threshold'],
            common_terms=phrase_model['common_terms'],
            scorer=None,
            phrasegrams=phrase_model['phrasegrams'])  # we will use our score_item function redefinition
        new_s = []
        for words, score in bigrams:
            if score is not None:
                words = delimiter.join(words)
            new_s.append(words)
        return [utils.to_unicode(w) for w in new_s]

def split_ngrams(seq, n):
    """
    'AGAMQSASM' => [['AGA', 'MQS', 'ASM'], ['GAM','QSA'], ['AMQ', 'SAS']]
    """
    all_ngrams=[]
    for x in range(n):
        all_ngrams.append(zip(*[iter(seq[x:])]*n))
    str_ngrams = []
    for ngrams in all_ngrams:
        x = []
        for ngram in ngrams:
            x.append("".join(ngram))
        str_ngrams.append(x)
    return str_ngrams

def to_vecs(seq,phrase_model,kmer,word2vec_index):
    """
    convert sequence to three n-length vectors
    e.g. 'AGAMQSASM' => [ array([  ... * 100 ], array([  ... * 100 ], array([  ... * 100 ] ]
    """
    ngram_patterns = split_ngrams(seq, kmer)

    protvecs = []
    for ngrams in ngram_patterns:
        ngram_vecs = []

        if phrase_model=='none':
            ngramss = ngrams
        else:
            ngramss=get_phrase(get_phrase(ngrams,phrase_model),phrase_model)
            
        for ngram in ngramss:
            try:
                ngram_vecs.append(np.array(word2vec_index[ngram]))
            except  KeyError:
                continue
        protvecs.append(sum(ngram_vecs))
    return protvecs
