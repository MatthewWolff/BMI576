NEGATIVE_INFINITY = float("-inf")


# other cells in top row and left most column = negative infinity
def align_global_affine_cs_gaps(x, y, submatrix, gap_score_matrix, space_score):
    _x, _y = len(x) + 1, len(y) + 1  # len + 1, reg len
    M = [[float()] * _y for _ in range(_x)]
    I_x = [[float()] * _y for _ in range(_x)]
    I_y = [[float()] * _y for _ in range(_x)]
    traceback = dict()

    # initialize first cells
    I_x[0][0] = gap_score_matrix["", (x + " ")[0].strip()]  # add space + strip after indexing to avoid bad index
    I_y[0][0] = gap_score_matrix["", (y + " ")[0].strip()]

    # initialize the rest
    for i in range(len(y)):
        M[0][i + 1] = I_x[0][i + 1] = NEGATIVE_INFINITY
        I_y[0][i + 1] = I_y[0][i] + space_score
        traceback[("I_y", 0, i + 1)] = ("I_y", 0, i)
    for i in range(len(x)):
        M[i + 1][0] = I_y[i + 1][0] = NEGATIVE_INFINITY
        I_x[i + 1][0] = I_x[i][0] + space_score
        traceback[("I_x", i + 1, 0)] = ("I_x", i, 0)

    # compute 3 values and store connections
    for i in range(len(x)):
        for j in range(len(y)):
            # M matrix
            entries = {"M": M[i][j], "I_x": I_x[i][j], "I_y": I_y[i][j]}
            options = [(k, i, j) for k, v in entries.items() if v == max(entries.values())]  # select best
            traceback["M", i + 1, j + 1] = options[0] if len(options) is 1 else options  # tuple or list of tuples
            M[i + 1][j + 1] = max(entries.values()) + submatrix[(x[i], y[j])]  # take best

            # I_x matrix
            gap = gap_score_matrix[y[j], "" if len(y) == j + 1 else y[j + 1]]
            entries = {"M": M[i][j + 1] + space_score + gap, "I_x": I_x[i][j + 1] + space_score}
            options = [(k, i, j + 1) for k, v in entries.items() if v == max(entries.values())]
            traceback["I_x", i + 1, j + 1] = options[0] if len(options) is 1 else options  # tuple or list of tuples
            I_x[i + 1][j + 1] = max(entries.values())  # take best

            # I_y matrix
            gap = gap_score_matrix[(x[i], "" if len(x) == i + 1 else x[i + 1])]
            entries = {"M": M[i + 1][j] + space_score + gap, "I_y": I_y[i + 1][j] + space_score}
            options = [(k, i + 1, j) for k, v in entries.items() if v == max(entries.values())]
            traceback["I_y", i + 1, j + 1] = options[0] if len(options) is 1 else options  # tuple or list of tuples
            I_y[i + 1][j + 1] = max(entries.values())  # take best

    entries = [I_x[-1][-1], M[-1][-1], I_y[-1][-1]]
    matrices = ["I_x", "M", "I_y"]

    # traceback
    x_str = ""
    y_str = ""
    matrix, curr_x, curr_y = matrices[entries.index(max(entries))], len(x), len(y)
    while not curr_x == curr_y == 0:
        options = traceback[matrix, curr_x, curr_y]
        if type(options) is list:  # check for multiple - sort and check for I_x
            options.sort()
            options = options[0] if options[0][0] == "I_x" else options[len(options) - 1]
        x_str += x[options[1]] if (curr_x - 1) == options[1] else "-"
        y_str += y[options[2]] if (curr_y - 1) == options[2] else "-"
        matrix, curr_x, curr_y = options
    best_alignment = [x_str[::-1], y_str[::-1]]
    return max(entries), best_alignment


DNA = list("ACGT")
DNA_AND_EMPTY = DNA + [""]

# A simple match=+1 and mismatch=-1 matrix
basic_submatrix = {(a, b): 1 if a == b else -1 for a in DNA for b in DNA}

# A gap score matrix with the same penalty for all contexts
uniform_gap_score_matrix = {(a, b): -2 for a in DNA_AND_EMPTY for b in DNA_AND_EMPTY}

# A gap score matrix with smaller penalties for gaps flanked by A and/or T
biased_gap_score_matrix = {(a, b): -1 if 'T' in (a, b) or 'A' in (a, b) else -2
                           for a in DNA_AND_EMPTY for b in DNA_AND_EMPTY}

# A gap score matrix with all zeros.  Equivalent to a linear gap penalty function
zero_gap_score_matrix = {(a, b): 0 for a in DNA_AND_EMPTY for b in DNA_AND_EMPTY}


# A utility function for displaying a substitution or context-sensitive gap score matrix


def print_dict_matrix(m, width=5):
    """Prints the dict-based matrix m with the specified width for each column."""

    def print_row(fields):
        print("".join(["{:>{width}}".format(field, width=width) for field in fields]))

    labels = sorted({c for pair in m for c in pair})
    print_row([""] + labels)
    for a in labels:
        print_row([a] + [m.get((a, b), "") for b in labels])


def score_alignment(alignment, submatrix, gap_score_matrix, space_score):
    """Returns the score of the alignment given a substitution matrix, gap score matrix, and space score."""
    aligned_pair_scores = [submatrix[pair] for pair in zip(*alignment) if '-' not in pair]
    gap_x_scores, gap_y_scores = [gap_scores(s, gap_score_matrix, space_score) for s in alignment]
    return sum(aligned_pair_scores) + sum(gap_x_scores) + sum(gap_y_scores)


def gap_scores(s, gap_score_matrix, space_score):
    """Returns a list of the scores of the gaps within the aligned sequence s."""
    return [gap_score_matrix[context] + space_score * length for length, context in gaps(s)]


import re


def gaps(s):
    """Returns an iterator over the gaps of the string s.
    Each element of the iterator has the form (length, context)
    where length is the length of the gap and context is a tuple of the two characters flanking the gap."""
    gap_pattern = re.compile('-+')
    for match in gap_pattern.finditer(s):
        yield (match.end() - match.start(),
               (s[match.start() - 1: match.start()], s[match.end(): match.end() + 1]))


test_case_inputs = {
    'small_1': ("AGTA", "AGA", basic_submatrix, uniform_gap_score_matrix, -1),
    'small_2': ("AGT", "AGT", basic_submatrix, uniform_gap_score_matrix, -1),
    'small_3': ("AGT", "G", basic_submatrix, uniform_gap_score_matrix, -1),
    'small_4': ("G", "AGT", basic_submatrix, uniform_gap_score_matrix, -1),
    'small_5': ("AGT", "", basic_submatrix, uniform_gap_score_matrix, -1),
    'small_6': ("", "AGT", basic_submatrix, uniform_gap_score_matrix, -1),
    'small_7': ("A", "", basic_submatrix, uniform_gap_score_matrix, -1),
    'small_8': ("", "", basic_submatrix, uniform_gap_score_matrix, -1),
    'small_9': ("CTC", "TGCT", basic_submatrix, zero_gap_score_matrix, -1),
    'small_10': ("TATA", "TTAC", basic_submatrix, biased_gap_score_matrix, -1),
    'large_1': ("TCATTCTGTTTATACTATCTTACTGGTTACCTTAATAATACAATCAGAATCGTAATTCGTCCTGTTCGT",
                "TGTATAACTATGTCATCTAACCCCAAGCTTATCACTGCTTACGGAGGACAG",
                basic_submatrix, biased_gap_score_matrix, -1)
}

test_case_correct_outputs = {
    'small_1': (0, ['AGTA',
                    'AG-A']),
    'small_2': (3, ['AGT',
                    'AGT']),
    'small_3': (-5, ['AGT',
                     'G--']),
    'small_4': (-5, ['--G',
                     'AGT']),
    'small_5': (-5, ['AGT',
                     '---']),
    'small_6': (-5, ['---',
                     'AGT']),
    'small_7': (-3, ['A',
                     '-']),
    'small_8': (0, ['',
                    '']),
    'small_9': (-1, ['--CTC',
                     'TGCT-']),
    'small_10': (-1, ['TATA-',
                      'T-TAC']),
    'large_1': (-15, ['TCATTCTGTTTATACTATCTTACTGGTTACCTTAATAATACAATCAGAATCGTAATTCGTCCTGTTCGT',
                      'T------GTATA-ACTATGTCA-TC-TAACCCCAAGCTTA---TCACT-GCTTA---CGGA--GGACAG'])
}

import numbers


def check_valid_alignment_result(result, x, y):
    """Checks that the alignment result is valid for sequences x and y."""
    assert isinstance(result, tuple), "Output is not a tuple"
    assert len(result) == 2, "Output does not have exactly two elements"
    score, alignment = result
    assert isinstance(alignment, list), "Alignment is not a list"
    assert isinstance(score, numbers.Number), "Score is not a number"
    assert len(alignment) == 2, "Alignment does not have exactly two elements"
    assert all(isinstance(element, str) for element in alignment), "Alignment elements are not strings"
    assert len(alignment[0]) == len(alignment[1]), "Alignment strings do not have the same length"
    assert alignment[0].replace('-', '') == x, "First string of alignment is not x"
    assert alignment[1].replace('-', '') == y, "Second string of alignment is not y"


def check_valid_alignment_score(result, submatrix, gap_score_matrix, space_score):
    """Checks that the computed score of the alignment is equal to the score given in the result."""
    score, alignment = result
    computed_score = score_alignment(alignment, submatrix, gap_score_matrix, space_score)
    assert computed_score == score, "Computed score ({}) does not equal the returned score ({})".format(computed_score,
                                                                                                        score)


def check_test_case(case_name, test_name=None, valid_result=True, valid_score=True, correct_alignment=True):
    inputs = test_case_inputs[case_name]
    correct_output = test_case_correct_outputs[case_name]
    result = align_global_affine_cs_gaps(*inputs)
    if valid_result:
        check_valid_alignment_result(result, *inputs[:2])
    if valid_score:
        check_valid_alignment_score(result, *inputs[2:])
    if correct_alignment:
        assert result == correct_output
    print("SUCCESS:", test_name if test_name else case_name, "passed!")


def test_all():
    list(map(check_test_case, list(map(lambda x: x[0] + str(x[1]), zip(["small_"] * 10, range(1, 11)))) + ["large_1"]))


check_test_case("small_1", test_name="small_1 valid output", valid_score=False, correct_alignment=False)
test_all()
