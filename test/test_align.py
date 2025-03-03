# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    
    # Run alignment
    nw.align(seq1, seq2)
    
    # Expected matrices based on BLOSUM62 scoring
    # For sequence MYQR vs MQR:
    # M matches M (score 5)
    # Y gets gap (-11)
    # Q matches Q (score 5)
    # R matches R (score 5)
    
    expected_align = np.array([
        [  0, -11, -12, -13],
        [-11,   5,  -6,  -7],
        [-12,  -6,   4,  -7],
        [-13,  -7,  -1,   5],
        [-14,  -8,  -6,   4]
    ], dtype=np.float32)
    
    # Assert alignment matrix matches expected
    print(nw._align_matrix)
    np.testing.assert_array_almost_equal(nw._align_matrix, expected_align)
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    
    score, align1, align2 = nw.align(seq3, seq4)
    
    expected_score = 17
    expected_align1 = "MAVHQLIRRP"
    expected_align2 = "M---QLIRHP"
    
    assert score == expected_score
    assert align1 == expected_align1
    assert align2 == expected_align2
