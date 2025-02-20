# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        if not seqA or not seqB:
            raise ValueError("You gotta give me two sequences, bruh!")
        
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        lenA, lenB = len(seqA), len(seqB)

        if lenA == 0 or lenB == 0:
            raise ValueError("One of your sequences is empty, bruh!")
        
        align_mat_shape = (lenA + 1, lenB + 1)
        self._align_matrix = np.zeros(align_mat_shape)
        self._back = np.empty(align_mat_shape, dtype=object)
        self._gaps = np.zeros(align_mat_shape, dtype=int)
        
        # first row, first column, boooorrrriiinnggg
        for i in range(1, lenA + 1):
            self._align_matrix[i, 0] = self.gap_open + i * self.gap_extend
            self._gaps[i, 0] = 1 # gaparooni
            self._back[i, 0] = None 

        for j in range(1, lenB + 1):
            self._align_matrix[0, j] = self.gap_open + j * self.gap_extend
            self._gaps[0, j] = 1 # gaparooni
            self._back[0, j] = None

        # Implement global alignment here
        for i in range(1, lenA + 1):
            for j in range(1, lenB + 1):
                match_score = self.sub_dict.get((seqA[i - 1], seqB[j - 1]))  # Substitution score

                # Calculate scores for the three possible directions
                move_diagonal = self._align_matrix[i - 1, j - 1] + match_score

                move_up = self._align_matrix[i - 1, j] + (
                    self.gap_open + self.gap_extend 
                    if self._gaps[i - 1, j] == 0 else self.gap_extend
                )

                move_left = self._align_matrix[i, j - 1] + (
                    self.gap_open + self.gap_extend 
                    if self._gaps[i, j - 1] == 0 else self.gap_extend
                )

                # Choose the slayest (read: highest) score
                slayest_score = max(move_diagonal, move_up, move_left)

                # Update the backtrace matrix and direction
                if slayest_score == move_diagonal:
                    self._back[i, j] = (i - 1, j - 1)  # Diagonal
                    self._gaps[i, j] = 0  # No gap
                elif slayest_score == move_up:
                    self._back[i, j] = (i - 1, j)  # Up
                    self._gaps[i, j] = 1  # Gap
                elif slayest_score == move_left:
                    self._back[i, j] = (i, j - 1)  # Left
                    self._gaps[i, j] = 1  # Gap
                else:
                    raise RuntimeError("Something about my backtracing scores is off. Try again?")

                self._align_matrix[i, j] = slayest_score

        self.alignment_score = self._align_matrix[lenA, lenB]      		
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Start from the bottom right corner
        i, j = len(self._seqA), len(self._seqB)
        
        # Initialize empty alignment strings
        seqA_align = []
        seqB_align = []
        
        # Trace back until we hit the top left corner
        while i > 0 or j > 0:
            # Get the next position from the backtrace matrix
            next_pos = self._back[i, j]
            
            if next_pos is None:
                # We're at the edge of the matrix, handle remaining sequence
                while i > 0:
                    seqA_align.append(self._seqA[i-1])
                    seqB_align.append('-')
                    i -= 1
                while j > 0:
                    seqA_align.append('-')
                    seqB_align.append(self._seqB[j-1])
                    j -= 1
                break
                
            prev_i, prev_j = next_pos
            
            if i == prev_i:
                # Horizontal move (gap in seqA)
                seqA_align.append('-')
                seqB_align.append(self._seqB[j-1])
            elif j == prev_j:
                # Vertical move (gap in seqB)
                seqA_align.append(self._seqA[i-1])
                seqB_align.append('-')
            else:
                # Diagonal move (match/mismatch)
                seqA_align.append(self._seqA[i-1])
                seqB_align.append(self._seqB[j-1])
                
            i, j = prev_i, prev_j
        
        # Reverse the alignments (since we built them backwards)
        self.seqA_align = ''.join(reversed(seqA_align))
        self.seqB_align = ''.join(reversed(seqB_align))
        
        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
