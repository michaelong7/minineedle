from typing import Optional, Sequence, Any

from minineedle.needle import NeedlemanWunsch
from copy import deepcopy
from minineedle.typesvars import ItemToAlign

class Alignment_Data():
    def __init__(self, _alseq1, _alseq2, smatrix, _score, _identity, _nmatrix, _pmatrix, _gap_character, _used_indices, _seq1_start, _seq2_start) -> None:
        self._alseq1 = _alseq1
        self._alseq2 = _alseq2
        self.smatrix = smatrix
        self._score = _score
        self._identity = _identity
        self._nmatrix = _nmatrix
        self._pmatrix = _pmatrix
        self._gap_character = _gap_character
        self._used_indices = _used_indices
        self._seq1_start = _seq1_start
        self._seq2_start = _seq2_start


class SemiGlobal(NeedlemanWunsch[ItemToAlign]):
    """
    Semi-global algorithm that ignores gaps at the start and end of the first (reference) string.
    """

    def __init__(self, seq1: Sequence[ItemToAlign], seq2: Sequence[ItemToAlign]) -> None:
        super().__init__(seq1, seq2)
        self._used_indices = []
        self.alignments = {}
        self._seq1_start = 0
        self._seq2_start = 0

    def k_align(self, kbest: int) -> None:
        """
        Performs k semiglobal alignments with the given sequences and the
        corresponding ScoreMatrix.
        """
        if kbest not in self.alignments.keys():
            for i in range(kbest):
                self._align()
                self.alignments[i] = Alignment_Data(
                self._alseq1.copy(),
                self._alseq2.copy(),
                self.smatrix,
                self._score,
                self._identity,
                deepcopy(self._nmatrix),
                deepcopy(self._pmatrix),
                self._gap_character,
                self._used_indices.copy(),
                self._seq1_start,
                self._seq2_start)

    
    def _align(self) -> None:
        """
        Performs a semiglobal alignment with the given sequences and the
        corresponding ScoreMatrix.
        """
        self._add_initial_pointers()
        self._add_gap_penalties()
        self._fill_matrices()

        imax, jmax = self._get_last_cell_position()
        self._get_alignment_score(imax, jmax)
        self._trace_back_alignment(imax, jmax)

    def _add_gap_penalties(self) -> None:
        """
        Fills number matrix first row and first column with the gap penalties.
        """
        for j in range(1, len(self.seq1) + 1):
            self._nmatrix[0][j] = 0

        for i in range(1, len(self.seq2) + 1):
            self._nmatrix[i][0] = self._nmatrix[i - 1][0] + self.smatrix.gap

    def _get_last_cell_position(self) -> tuple[int, int]:
        """
        Returns the cell row and column of the last cell in the matrix in which
        the alignment ends. For semiglobal, this is the cell in the last row 
        with the highest score.
        """
        jmax = 0
        max_score = 0
        imax = len(self._nmatrix) - 1
        for jcol in range(0, len(self._nmatrix[0])):
            score = self._nmatrix[imax][jcol]
            if score > max_score:
                jmax = jcol
                max_score = score
        return imax, jmax
    
    def _check_best_score(self, diagscore: int, topscore: int, leftscore: int, irow: int, jcol: int) -> None:
        best_pointer = str()
        best_score = int()
        if (irow + 1, jcol + 1) in self._used_indices:
            best_pointer = ""
            best_score = self._nmatrix[len(self._nmatrix) - 1][jcol] - 1000
        else:
            if diagscore >= topscore:
                if diagscore >= leftscore:
                    best_pointer, best_score = ("diag", diagscore)
                else:
                    best_pointer, best_score = ("left", leftscore)
            else:
                if topscore > leftscore:
                    best_pointer, best_score = ("up", topscore)
                else:
                    best_pointer, best_score = ("left", leftscore)

        self._pmatrix[irow + 1][jcol + 1] = best_pointer
        self._nmatrix[irow + 1][jcol + 1] = best_score

    def _trace_back_alignment(self, irow: int, jcol: int) -> None:
        self._alseq1, self._alseq2 = [], []
        
        while True:
            self._used_indices.append((irow, jcol))
            if self._pmatrix[irow][jcol] == "diag":
                self._alseq1.append(self.seq1[jcol - 1])
                self._alseq2.append(self.seq2[irow - 1])
                if self.seq1[jcol - 1] == self.seq2[irow - 1]:
                    self._identity += 1
                irow -= 1
                jcol -= 1
            elif self._pmatrix[irow][jcol] == "up":
                self._alseq1.append(Gap(self._gap_character))
                self._alseq2.append(self.seq2[irow - 1])
                irow -= 1
            elif self._pmatrix[irow][jcol] == "left":
                self._alseq1.append(self.seq1[jcol - 1])
                self._alseq2.append(Gap(self._gap_character))
                jcol -= 1
            else:
                break
        self._alseq1 = list(reversed(self._alseq1))
        self._alseq2 = list(reversed(self._alseq2))
        self._identity = (self._identity / len(self._alseq1)) * 100
        self._seq1_start = jcol
        self._seq2_start = irow

    def _remove_used_indices(self) -> None:
        pass

class Gap(object):
    def __init__(self, character: str = "-") -> None:
        self.character = character

    def __str__(self) -> str:
        return str(self.character)

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, Gap)