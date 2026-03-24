"""Algorithm X with Dancing Links (DLX) for exact cover problems.

Array-based implementation following Knuth's design.
Supports primary columns (must be covered exactly once) and
secondary columns (may be covered at most once).

Node layout:
    Node 0          -- root header
    Nodes 1..P      -- primary column headers (P = num_primary)
    Nodes P+1..P+S  -- secondary column headers (S = num_secondary)
    Subsequent nodes -- row elements added by add_row()

Each node stores: L, R, U, D (link indices), C (column index), row_id.
Column headers additionally use the S array for column size.
"""

_INITIAL_CAPACITY = 4096


class DLX:
    """Solve exact cover problems using Dancing Links."""

    __slots__ = (
        "_num_primary",
        "_num_cols",
        "_size",
        "_cap",
        # Per-node parallel arrays
        "L",       # left link
        "R",       # right link
        "U",       # up link
        "D",       # down link
        "C",       # column header index
        "S",       # column size (only meaningful for header nodes)
        "ROW",     # row_id stored on every node in that row
        # Solver state
        "_solutions",
        "_solution_stack",
        "_done",   # flag to stop search early
    )

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, num_primary: int, num_secondary: int = 0):
        """Initialise the DLX structure.

        Columns 0..num_primary-1 are primary (must be covered).
        Columns num_primary..num_primary+num_secondary-1 are secondary.
        """
        if num_primary < 1:
            raise ValueError("Need at least one primary column")

        self._num_primary = num_primary
        self._num_cols = num_primary + num_secondary
        num_headers = 1 + self._num_cols  # root + column headers

        # Pre-allocate flat arrays.
        cap = max(_INITIAL_CAPACITY, num_headers * 2)
        self._cap = cap
        self._size = num_headers

        self.L = [0] * cap
        self.R = [0] * cap
        self.U = [0] * cap
        self.D = [0] * cap
        self.C = [0] * cap
        self.S = [0] * cap
        self.ROW = [None] * cap

        # --- Wire up root (node 0) and column headers (nodes 1..num_cols) ---

        # Primary columns: 1..num_primary  (linked to root in header row)
        # Secondary columns: num_primary+1..num_cols  (linked among themselves,
        #   but NOT reachable from root -- the solver never selects them)

        # Each column header is its own vertical circular list initially.
        for i in range(num_headers):
            self.U[i] = i
            self.D[i] = i
            self.C[i] = i
            self.S[i] = 0

        # Link primary headers: root <-> 1 <-> 2 <-> ... <-> P <-> root
        self.R[0] = 1
        self.L[0] = num_primary
        for i in range(1, num_primary + 1):
            self.L[i] = i - 1
            self.R[i] = i + 1 if i < num_primary else 0

        # Link secondary headers among themselves.
        # They form a separate circular list so cover/uncover still works
        # on them, but the root never reaches them.
        if num_secondary > 0:
            first_sec = num_primary + 1
            last_sec = self._num_cols
            for i in range(first_sec, last_sec + 1):
                self.L[i] = i - 1 if i > first_sec else last_sec
                self.R[i] = i + 1 if i < last_sec else first_sec

        self._solutions: list[list] = []
        self._solution_stack: list = []
        self._done = False

    # ------------------------------------------------------------------
    # Building the matrix
    # ------------------------------------------------------------------

    def add_row(self, columns: list[int], row_id=None):
        """Add a row covering the given column indices.

        *columns* should be a list of 0-based column indices.
        *row_id* is an opaque identifier returned in solutions.
        """
        if not columns:
            return

        # Map from 0-based external column index to internal header node
        # index (header node = column_index + 1).
        first_node = self._size

        for col in columns:
            if col < 0 or col >= self._num_cols:
                raise ValueError(
                    f"Column {col} out of range [0, {self._num_cols})"
                )
            node = self._alloc()
            header = col + 1

            self.C[node] = header
            self.ROW[node] = row_id

            # Insert into the column's vertical circular list (above header).
            up = self.U[header]
            self.U[node] = up
            self.D[node] = header
            self.D[up] = node
            self.U[header] = node

            self.S[header] += 1

        last_node = self._size - 1

        # Link all nodes of this row horizontally in a circular list.
        for i in range(first_node, last_node + 1):
            self.L[i] = i - 1 if i > first_node else last_node
            self.R[i] = i + 1 if i <= last_node - 1 else first_node

    # ------------------------------------------------------------------
    # Solving
    # ------------------------------------------------------------------

    def solve(self, max_solutions: int = 0) -> int:
        """Find all solutions (or up to *max_solutions* if > 0).

        Returns the number of solutions found.  Retrieve them via the
        *solutions* property.

        The solver may be called multiple times on the same instance;
        each call resets solutions and restores the matrix fully.
        """
        self._solutions = []
        self._solution_stack = []
        self._done = False
        self._search(max_solutions)
        return len(self._solutions)

    @property
    def solutions(self) -> list[list]:
        """Return the solutions found by the last call to *solve*.

        Each solution is a list of *row_id* values.
        """
        return self._solutions

    # ------------------------------------------------------------------
    # Core Algorithm X recursion
    # ------------------------------------------------------------------

    def _search(self, max_solutions: int):
        """Recursive backtracking search.

        Sets *_done* flag when *max_solutions* is reached.  Every level
        of recursion checks the flag after returning from a child call
        and unwinds cleanly, so the matrix is always fully restored.
        """
        R = self.R

        # If no primary columns remain, record a solution.
        if R[0] == 0:
            self._solutions.append(list(self._solution_stack))
            if max_solutions > 0 and len(self._solutions) >= max_solutions:
                self._done = True
            return

        # MRV heuristic: choose the primary column with fewest rows.
        col = self._choose_column()

        self._cover(col)

        # Try each row in this column.
        D = self.D
        node = D[col]
        while node != col:
            self._solution_stack.append(self.ROW[node])

            # Cover every other column in this row.
            j = R[node]
            while j != node:
                self._cover(self.C[j])
                j = R[j]

            self._search(max_solutions)

            # Undo: uncover in reverse order.
            self._solution_stack.pop()
            j = self.L[node]
            while j != node:
                self._uncover(self.C[j])
                j = self.L[j]

            if self._done:
                break

            node = D[node]

        self._uncover(col)

    # ------------------------------------------------------------------
    # Column selection (MRV)
    # ------------------------------------------------------------------

    def _choose_column(self) -> int:
        """Return the primary column header with the smallest size."""
        R = self.R
        S = self.S
        best = R[0]
        best_size = S[best]
        col = R[best]
        while col != 0:
            if S[col] < best_size:
                best = col
                best_size = S[col]
                if best_size == 0:
                    break  # can't do better than 0
            col = R[col]
        return best

    # ------------------------------------------------------------------
    # Cover / Uncover
    # ------------------------------------------------------------------

    def _cover(self, col: int):
        """Cover column *col*: remove it from the header list and remove
        all rows that intersect it from their other columns."""
        L = self.L
        R = self.R
        U = self.U
        D = self.D
        C = self.C
        S = self.S

        # Remove column header from horizontal list.
        R[L[col]] = R[col]
        L[R[col]] = L[col]

        # Walk down the column.
        i = D[col]
        while i != col:
            # Walk right across the row.
            j = R[i]
            while j != i:
                # Remove node j from its column's vertical list.
                D[U[j]] = D[j]
                U[D[j]] = U[j]
                S[C[j]] -= 1
                j = R[j]
            i = D[i]

    def _uncover(self, col: int):
        """Uncover column *col* -- exact reverse of *_cover*."""
        L = self.L
        R = self.R
        U = self.U
        D = self.D
        C = self.C
        S = self.S

        # Walk up the column (reverse of cover's downward walk).
        i = U[col]
        while i != col:
            # Walk left across the row (reverse of cover's rightward walk).
            j = L[i]
            while j != i:
                S[C[j]] += 1
                U[D[j]] = j
                D[U[j]] = j
                j = L[j]
            i = U[i]

        # Restore column header in horizontal list.
        R[L[col]] = col
        L[R[col]] = col

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _alloc(self) -> int:
        """Allocate a new node, growing arrays if necessary."""
        idx = self._size
        if idx >= self._cap:
            self._grow()
        self._size += 1
        return idx

    def _grow(self):
        """Double array capacity."""
        old = self._cap
        new = old * 2
        extend = [0] * old
        self.L.extend(extend)
        self.R.extend(extend)
        self.U.extend(extend)
        self.D.extend(extend)
        self.C.extend(extend)
        self.S.extend(extend)
        self.ROW.extend([None] * old)
        self._cap = new
