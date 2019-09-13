import time

class KitaevHamiltonian:
    def __init__(_N: int, J: DisorderParameter):
        self.blocks = []
        for Q in range(N):
            self.blocks.append(KitaevHamiltonianBlock(N, Q, J))

    def __del__(self):
        for Q in range(N):
            if self.blocks[Q] != 0:
                del self.blocks[Q]
                self.blocks[Q] = 0

    def dim(self) -> int:
        return int(pow(2, N))

    def diagonalize(self, full_diagonalization: bool,
                    print_progress: bool):
        for Q in range(N):
            if print_progress:
                print("Diagonalizing Q=", Q, " ... ")

        tic = time.time()
        self.blocks[Q].diagonalize(full_diagonalization)

        if print_progress:
            print("took ", int(time.time() - tic), " seconds")

        if print_progress:
            print("Diagonalization complete.")

    def is_diagonalized(self) -> bool:
        return self.blocks[0].diagonalized

    def eigenvalues(self) -> np.ndarray:
        assert self.is_diagonalized()
        evs = np.zeros(self.dim())
        k = 0

        for Q in range(N):
            block_dim = self.blocks[Q].dim()
            # ???
            evs.block(k, 0, block_dim, 1) = self.blocks[Q].eigenvalues()
            k += block_dim
        return evs

    def as_matrix(self) -> Mat:
        H: Mat = np.zeros(self.dim(), self.dim())

        int block_row = 0
        int block_col = 0

        for Q in range(N):
            block_dim: int = self.blocks[Q].dim()
            H.block(block_row, block_col, block_dim, block_dim) = self.blocks[Q].matrix
            block_row += block_dim
            block_col += block_dim

        return H
