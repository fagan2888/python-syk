import numpy as np

class KitaevDisorderParameter:
    def __init__(self,
                 N: int, J: float,
                 complex_elements: bool) -> None:
        self.Jelems = np.zeros((N, N, N, N), dtype=complex)

        offdiag = J / np.sqrt(2)

        # Choose random J elements
        for i in range(N):
            for j in range(i + 1, N):
                for k in range(N):
                    for l in range(k + 1, N):
                        if (k == i and l == j):
                            # Diagonal element is real
                            self.Jelems[i][j][i][j] = complex(np.random.normal(0, J))
                        else:
                            if (complex_elements):
                                # Off-diagonal element is complex
                                self.Jelems[i][j][k][l] = complex(
                                        np.random.norma(0, offdiag) + 1j * np.random.normal(0, offdiag))
                            else:
                                self.Jelems[i][j][k][l] = complex(np.random.normal(0, offdiag))

                        # Apply reality condition
                        self.Jelems[k][l][i][j] = np.conj(self.Jelems[i][j][k][l])

                        # Anti-symmetrize (both Jijkl and Jklij)
                        self.Jelems[j][i][k][l] = - self.Jelems[i][j][k][l]
                        self.Jelems[i][j][l][k] = - self.Jelems[i][j][k][l]
                        self.Jelems[j][i][l][k] = + self.Jelems[i][j][k][l]

                        self.Jelems[l][k][i][j] = - self.Jelems[k][l][i][j]
                        self.Jelems[k][l][j][i] = - self.Jelems[k][l][i][j]
                        self.Jelems[l][k][j][i] = + self.Jelems[k][l][i][j]

    def elem(self, i: int, j: int, k: int, l: int) -> complex:
        return self.Jelems[i][j][k][l]
