import numpy as np
from Crypto.Util import number

import shamir.polynomial as polynomial
import shamir.player as player


def get_det(M):
    M = [row[:] for row in M] # make a copy to keep original M unmodified
    N, sign, prev = len(M), 1, 1
    for i in range(N-1):
        if M[i][i] == 0: # swap with another row having nonzero i's elem
            swapto = next( (j for j in range(i+1,N) if M[j][i] != 0), None )
            if swapto is None:
                return 0 # all M[*][i] are zero => zero determinant
            M[i], M[swapto], sign = M[swapto], M[i], -sign
        for j in range(i+1,N):
            for k in range(i+1,N):
                assert ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) % prev == 0
                M[j][k] = ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) // prev
        prev = M[i][i]
    return sign * M[-1][-1]


class Dealer:
    """Dealer model.
    """

    def __init__(self, p: int, k: int, n: int, t: int) -> None:
        """Initialize new Dealer.

        :param p: Modulus for finite field.
        :type p: int
        :param k: Secret to share. Must be less than p.
        :type k: int
        :param n: Number of players.
        :type n: int
        :param t: Number of players enough to restore the secret.
        :type t: int
        """

        self.p: int = p
        self.k: int = k
        self.n: int = n
        self.t: int = t

    def generate_polynomial(self) -> None:
        """Generates polynomial with random coefficients."""

        coeffs: list[int] = \
            [number.getRandomRange(0, self.p) for _ in range(self.t)]
        coeffs[0] = self.k

        self.poly: polynomial.Polynomial = \
            polynomial.Polynomial(self.t, coeffs, self.p)

    def generate_players(self) -> None:
        """Generates Players to be used for secret restoration."""

        self.players: list[player.Player] = \
            [player.Player(i + 1, self.poly.evaluate(i + 1)) for i in range(self.n)]

    def generate_matrix(self, players_subset: list[player.Player]) -> None:
        """Generates matrix to be used when restoring secret."""

        result = []
        for p in players_subset:
            result += [p.x ** i for i in range(self.t)]
        self.matrix: np.ndarray = np.array(result, dtype=int).reshape(self.t, self.t)

    def generate_vector(self, players_subset: list[player.Player]) -> None:
        """Generates vector of values to be used when restoring secret. Sets up
        b in Ax = b."""

        self.vector: np.ndarray = np.array([p.y for p in players_subset], dtype=int)

    def generate_lagrange_coeffs(self, players_subset: list[player.Player]):
        """Generate lagrange coefficients to use when restoring secret."""

        for i in range(self.t):
            players_subset[i].b = np.prod([pk.x * number.inverse(pk.x - players_subset[i].x, self.p)
                for pk in np.delete(players_subset, i)]) % self.p

    def restore_secret_lagrange(self, players_subset: list[player.Player]) -> int:
        """Restore secret using lagrange polynomial.

        :param players_subset: Players' subset to restore secret with.
        :type players_subset: list[player.Player]

        :returns: Restored key.
        :rtype: int
        """

        self.generate_lagrange_coeffs(players_subset)
        return sum([p.b * p.y for p in players_subset]) % self.p
    
    def restore_secret_matrix(self, players_subset: list[player.Player]) -> int:
        """Restore secret using linear equations.

        :param players_subset: Players' subset to restore secret with.
        :type players_subset: list[player.Player]

        :returns: Restored key.
        :rtype: int
        """

        self.generate_matrix(players_subset)
        self.generate_vector(players_subset)

        matrix_kramer_1: np.ndarray = self.matrix.copy()
        matrix_kramer_1[:, 0] = self.vector.copy()

        det: int = get_det(self.matrix)
        det_1: int = get_det(matrix_kramer_1)

        k: int = (det_1 * number.inverse(det, self.p)) % self.p

        return k
