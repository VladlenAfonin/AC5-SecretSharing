import numpy as np
from Crypto.Util import number

import pedersen.polynomial as polynomial
import pedersen.player as player


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

    def __init__(self, qBits, pBits, n: int, t: int) -> None:
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

        self.generate_params(qBits, pBits)
        self.n: int = n
        self.t: int = t

    def generate_params(self, qBits: int, pBits: int) -> None:
        """Additional parameters generation.

        :param qBits: Number of bits q must be.
        :type qBits: int
        :param pBits: Number of bits p must be.
        :type pBits: int
        """

        self.q = number.getPrime(qBits)
        self.p = number.getRandomInteger(pBits - qBits) * self.q + 1
        while not number.isPrime(self.p):
            self.q = number.getPrime(qBits)
            self.p = number.getRandomInteger(pBits - qBits) * self.q + 1

        self.k = number.getRandomRange(0, self.q)
        self.g = number.getRandomNBitInteger(pBits) ** ((self.p - 1) // self.q) % self.p
        self.h = number.getRandomNBitInteger(pBits) ** ((self.p - 1) // self.q) % self.p

    def generate_pedersen(self):
        """Generate coefficients for pedersen check.

        :returns: Pedersen coefficients.
        :rtype: list[int]
        """

        return [self.g ** self.delta.coeffs[i] * self.h ** self.gamma.coeffs[i] % self.p for i in range(self.t)]

    def check_pedersen(self, players: list[player.Player]) -> list[bool]:
        """Perform Pedersen check.

        :param players: List of players who will verify.
        :type players: list[player.Player]

        :returns: List of check results.
        :rtype: list[bool]
        """

        pedersen_params: list[int] = self.generate_pedersen()

        # print(pedersen_params)
        # for p in players:
        #     print(
        #         f'{(self.g ** p.u * self.h ** p.w) % self.p} == {np.prod([ pp ** (p.x ** i) % self.p for i, pp in enumerate(pedersen_params) ], dtype=int) % self.p}'
        #     )

        results = [
            (self.g ** p.u * self.h ** p.w) % self.p ==
                np.prod([
                    pp ** (p.x ** i) % self.p for i, pp in enumerate(pedersen_params)
                ], dtype=int) % self.p
            for p in players
        ]

        return results

    def generate_polynomial(self) -> None:
        """Generates polynomial with random coefficients."""

        coeffs: list[int] = \
            [number.getRandomRange(0, self.q) for _ in range(self.t)]
        coeffs[0] = self.k

        self.delta: polynomial.Polynomial = \
            polynomial.Polynomial(self.t, coeffs, self.q)

        coeffs: list[int] = \
            [number.getRandomRange(0, self.q) for _ in range(self.t)]

        self.gamma: polynomial.Polynomial = \
            polynomial.Polynomial(self.t, coeffs, self.q)

    def generate_players(self) -> None:
        """Generates Players to be used for secret restoration."""

        self.players: list[player.Player] = \
            [player.Player(i + 1, self.delta.evaluate(i + 1), self.gamma.evaluate(i + 1)) for i in range(self.n)]

    def generate_matrix(self, players_subset: list[player.Player]) -> None:
        """Generates matrix to be used when restoring secret."""

        result = []
        for p in players_subset:
            result += [p.x ** i for i in range(self.t)]
        self.matrix: np.ndarray = np.array(result, dtype=int).reshape(self.t, self.t)

    def generate_vector(self, players_subset: list[player.Player]) -> None:
        """Generates vector of values to be used when restoring secret. Sets up
        b in Ax = b."""

        self.vector: np.ndarray = np.array([p.u for p in players_subset], dtype=int)

    def generate_lagrange_coeffs(self, players_subset: list[player.Player]):
        """Generate lagrange coefficients to use when restoring secret."""

        for i in range(self.t):
            players_subset[i].b = np.prod([pk.x * number.inverse(pk.x - players_subset[i].x, self.q)
                for pk in np.delete(players_subset, i)]) % self.q

    def restore_secret_lagrange(self, players_subset: list[player.Player]) -> int:
        """Restore secret using lagrange polynomial.

        :param players_subset: Players' subset to restore secret with.
        :type players_subset: list[player.Player]

        :returns: Restored key.
        :rtype: int
        """

        self.generate_lagrange_coeffs(players_subset)
        return sum([p.b * p.u for p in players_subset]) % self.q
    
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

        k: int = (det_1 * number.inverse(det, self.q)) % self.q

        return k
