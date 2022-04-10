import numpy as np
from Crypto.Util import number

from shamir import player
from shamir import dealer


def feldman_check(pl: player.Player, g: int, q: int, feldman_params: list[int]) -> bool:
    """Perform Feldman check.

    :param g: Base.
    :type g: int
    :param p: Modulus.
    :type p: int
    :param feldman_params: List of Feldman parameters.
    :type feldman_params: list[int]

    :returns: True if check was successful.
    :rtype: bool
    """

    # print(f'{g ** pl.y % q} == {np.prod([f ** (pl.x ** i) % q for i, f in enumerate(feldman_params)], dtype=int) % q}')
    return g ** pl.y % q == np.prod([f ** (pl.x ** i) % q for i, f in enumerate(feldman_params)], dtype=int) % q


class Dealer(dealer.Dealer):
    def generate_params(self, pBits: int, qBits: int) -> None:
        """Generate parameters suitable for Feldman check.

        :param qBits: Number of bits q must be.
        :type qBits: int
        :param pBits: Number of bits p must be.
        :type pBits: int
        """

        self.p = number.getPrime(pBits)
        self.q = number.getRandomInteger(qBits - pBits) * self.p + 1

        while not number.isPrime(self.q):
            self.p = number.getPrime(pBits)
            self.q = number.getRandomInteger(qBits - pBits) * self.p + 1

        self.k = number.getRandomRange(0, self.p)
        self.g = number.getRandomNBitInteger(qBits) ** ((self.q - 1) // self.p) % self.q

    def generate_feldman_params(self) -> list[int]:
        """Generate Feldman parameters.

        :returns: List of Feldman parameters.
        :rtype: list[int]
        """

        # print([self.g ** c % self.q for c in self.poly.coeffs])
        return [self.g ** c % self.q for c in self.poly.coeffs]


def Demo(fail: bool=False):
    n: int = 5
    t: int = 3

    d: Dealer = Dealer(0, 0, n, t)

    # Here the key and modulus are generated.
    d.generate_params(8, 10)
    print(f'{d.p = }, {d.q = }, {d.k = }, {d.g = }')

    d.generate_polynomial()
    print(f'a(x) = {d.poly}')

    d.generate_players()
    print(f'Players: {[str(p) for p in d.players]}')

    if fail: d.players[0].y = 1

    checks: list[bool] = [feldman_check(p, d.g, d.q, d.generate_feldman_params()) for p in d.players[:t]]

    print(f'Check results: {checks}')

    k_restored_linear_equations: int = d.restore_secret_matrix(d.players[:d.t])
    print(f'Linear equations: {k_restored_linear_equations}')

    k_restored_lagrange_polynomial: int = d.restore_secret_lagrange(d.players[:d.t])
    print(f'Lagrange polynomial: {k_restored_lagrange_polynomial}')
