from Crypto.Util import number
from pedersen.dealer import Dealer


def Demo() -> None:
    qBits: int = 8
    pBits: int = 10
    n: int = 5
    t: int = 3

    d: Dealer = Dealer(qBits, pBits, n, t)

    print(f'{d.p = }, {d.q = }, {d.k = }, {d.g = }, {d.h = }')

    d.generate_polynomial()
    print(f'delta(x) = {d.delta}')
    print(f'gamma(x) = {d.gamma}')

    d.generate_players()
    print(f'Players: {[str(p) for p in d.players]}')

    results = d.check_pedersen(d.players[:t])
    print(f'Check results: {results}')

    k_restored_linear_equations: int = d.restore_secret_matrix(d.players[:d.t])
    print(f'Linear equations: {k_restored_linear_equations}')

    k_restored_lagrange_polynomial: int = d.restore_secret_lagrange(d.players[:d.t])
    print(f'Lagrange polynomial: {k_restored_lagrange_polynomial}')
