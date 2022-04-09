from Crypto.Util import number
from shamir.dealer import Dealer

p: int = number.getPrime(8)
k: int = number.getRandomRange(0, p)
n: int = 5
t: int = 3

print(p, k)
d: Dealer = Dealer(p, k, n, t)

d.generate_polynomial()
print(d.poly)

d.generate_players()
print([str(p) for p in d.players])

k_restored_linear_equations: int = d.restore_secret_matrix(d.players[:d.t])
print(f'Linear equations: {k_restored_linear_equations = }')

k_restored_lagrange_polynomial: int = d.restore_secret_lagrange(d.players[:d.t])
print(f'Lagrange polynomial: {k_restored_lagrange_polynomial = }')
