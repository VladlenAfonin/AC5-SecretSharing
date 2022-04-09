class Polynomial:
    """Polynomial model."""

    def __init__(self, degree: int, coeffs: list[int], p: int) -> None:
        """Initializes a new Polynomial of degree given with coefficients given.

        :param degree: Polynomial degree.
        :type degree: int
        :param coeffs: Polynomial coefficients.
        :type degree: list[int]
        :param p: Modulus.
        :type p: int
        """

        self.degree: int = degree
        self.coeffs: list[int] = coeffs.copy()
        self.p: int = p

    def evaluate(self, x: int) -> int:
        """Evaluate polynomial on the input x.

        :param x: Input value.
        :type x: int

        :returns: Polynomial value on the input x.
        :rtype: int
        """

        result: int = \
            sum([(c * x ** i) for i, c in enumerate(self.coeffs)]) % self.p
        return result

    def __str__(self) -> str:
        """Returns string representation of a polynomial.

        :returns: String representation of a polynomial.
        """

        return ''.join([f'{c}x^{i}  ' for i, c in enumerate(self.coeffs)])
