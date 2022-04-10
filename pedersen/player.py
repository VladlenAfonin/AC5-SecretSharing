class Player:
    """Player model."""

    def __init__(self, x: int, u: int, w: int) -> None:
        """Initialize new Player.

        :param x: Id basically.
        :type x: int
        :param u: Delta polynomial value.
        :type u: int
        :param w: Gamma polynomial value.
        :type w: int
        """

        self.x: int = x
        self.u: int = u
        self.w: int = w
        self.b: int # Lagrange coefficient.

    def __str__(self) -> str:
        """Return string representation of a Player.

        :returns: String representation of a Player.
        :rtype: str
        """

        return f'(x = {self.x}, (u, w) = {self.u, self.w})'
