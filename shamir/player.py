class Player:
    """Player model."""

    def __init__(self, x: int, y: int) -> None:
        """Initialize new Player.

        :param x: Id basically.
        :type x: int
        :param y: Polynomial value.
        :type y: int
        """

        self.x: int = x
        self.y: int = y
        self.b: int # Lagrange coefficient.

    def __str__(self) -> str:
        """Return string representation of a Player.

        :returns: String representation of a Player.
        :rtype: str
        """

        return f'(x = {self.x}, y = {self.y})'
