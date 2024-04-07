class Curve:
    """
    A class used to represent an Elliptic Curve

    ...

    Attributes
    ----------
    a : int
        The 'a' coefficient of the elliptic curve equation
    b : int
        The 'b' coefficient of the elliptic curve equation
    p : int
        The prime number used for the finite field

    Methods
    -------
    add(P, Q):
        Adds two points P and Q on the elliptic curve
    """

    def __init__(self, a, b, p, G):
        """
        Parameters
        ----------
        a : int
            The 'a' coefficient of the elliptic curve equation
        b : int
            The 'b' coefficient of the elliptic curve equation
        p : int
            The prime number used for the finite field
        G : (int, int)
            The generator point of the elliptic curve
        """

        self.a = a
        self.b = b
        self.p = p
        self.G = G
        self.__order__()

    def modinv(self, n):
        """
        Calculate the modular inverse of n mod p

        Parameters
        ----------
        n : int
            The number to find the modular inverse of

        Returns
        -------
        int
            The modular inverse of n mod p
        """

        return pow(n, self.p - 2, self.p)

    def __correct__(self, n):
        """
        Correct the value of n to be within the range of 0 and p

        Parameters
        ----------
        n : int
            The number to correct

        Returns
        -------
        int
            The corrected value of n
        """

        return (n % self.p + self.p) % self.p

    def add(self, P, Q):
        """
        Adds two points P and Q on the elliptic curve

        Parameters
        ----------
        P : tuple
            A point on the elliptic curve
        Q : tuple
            Another point on the elliptic curve

        Returns
        -------
        tuple
            The result of the addition of P and Q
        """

        # If P is the point at infinity, return Q
        if P is None:
            return Q
        # If Q is the point at infinity, return P
        if Q is None:
            return P

        x1, y1 = P
        x2, y2 = Q

        # If P = Q (point doubling)
        if P == Q:
            if y1 == 0:
                # This would mean that 2P is the point at infinity
                return None
            # Calculate lambda for point doubling
            lamda = (3 * x1 * x1 + self.a) * self.modinv(2 * y1) % self.p
        else:
            # If P and Q are inverses of each other
            if x1 == x2 and y1 != y2:
                return None
            # Calculate lambda for the addition of two different points
            lamda = self.__correct__(y2 - y1) * self.modinv(x2 - x1) % self.p

        x3 = (lamda**2 - x1 - x2) % self.p
        y3 = (lamda * (x1 - x3) - y1) % self.p

        return (x3, y3)

    def __order__(self):
        """
        Calculate the order of a point G on the elliptic curve

        Parameters
        ----------
        G : tuple
            A point on the elliptic curve

        Returns
        -------
        int
            The order of the point G
        """

        P = self.G
        order = 1
        while P is not None:
            P = self.add(P, self.G)
            order += 1

        self.order = order

    def product(self, n, P=None):

        n = n % self.order

        """
        Calculate the product of a point P and an integer n

        Parameters
        ----------
        P : tuple
            A point on the elliptic curve
        n : int
            The integer to multiply the point P by

        Returns
        -------
        tuple
            The result of the multiplication of P by n
        """

        Q = None
        R = self.G if P is None else P
        while n > 0:
            if n % 2 == 1:
                Q = self.add(Q, R)
            R = self.add(R, R)
            n = n // 2

        return Q
