# coding=utf-8
# for python 3.X
import math


def insert_data():
    print("Notice that: alfa^x = beta (mod p)\n")
    alfa = int(input("Insert a generator element: alfa = "))
    beta = int(input("Insert an element of the chosen group (Zp - {0}): beta = "))
    p = int(input("Inserting a prime number: p = "))
    alfa = alfa % p
    beta = beta % p
    print("")
    return int(p), int(alfa), int(beta)


def eo(p, beta):
    """
    Even or Odd algorithm
    :param p: prime
    :param beta: an element
    :return: nothing
    """
    constant = beta

    print("i = 1, a = " + str(beta))

    for i in range(2, int((int(p)-1)/2 + 1)):
        beta = (beta * constant) % p
        if beta == 0:
            beta = 1
        print("i = " + str(i) + ", a = " + str(beta))

    print("")

    if beta == p - 1:
        beta = -1

    if beta == 1:
        print("The discrete logarithm is even!\n")
    elif beta == -1:
        print("The discrete logarithm is odd!\n")
    else:
        print("The discrete logarithm could not exist...\n")


def er(p, alfa, beta):
    """
    Exhaustive research
    :param p: prime
    :param alfa: generator
    :param beta: element of the group
    :return: the discrete logarithm
    """
    if alfa == beta:
        return 1

    constant = alfa
    alfa = 1
    for i in range(1, p):

        alfa = (alfa * constant) % p
        if alfa == 0:
            alfa = 1
        print("i = " + str(i) + ", alfa = " + str(alfa))

        if alfa == beta:
            return i

    return -1


def insert_value(table, value):
    """
    Inserting a value in the table
    :param table: the table where insert the value
    :param value: the value
    :return: the table updated
    """

    table.append(value)


def egcd(a, b):
    """
    Extended Greatest Common Divisor algorithm
    :param a: first element
    :param b: second element
    :return: (g, x, y) such that a^x + b^y = g = gcd(a, b)
    """

    if a == 0:
        return b, 0, 1
    else:
        b_div_a, b_mod_a = divmod(b, a)
        g, x, y = egcd(b_mod_a, a)
        print(str(g) + ", " + str(y - b_div_a * x) + ", " + str(x))

        return g, y - b_div_a * x, x


def mod_inv(a, b):
    """
    :param a: element of the group
    :param b: the order of the group
    :return: x such that (x * a) % b == 1
    """

    g, x, _ = egcd(a, b)

    if g != 1:
        raise Exception

    return x % b


def extended_euclide_algorithm(value, prime):
    """
    :param value: element
    :param prime: prime element
    :return: the inverse of value element in modulo arithmetic
    """
    x = -1
    print("Computing the Extended Euclide algorithm for the inverse of " + str(value) + " (mod " + str(prime) + ")")
    try:
        x = mod_inv(value, prime)
    except Exception:
        print("gcd(a,b) != 1")
    return x


def search_value_in_my_table(table, value):
    """
    Searching an element in a list of elements
    :param table: list of list of two elements
    :param value: the value to find at the first "sub" element
    :return: if the value is present, the relative index
    """
    ok = -1
    index = -1
    for index in range(len(table)):
        if value == table[index][0]:
            ok = 1
            break

    return ok, index


def bsgs(p, alfa, beta):
    """
    Baby-step Giant-step algorithm
    :param p: the prime number
    :param alfa: the generator
    :param beta: an element of the group
    :return: the discrete logarithm
    """
    if alfa == beta:
        return 1

    constant = alfa

    """1) set m = ceil(square root of p)"""

    m = int(math.ceil(math.sqrt(p)))
    print("set m = " + str(m) + "\n")

    """2) building the table"""

    my_table = []
    insert_value(my_table, [1, 0])
    alfa = alfa % p
    insert_value(my_table, [alfa, 1])

    for j in range(2, m):
        alfa = (alfa * constant) % p
        insert_value(my_table, [alfa, j])

    """sorting the table"""

    my_table.sort()
    print("[alfa^j, j]")
    for o in range(len(my_table)):
        print(my_table[o])

    print("")

    """3) inverse of alfa^m"""

    alfa = constant
    for k in range(2, m+1):
        alfa = (alfa * constant) % p
    print("alfa^m = " + str(alfa) + "\n")
    alfa_m_inverted = extended_euclide_algorithm(alfa, p)

    if alfa_m_inverted == -1:
        print("[Failure] alfa_m_inverse does not exists...")
        return -1

    c = beta
    print("\nInverse of alfa^m: " + str(alfa_m_inverted) + "\n")

    """4) searching c in the table"""

    for i in range(0, m):

        ok, index = search_value_in_my_table(my_table, c)
        print("i = " + str(i) + " c = " + str(c) + " index = " + str(index) + "\n")

        if ok == 1:
            return (i * m + my_table[index][1]) % p

        c = (c * alfa_m_inverted) % p
        if c == 0:
            c = 1

    return -1


def xab_sequence(x, a, b, alfa, beta, p, q):
    """
    Iteration of Pollard's Rho
    :param x: element of the x sequence
    :param a: element of the a sequence
    :param b: element of the b sequence
    :param alfa: the generator
    :param beta: an element of the group chosen by the user
    :param p: the prime number
    :param q: (p - 1) / 2
    :return: x, a, b
    """

    sets = x % 3
    # S2
    if sets == 0:
        x = (alfa * x) % p
        a = (a + 1) % q
        # b = b
    # S1
    elif sets == 1:
        x = (beta * x) % p
        # a = a
        b = (b + 1) % q
    # S3
    elif sets == 2:
        x = (x * x) % p
        a = (2 * a) % q
        b = (2 * b) % q

    return int(x), int(a), int(b)


def rho(p, alfa, beta):
    """
    Pollard's Rho algorithm
    :param p: the prime number chosen by the user
    :param alfa: the generator element
    :param beta: an element of the group
    :return: the discrete logarithm
    """
    if alfa == beta:
        return 0
    q = (p - 1) / 2
    q = int(q)
    x = alfa * beta
    a = 1
    b = 1

    xx = x
    aa = a
    bb = b

    # Floyd's algorithm for cycle detention
    for i in range(1, p):

        # tortoise
        x, a, b = xab_sequence(x, a, b, alfa, beta, p, q)
        print("Tortoise = " + str(x) + ", " + str(a) + ", " + str(b))

        # rabbit
        xx, aa, bb = xab_sequence(xx, aa, bb, alfa, beta, p, q)
        xx, aa, bb = xab_sequence(xx, aa, bb, alfa, beta, p, q)
        print("Rabbit = " + str(xx) + ", " + str(aa) + ", " + str(bb))

        if xx == x:
            break

    r = (bb - b) % q
    print("r = " + str(r))

    if r == 0:
        print("[Failure] B - b = 0")
        return -1
    else:
        r_inverse = extended_euclide_algorithm(r, q)
        if r_inverse == -1:
            print("[Failure] r_inverse does not exist...")
            return -1

    print("r_inverse = " + str(r_inverse) + "\n")
    dl = int(r_inverse * (a - aa) % q)

    if pow(int(alfa), int(dl), int(p)) == beta:
        return dl
    if pow(int(alfa), int(dl + q), int(p)) == beta:
        return dl + q

    return -1


def get_factorization_list(number_n):
    number = number_n
    list_factor = []

    while number % 2 == 0:
        list_factor.append(2)
        number = number / 2

    for i in range(3, int(math.ceil(number_n/2)), 2):
        while number % i == 0:
            list_factor.append(i)
            number = number / i

    if number > 2:
        list_factor.append(int(number))

    return list_factor


def gauss_algorithm(x, factor_exp_factorization_list, p):
    y = 0
    print("Gauss' algorithm:" + "\n")
    print("x = ")
    print(x)
    print("factor_exp_factorization_list = ")
    print(factor_exp_factorization_list)
    print("n = " + str(p - 1))

    for i in range(len(x)):
        ni = pow(factor_exp_factorization_list[i][0], factor_exp_factorization_list[i][1])
        print("ni = " + str(ni))
        big_n_i = int((p - 1) / ni)
        print("Ni = " + str(big_n_i))
        big_m_i = int(extended_euclide_algorithm(big_n_i, ni)) % ni
        print("Mi = " + str(big_m_i))
        y += int(x[i]) * big_n_i * big_m_i
        print("y = " + str(y))
    return y % (p - 1)


def ph(p, alfa, beta):
    """
    Pohlig-Hellman
    :param p: the prime number
    :param alfa: the generator
    :param beta: an element of the group
    :return: the discrete logarithm
    """
    if alfa == beta:
        return 1
    """1) Computing the factorization of p"""

    factorization_list = get_factorization_list(p - 1)
    # organizing the list
    factor_exp_factorization_list = []
    i = 0
    last_number = factorization_list[0]

    for k in range(len(factorization_list)):
        if factorization_list[k] == last_number:
            i = i + 1
        else:
            factor_exp_factorization_list.append([last_number, i])
            last_number = factorization_list[k]
            i = 1
        if k == len(factorization_list) - 1:
            factor_exp_factorization_list.append([last_number, i])

    print("factorization of n (= " + str(p) + " - 1):")
    print(factor_exp_factorization_list)

    """2) Computing xi = lo + ... + lei-1 * pi^ei-1"""

    x = []
    for r in range(0, len(factor_exp_factorization_list)):

        """2.a)"""

        q = factor_exp_factorization_list[r][0]
        e = factor_exp_factorization_list[r][1]

        print("q = " + str(q))
        print("e = " + str(e))
        """2.b)"""

        gamma = 1
        l = []

        """2.c)"""

        alfa_bar = pow(int(alfa), int((p - 1)/q), int(p))
        print("alfa_bar = " + str(alfa_bar))

        """2.d) calculating lj"""

        for j in range(0, e):

            """2.d.i)"""
            if j == 0:
                gamma = 1
            else:
                gamma = gamma * pow(alfa, l[j-1] * pow(q, j-1)) % p
            print("gamma = " + str(gamma))

            if j == 0:
                gamma_inverse = 1
            else:
                gamma_inverse = extended_euclide_algorithm(gamma, p)
            print("gamma_inverse = " + str(gamma_inverse))

            beta_bar = pow(beta * gamma_inverse, (p - 1) / pow(q, j+1)) % p
            print("beta_bar = " + str(beta_bar))

            """2.d.ii)"""

            print("p = " + str(p) + " alfa_bar = " + str(alfa_bar) + " beta_bar = " + str(beta_bar))
            if p < 100000:
                number = er(p, alfa_bar, beta_bar)
            else:
                number = bsgs(p, alfa_bar, beta_bar)

            print("number = " + str(number))
            if number == -1:
                return -1
            l.append(number)
            print("l[" + str(j) + "] = " + str(l[j]))

        """2.e)"""

        x.append(0)
        for giro in range(e):
            x[r] += int(l[giro] * pow(q, giro))
        if x[r] == 0:
            x[r] = 1
        print("x[" + str(r) + "] = " + str(x[r]))

    """3) Gauss' algorithm to compute x = xi (mod pi^ei)"""
    dl = gauss_algorithm(x, factor_exp_factorization_list, p)

    return dl


def algorithm(number):
    dl = -1
    if number == 0:
        print("You choose to exit! Goodbye!")
        exit(0)
    elif number == 1:
        print("You choose Even or Odd\n")
        (p, alfa, beta) = insert_data()
        eo(p, beta)
        dl = -2
    elif number == 2:
        print("You choose Exhaustive Research\n")
        (p, alfa, beta) = insert_data()
        dl = er(p, alfa, beta)
    elif number == 3:
        print("You choose Baby-Step Giant-Step\n")
        (p, alfa, beta) = insert_data()
        dl = bsgs(p, alfa, beta)
    elif number == 4:
        print("You choose Pohlig-Hellman\n")
        (p, alfa, beta) = insert_data()
        dl = ph(p, alfa, beta)
    elif number == 5:
        print("You choose Pollard's Rho\n")
        print("Keep on mind that q = (p - 1) / 2 has to be a prime number !!\n")
        (p, alfa, beta) = insert_data()
        dl = rho(p, alfa, beta)
    else:
        print("Something was wrong!\n")
        dl = -2

    if dl == -1:
        print("The discrete logarithm could not exist..." + "\n")
    elif dl != -2:
        print("The discrete logarithm is x = " + str(dl) + "\n")


if __name__ == '__main__':
    print("Algorithms for the problem of discrete logarithm\n<<<<<<  Logica e Algebra 2  >>>>>>\n")

    while True:
        print("These are the implemented algorithms:\n")
        print('0: Exit')
        print("1: Even or Odd")
        print("2: Exhaustive Research")
        print("3: Baby-Step Giant-Step")
        print("4: Pohlig-Hellman")
        print("5: Pollard's Rho\n")

        choice = input("Which algorithm would you like to choose? ")
        print("")
        try:
            algorithm(int(choice))
        except ValueError:
            print("\nChars are not admitted...\n")
