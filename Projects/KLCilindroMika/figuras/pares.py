"""pares.py: os pares de pontos aleatorios do KLTaludeCeron (talude.cpp, struct Sorteio) refeitos em Python:
std::mt19937_64 com a semente 1, U = (g >> 11) 2^-53, 500 pares uniformes no talude (rejeicao no retangulo
[0,25] x [0,10]) e 500 pares proximos (o segundo ponto a ate Lx/10 em x e Ly/4 em y do primeiro), na mesma ordem."""
M64 = 0xFFFFFFFFFFFFFFFF


class MT64:
    def __init__(self, seed):
        self.mt = [0] * 312
        self.i = 312
        self.mt[0] = seed & M64
        for k in range(1, 312):
            self.mt[k] = (6364136223846793005 * (self.mt[k - 1] ^ (self.mt[k - 1] >> 62)) + k) & M64

    def __call__(self):
        if self.i >= 312:
            mt = self.mt
            for k in range(312):
                x = (mt[k] & 0xFFFFFFFF80000000) | (mt[(k + 1) % 312] & 0x7FFFFFFF)
                xa = x >> 1
                if x & 1:
                    xa ^= 0xB5026F5AA96619E9
                mt[k] = mt[(k + 156) % 312] ^ xa
            self.i = 0
        y = self.mt[self.i]
        self.i += 1
        y ^= (y >> 29) & 0x5555555555555555
        y ^= (y << 17) & 0x71D67FFFEDA60000
        y ^= (y << 37) & 0xFFF7EEE000000000
        y ^= y >> 43
        return y & M64


def no_talude(x, y):
    if x < 0 or x > 25 or y < 0 or y > 10:
        return False
    return y <= 5 or x + y <= 20


def pontos(npares=500, semente=1, Lx=20., Ly=2.):
    """4 npares pontos: os pares uniformes (0..2n-1) e os proximos (2n..4n-1), (a, b) consecutivos"""
    g = MT64(semente)

    def U():
        return (g() >> 11) * (1.0 / 9007199254740992.0)

    def unif():
        while True:
            x = 25 * U()
            y = 10 * U()
            if no_talude(x, y):
                return (x, y)
    pts = []
    for _ in range(npares):
        pts += [unif(), unif()]
    for _ in range(npares):
        a = unif()
        while True:
            b = (a[0] + Lx / 10 * (2 * U() - 1), a[1] + Ly / 4 * (2 * U() - 1))
            if no_talude(*b):
                break
        pts += [a, b]
    return pts


if __name__ == "__main__":
    g = MT64(5489)
    for _ in range(9999):
        g()
    print("mt19937_64: a 10000a saida com a semente 5489 deve ser 9981545732273789042:", g())
