from mpmath import *
import matplotlib.pyplot as plt


c = 3 * (10 ** 8)  # speed of light
l = 400 * (10 ** (-9))  # start wavelength
delta = 1 * 10 ** (-9)  # delta of wavelength
R = 100 * (10 ** (-9))  # Radius of sphere
i = 0  # counter
ind = 1

f = open('si_n.txt', 'r')
u = open('si_k.txt', 'r')
thx = open('theor_x.txt', 'r')
thy = open('theor_y.txt', 'r')
rsl = open ('res.txt', 'w')
rsl1 = open ('res_x.txt', 'r')
rsl2 = open ('res_y.txt', 'r')

result_y = []
result_x = []
theory_y = []
theory_x = []


for line in thx:
    theory_x.append(float(line))
for line in thy:
    theory_y.append(float(line))


if (ind == 0):
    while i < 801:
        if (i % 10 == 0 and i != 800):
            if (i == 0):  # interpolation
                n1 = float(f.readline()) + j * float(u.readline())  # reading complex refractive index
                nr = n1
                n1n = float(f.readline()) + j * float(u.readline())
            else:
                n1 = n1n
                nr = n1
                n1n = float(f.readline()) + j * float(u.readline())
        elif (i == 800):
            n1 = n1n
        else:
            n1 = n1 + (n1n - nr) / 10

        w = 2 * pi * c / l  # frequency
        k = w * 1 / c  # wave number in air
        k1 = w * n1 / c  # wave number in material
        r = k * R
        r1 = k1 * R

        jn = lambda num, x: sqrt(pi / (2 * x)) * besselj(num + 0.5, x)  # spherical Bessel functions of the first kind
        jns = lambda num, x: sqrt(pi / (2 * x)) * besselj(num + 0.5, x, 1) - (sqrt(2 * pi / x) / (4 * x)) * besselj(
        num + 0.5, x)  # derivative of the spherical Bessel function of the first kind

        yn = lambda num, x: sqrt(pi / (2 * x)) * bessely(num + 0.5, x)  # spherical Bessel functions of the second kind
        yns = lambda num, x: sqrt(pi / (2 * x)) * bessely(num + 0.5, x, 1) - (sqrt(2 * pi / x) / (4 * x)) * bessely(
        num + 0.5, x)  # derivative of the spherical Bessel function of the second kind

        hn = lambda num, x: sqrt(pi / (2 * r)) * besselj(num + 0.5, r) + j * sqrt(pi / (2 * r)) * bessely(num + 0.5, r)
        # spherical Hankel functions of the first and second kind Bessel functions
        hns = lambda num, x: sqrt(pi / (2 * r)) * besselj(num + 0.5, r, 1) - (sqrt(2 * pi / r) / (4 * r)) * besselj(
        num + 0.5, r) + j * (sqrt(pi / (2 * r)) * bessely(num + 0.5, r, 1) - (sqrt(2 * pi / r) / (4 * r)) * bessely(
        num + 0.5, r))  # derivative of the spherical Hankel functions of the first and second kind Bessel functions

        an = lambda num: ((n1 ** 2) * (jn(num, r) + r * jns(num, r)) * jn(num, r1) - (jn(num, r1) + r1 * jns(num, r1)) *
                      jn(num, r)) / ((n1 ** 2) * (hn(num, r) + r * hns(num, r)) * jn(num, r1) - (
                jn(num, r1) + r1 * jns(num, r1)) * hn(
        num, r))  # scattering coefficient

        bn = lambda num: ((jn(num, r) + r * jns(num, r)) * jn(num, r1) - (jn(num, r1) + r1 * jns(num, r1)) * jn(num, r)) / (
            (hn(num, r) + r * hns(num, r)) * jn(num, r1) - (jn(num, r1) + r1 * jns(num, r1)) * hn(num, r))
        # scattering coefficient

        with extradps(50):
            s = (2 * pi / k ** 2) * nsum(lambda n: (2 * n + 1) * (fabs(an(n)) ** 2 + fabs(bn(n)) ** 2), [1, inf])
            # scattering cross section
            result_y.append(s*10**13)  # print result for this wavelength
        result_x.append(l*10**9)
        rsl.write(str(l) + ' ' + str(s) + '\n')
        l = l + delta  # change wavelength
        i = i + 1  # change counter


if (ind == 1):
    for line in rsl1:
        result_x.append(float(line))
    for line in rsl2:
        result_y.append(float(line))


plt.plot(result_x, result_y, label = 'Theory')
plt.plot(theory_x, theory_y, 'r+', label = 'COMSOL model')
plt.ylabel('Scattering cross-section, m^2')
plt.xlabel('Wavelength, nm')
plt.legend()
plt.grid()
plt.show()
