from __future__ import division
import itertools
import math

from cycler import cycler
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import fssa
"""https://pyfssa.readthedocs.io/en/stable/tutorial.html"""


def mock_scaling_f(x):
    """Mock scaling function"""
    return np.exp(-(x + 1.0)**2)


def mock_scaled_data(l, rho, rho_c=0.5, nu=2.5, zeta=1.5):
    """Generate scaled data from mock scaling function"""
    return np.transpose(
        np.power(l, zeta / nu) *
        mock_scaling_f(
            np.outer(
                (rho - rho_c), np.power(l, 1 / nu)
            )
        )
    )


x = np.linspace(-4.0, 2.0, num=200)

fig, ax = plt.subplots()
ax.plot(x, mock_scaling_f(x), label=r'$\tilde{f}(x)$', rasterized=True)
ax.set_xbound(x.min(), x.max())
ax.set_ybound(0.0, 1.1)
ax.set_xlabel(r'$x$')
ax.legend()
plt.show()


rhos = np.linspace(-0.5, 0.8, num=200)
ls = [math.floor(l) for l in np.logspace(1, 3, num=5)]
print(ls)

palette = sns.cubehelix_palette(
    n_colors=len(ls), start=2.0, rot=0.35, gamma=1.0, hue=1.0, light=0.6, dark=0.2,
)
sns.palplot(palette)

a = mock_scaled_data(ls, rhos)

fig, ax = plt.subplots()
ax.set_prop_cycle(cycler('color', palette))
for l_index, l in enumerate(ls):
    ax.plot(
        rhos, a[l_index, :],
        '.',
        label=r'${}$'.format(l),
        rasterized=True,
    )
ax.set_xbound(rhos.min(), rhos.max())
ax.set_xlabel(r'$\rho$')
ax.legend(title=r'$L$', loc='upper left')
plt.show()


da = a * 0.1

rho_c = np.tile(0.5, (3, 3))
nu = np.tile(2.5, (3, 3))
zeta = np.tile(1.5, (3, 3))
rho_c[0, :] = [0.25, 0.5, 0.75]
nu[1, :] = [2.0, 2.5, 3.0]
zeta[2, :] = [1.0, 1.5, 2.0]

# re-scale data (manually)
scaled_data = list()
quality = list()
for i in range(3):
    my_scaled_data = list()
    my_quality = list()
    for j in range(3):
        my_scaled_data.append(
            fssa.scaledata(
                ls, rhos, a, da,
                rho_c[i, j], nu[i, j], zeta[i, j]
            )
        )
        my_quality.append(fssa.quality(*my_scaled_data[-1]))
    scaled_data.append(my_scaled_data)
    quality.append(my_quality)

# plot manually re-scaled data
fig, axes = plt.subplots(
    nrows=3, ncols=3, squeeze=True,
    #figsize=(8, 7),
    sharex=True, sharey=True,
)

for (i, j) in itertools.product(range(3), range(3)):
    ax = axes[i, j]
    ax.set_prop_cycle(cycler('color', palette))
    my_scaled_data = scaled_data[i][j]
    for l_index, l in enumerate(ls):
        ax.plot(
            my_scaled_data.x[l_index, :], my_scaled_data.y[l_index, :],
            '.',
            label=r'${}$'.format(l),
            rasterized=True,
        )
    ax.set_xbound(-5, 2)
    if i == 0:
        ax.set_title(
            r'$\rho_c = {}$'.format(rho_c[i, j]),
            position=(0.25, 0.65),
        )
    elif i == 1:
        ax.set_title(
            r'$\nu = {}$'.format(nu[i, j]),
            position=(0.25, 0.65),
        )
    elif i == 2:
        ax.set_title(
            r'$\zeta = {}$'.format(zeta[i, j]),
            position=(0.25, 0.65),
        )
    if i == 2:
        ax.set_xlabel(r'$x$')
        ax.set_xticks([-4, -2, 0, ])
    if j == 0:
        ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.text(
        0.1, 0.5,
        r'$S={:.1f}$'.format(quality[i][j]),
        transform=ax.transAxes,
    )

plt.show()


print("###############################################")
print("###############################################")
print("fss analysis for")
print("ls ", ls)  # wielkości systemu
print("rhos ", rhos)  # argumenty (osobno dla każdego L)
print("a ", a)  # chyba wartości funkcji (osobno dla każdego L)
print("da", da)  # chyba błędy standardowe funkcji (osobno dla każdego L i każdej wartości)
print("###############################################")
print("###############################################")
ret = fssa.autoscale(ls, rhos, a, da, 0.4, 1.8, 2.2)
print(ret)

auto_scaled_data = fssa.scaledata(ls, rhos, a, da, ret.rho, ret.nu, ret.zeta)
print(ret.rho, ret.drho)
print(ret.nu, ret.dnu)
print(ret.zeta, ret.dzeta)
print(ret.fun)

fig, ax = plt.subplots()
ax.set_prop_cycle(cycler('color', palette))
ax.plot(
    auto_scaled_data.x.T, auto_scaled_data.y.T,
    '.',
)
ax.set_xbound(-4, 2)
ax.set_xlabel(r'$x$')
plt.show()

# configure plotting
"""%config InlineBackend.rc = {'figure.dpi': 300, 'savefig.dpi': 300, \
                            'figure.figsize': (6, 6 / 1.6), 'font.size': 12, \
                            'figure.facecolor': (1, 1, 1, 0)}
%matplotlib inline"""


