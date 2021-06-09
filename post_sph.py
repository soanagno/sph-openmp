import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

head = 1; skip = 1; fps = 30
pathx = "x.txt"
pathy = "y.txt"
x = np.loadtxt(pathx, dtype=float)
y = np.loadtxt(pathy, dtype=float)
# print('Data shape:', x.shape)

# Parameters
dt = x[0, 0]; steps = int(x[0, 1]); t_total = x[0, 2]
Dx = x[0, 3]; Dy = x[0, 4]; dx = x[0, 5]
bt = int(x[0, 6]); L = x[0, 7]; H = x[0, 8]; n = int(x[0, 9]); nb = int(x[0, 10]); n1 = int(x[0, 11])
rect = patches.Rectangle((0, 0), L+2*Dx, H+2*Dy, linewidth=1, edgecolor='r', facecolor='none')

# Plotting
fig, ax = plt.subplots()
major_ticks = np.arange(0, L+2*Dx, Dx)
point_size = int(10000/n)
sz = 5
def plotSPH(i):
    # plt.scatter(x[i], y[i], s=5)
    plt.scatter(x[i, :n1], y[i, :n1], s=sz, c='slateblue')
    plt.scatter(x[i, n1:n], y[i, n1:n], s=sz, c='deeppink')
    plt.scatter(x[i, n:], y[i, n:], s=sz, c='midnightblue')
    plt.ylim([0, 2+L+2*bt*dx])
    plt.xlim([0, 2+H+2*bt*dx])
    # ax.add_patch(rect)
    plt.axis('equal')
    # ax.set_xticks(major_ticks)
    # ax.set_yticks(major_ticks)
    # plt.grid(which='major')

# Animation
plt.ion()
for i in range(head, steps+head, skip):
    try:
        plotSPH(i)
        plt.draw()
        # plt.pause(0.5)
        plt.pause(0.01*1/fps)
        plt.clf()
    except Exception:
        break

plt.ioff()
plotSPH(i-skip)
plt.show()
