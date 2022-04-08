import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy import double

epochs, n, dim, lb, ub = map(double, input().split(" "))
epochs = int(epochs)
n = int(n)
dim = int(dim)

x = []
y = []
fig, ax = plt.subplots()

def animate(frame):
    ax.clear()
    # ax.plot(700,700,"bo") // For food
    for i in range(n):
        x, y = map(double, input().split())
        ax.plot(x, y, "ro")
        ax.set_xlim([lb,ub])
        ax.set_ylim([lb,ub])

ani = FuncAnimation(fig, animate, frames=epochs, interval=500, repeat=False)

ani.save("video.mp4")
