import matplotlib.pyplot as plt
import os
from collections import defaultdict
import pandas as pd

MARKERS = ["o", "s", "^", "P", "p", "*", "h", "D", "d"]


def plot_dfr(xdata, ydata, label):
    plot_multi_dfr(xdata, [ydata], [label])


def plot_multi_dfr(xdata, ydatas, labels):
    plt.xlabel('Počet chýb')
    plt.ylabel('DFR')
    plt.yscale("log")
    plt.ylim(ymin=(10 ** (-4)), ymax=1)
    for i, ydata in enumerate(ydatas):
        plt.plot(xdata, ydata, marker=MARKERS[i % len(MARKERS)], label=labels[i])
    plt.legend()
    plt.show()


def plot_avg_syndrome_weight():
    data = [list() for _ in range(89)]  # max achieved iteration in the dataset was 87
    for textfile in os.listdir("syndrome"):
        with open(os.path.join("syndrome", textfile), "r") as f:
            line = f.readline()
            line = line.split(";")
            line = [int(v) for v in line]
            for i, val in enumerate(line):
                data[i].append(val)
    yvals = [sum(l)/len(l) for l in data]
    xvals = [i for i in range(len(data))]
    plt.xlabel('Iterácia')
    plt.ylabel('Priemerná váha syndrómu')
    plt.plot(xvals, yvals, marker="o")
    plt.show()


def plot_avg_max_sigma():
    sums = [0 for _ in range(88)]
    counts = [0 for _ in range(88)]

    for experiment in range(100):
        p = os.path.join("exp", f"exp{experiment}")
        for iteration in range(88):
            path = os.path.join(p, f"sigmas-exp_{experiment}-iter_{iteration}.txt")
            if os.path.exists(path):
                with open(path, "r") as f:
                    lines = f.readlines()
                lines = [line.replace("\n", "") for line in lines]
                lines = [line.replace(";", " ").strip() for line in lines]
                lines = [line.split() for line in lines]
                maximums = []
                for i in range(len(lines[0])):
                    a, b, c = lines[0][i], lines[1][i], lines[2][i]
                    maximums.append(max([int(a), int(b), int(c)]))
                sums[iteration] += max(maximums)
                counts[iteration] += 1
    avg = [float(s) / c for (s, c) in zip(sums, counts)]
    x = [i for i in range(88)]
    plt.xlabel('Iterácia')
    plt.ylabel('Priemerná hodnota $\sigma_{max}$')
    plt.plot(x, avg, marker="o")

    plt.show()


def plot_same_as_max():
    counts = defaultdict(list)
    for experiment in range(100):
        p = os.path.join("exp", f"exp{experiment}")
        for iteration in range(88):
            path = os.path.join(p, f"sigmas-exp_{experiment}-iter_{iteration}.txt")
            if not os.path.exists(path):
                continue
            with open(path, "r") as f:
                lines = f.readlines()
            ls = []
            for line in lines:
                line = line.replace("\n", "").replace(";", " ").strip().split(" ")
                line = [int(x) for x in line]
                ls.append(line)
            lines = ls
            maximums = []
            for a, b, c in zip(*lines):
                maximums.append(max([a, b, c]))
            m = max(maximums)
            count = 0
            for a, b, c in zip(*lines):
                if a == m or b == m or c == m:
                    count += 1
            counts[iteration].append(count)

    plt.figure(figsize=(30, 10))
    plt.xlabel('Iterácia')
    plt.ylabel('Počet pozícií dosahujúcich $\sigma_{max}$')
    plt.yticks([i for i in range(11)])
    plt.boxplot(counts.values(), labels=counts.keys(), notch=True)
    plt.show()


if __name__ == "__main__":
    """
    x = [88, 90, 92, 94, 96]
    # y = [2, 38, 232, 760, 1943]  # 100 100 2293 37 88-96 0
    ys = [[0.00001, 8, 56, 247, 816], [2, 7, 35, 229, 766]]
    ys = [[val / 10_000.0 for val in y] for y in ys]
    labels = ["k=2339, dec=0", "k=2339, dec=1"]
    # plot_dfr(x, y)
    plot_multi_dfr(x, ys, labels)
    """
    # plot_avg_syndrome_weight()
    # plot_avg_max_sigma()
    plot_same_as_max()