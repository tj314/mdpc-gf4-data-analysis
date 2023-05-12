import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os
from collections import defaultdict, Counter
import numpy as np

MARKERS = ["o", "s", "^", "P", "p", "*", "h", "D", "d"]


def plot_dfr(xdata, ydata, label, fname=None):
    plot_multi_dfr(xdata, [ydata], [label], fname)


def plot_multi_dfr(xdata, ydatas, labels, fname=None):
    plt.figure(figsize=(5, 5))
    plt.xlabel('Počet chýb')
    plt.ylabel('DFR')
    plt.yscale("log")
    plt.ylim(ymin=(10 ** (-4)), ymax=1)
    plt.xticks(ticks=xdata)
    plt.grid()
    plt.grid(which="minor", linestyle="--")
    for i, ydata in enumerate(ydatas):
        plt.plot(xdata, ydata, marker=MARKERS[i % len(MARKERS)], label=labels[i])
    plt.legend(loc='lower right')
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname + ".png", bbox_inches="tight")


def plot_avg_syndrome_weight():
    data = [list() for _ in range(89)]  # max achieved iteration in the dataset was 87
    for textfile in os.listdir("syndrome"):
        with open(os.path.join("syndrome", textfile), "r") as f:
            line = f.readline()
            line = line.split(";")
            line = [int(v) for v in line]
            for i, val in enumerate(line):
                data[i].append(val)
    yvals = [sum(l) / len(l) for l in data]
    xvals = [i for i in range(len(data))]
    plt.xlabel('Iterácia')
    plt.ylabel('Priemerná váha syndrómu')
    plt.plot(xvals, yvals, marker="o")
    plt.show()


def plot_same_as_max(fname=None):
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

    """
    plt.figure(figsize=(30, 10))
    plt.xlabel('Iterácia')
    plt.ylabel('Počet pozícií dosahujúcich $\sigma_{max}$')
    plt.yticks([i for i in range(11)])
    plt.boxplot(counts.values(), labels=counts.keys(), notch=True)
    """
    x = [i for i in sorted(counts.keys())]
    y = [float(sum(counts[i])) / len(counts[i]) for i in x]
    plt.xlabel("Iterácia")
    plt.ylabel("Priemerný počet |{j: $\sigma_{max} = \sigma_j$}|")
    plt.plot(x, y)
    if fname is not None:
        plt.savefig(fname + "avg.png", bbox_inches="tight")
    plt.figure()
    # y = [float(sum(counts[pos])) / len(counts[pos]) for pos in x]
    y1 = [len([v for v in counts[pos] if v == 1]) for pos in x]
    y2 = [len([v for v in counts[pos] if v == 2]) for pos in x]
    y3 = [len([v for v in counts[pos] if v == 3]) for pos in x]
    y4 = [len([v for v in counts[pos] if v == 4]) for pos in x]
    plt.xlabel("Iterácia")
    plt.ylabel("Percento dekódovaní")
    plt.ylim(ymin=0, ymax=100)
    plt.grid()
    plt.plot(x, y1, c="r", label="C = 1")
    plt.plot(x, y2, c="g", label="C = 2")
    plt.plot(x, y3, c="b", label="C = 3")
    plt.plot(x, y4, c="orange", label="C = 4")
    plt.legend(loc='upper right', ncols=4)
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname + "by_count.png", bbox_inches="tight")


def plot_histogram_of_sigmas_per_syndrome_weight(fname=None):
    syndrome_p = "syndrome"
    exp_p = "exp"
    correct_sigmas_per_syndrome: dict[int, list[int]] = defaultdict(list)
    incorrect_sigmas_per_syndrome: dict[int, list[int]] = defaultdict(list)
    all_syndromes = []
    print("First stage!")
    for experiment in range(100):
        print(f"Progress: {experiment}/100")
        syndrome_path = os.path.join(syndrome_p, f"syndrome-exp_{experiment}.txt")
        exp_path = os.path.join(exp_p, f"exp{experiment}")
        errvec_path = os.path.join(exp_path, f"errorvec-exp_{experiment}.txt")
        with open(syndrome_path, "r") as f:
            line = f.readline()
        syndromes = [int(x) for x in line.replace("\n", "").split(";")]
        all_syndromes.extend(syndromes)
        with open(errvec_path, "r") as f:
            errvec = f.readline()
        errvec = errvec.replace("\n", "").replace(";", " ").strip().split(" ")
        errvec = [int(val) for val in errvec]
        for iteration in range(88):
            path = os.path.join(exp_path, f"sigmas-exp_{experiment}-iter_{iteration}.txt")
            if not os.path.exists(path):
                # print(f"path = {path} does not exist!")
                continue
            with open(path, "r") as f:
                lines = f.readlines()
            ls = []
            for line in lines:
                line = line.replace("\n", "").replace(";", " ").strip().split(" ")
                line = [int(x) for x in line]
                ls.append(line)
            lines = ls

            syndrome_this_iteration = syndromes[iteration]
            for i in range(4678):  # 2*2339
                a, b, c = lines[0][i], lines[1][i], lines[2][i]
                tmp = max([a, b, c])
                if errvec[i] == 0:
                    incorrect_sigmas_per_syndrome[syndrome_this_iteration].append(tmp)
                else:
                    correct_sigmas_per_syndrome[syndrome_this_iteration].append(tmp)

    # find thresholds per syndrome weight
    all_syndromes = [v for v in sorted(all_syndromes, reverse=True)][:15]
    for syndrome in all_syndromes:
        plt.hist(incorrect_sigmas_per_syndrome[syndrome], color="red", alpha=0.5, bins="doane")
        plt.hist(correct_sigmas_per_syndrome[syndrome], color="blue", alpha=0.5, bins="doane")
        plt.xlabel('Hodnota $\sigma_j$ pri w(s)=' + str(syndrome))
        plt.ylim(ymin=0, ymax=150)
        plt.ylabel('Počet výskytov')
        if fname is None:
            print(f"syndrome weight: {syndrome}")
            plt.show()
        else:
            plt.savefig(f"{fname}-{syndrome}.png", bbox_inches="tight")
            plt.clf()


def find_threshold_per_syndrome_weight_all(fname=None):
    syndrome_p = "syndrome"
    exp_p = "exp"
    correct_sigmas_per_syndrome: dict[int, list[int]] = defaultdict(list)
    incorrect_sigmas_per_syndrome: dict[int, list[int]] = defaultdict(list)
    all_syndromes = []
    print("First stage!")
    for experiment in range(100):
        print(f"Progress: {experiment}/100")
        syndrome_path = os.path.join(syndrome_p, f"syndrome-exp_{experiment}.txt")
        exp_path = os.path.join("exp", f"exp{experiment}")
        errvec_path = os.path.join(exp_path, f"errorvec-exp_{experiment}.txt")
        with open(syndrome_path, "r") as f:
            line = f.readline()
        syndromes = [int(x) for x in line.replace("\n", "").split(";")]
        all_syndromes.extend(syndromes)
        with open(errvec_path, "r") as f:
            errvec = f.readline()
        errvec = errvec.replace("\n", "").replace(";", " ").strip().split(" ")
        errvec = [int(val) for val in errvec]
        for iteration in range(88):
            path = os.path.join(exp_path, f"sigmas-exp_{experiment}-iter_{iteration}.txt")
            if not os.path.exists(path):
                # print(f"path = {path} does not exist!")
                continue
            with open(path, "r") as f:
                lines = f.readlines()
            ls = []
            for line in lines:
                line = line.replace("\n", "").replace(";", " ").strip().split(" ")
                line = [int(x) for x in line]
                ls.append(line)
            lines = ls

            syndrome_this_iteration = syndromes[iteration]
            for i in range(4678):  # 2*2339
                a, b, c = lines[0][i], lines[1][i], lines[2][i]
                tmp = max([a, b, c])
                if errvec[i] == 0:
                    incorrect_sigmas_per_syndrome[syndrome_this_iteration].append(tmp)
                else:
                    correct_sigmas_per_syndrome[syndrome_this_iteration].append(tmp)

    # find thresholds per syndrome weight
    all_syndromes = [v for v in sorted(all_syndromes, reverse=True)][:8000]
    ydata = []
    pos = 0
    print("Second stage!")
    length = len(all_syndromes)
    for syndrome in all_syndromes[:]:
        if pos % 100 == 0:
            print(f"Progress: {pos}/{length}")
        pos += 1
        all_sigmas_per_syndrome = set(incorrect_sigmas_per_syndrome[syndrome])
        for sigma in correct_sigmas_per_syndrome[syndrome]:
            all_sigmas_per_syndrome.add(sigma)
        all_sigmas_per_syndrome = sorted(list(all_sigmas_per_syndrome))
        num_correct_left = 0
        num_correct_total = len(correct_sigmas_per_syndrome[syndrome])
        num_incorrect_left = 0
        num_incorrect_total = len(incorrect_sigmas_per_syndrome[syndrome])
        correct_counts = Counter(correct_sigmas_per_syndrome[syndrome])
        incorrect_counts = Counter(incorrect_sigmas_per_syndrome[syndrome])
        bound = 0.99
        while True:
            if bound < 0.0:
                all_syndromes.remove(syndrome)
                break
            found = False
            for sigma in all_sigmas_per_syndrome:
                if float(num_incorrect_left) >= bound*num_incorrect_total \
                        and float(num_correct_left) < (1-bound)*num_correct_total:
                    ydata.append(sigma)
                    found = True
                    break
                else:
                    num_incorrect_left += incorrect_counts[sigma]
                    num_incorrect_left += correct_counts[sigma]
            if found:
                break
            else:
                bound -= 0.01

    # approximate
    model = np.polyfit(all_syndromes, ydata, 1)
    print(model)  # result coefficients
    predict = np.poly1d(model)
    predict_y = predict(all_syndromes)

    #plt.figure(figsize=(20, 10))
    plt.scatter(all_syndromes, ydata)
    plt.plot(all_syndromes, predict_y, color="red")
    plt.ylim(ymin=-20, ymax=20)
    plt.xlim(xmin=800, xmax=all_syndromes[0]+1)
    plt.ylabel("Hodnota $\sigma$")
    plt.xlabel("Váha syndrómu")
    plt.gca().xaxis.set_major_locator(MultipleLocator(100))
    plt.gca().yaxis.set_minor_locator(MultipleLocator(1))
    plt.gca().yaxis.set_major_locator(MultipleLocator(5))
    plt.gca().yaxis.grid()
    plt.gca().invert_xaxis()
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname + ".png", bbox_inches="tight")

def plot_iterations_delta(fname=None):
    iterations_per_decoder = []
    filenames = [
        "iterations/iteracie-dec_0-delta_0.txt",
        "iterations/iteracie-dec_2-delta_0.txt",
        "iterations/iteracie-dec_2-delta_1.txt",
        "iterations/iteracie-dec_2-delta_2.txt",
        "iterations/iteracie-dec_2-delta_3.txt",
        "iterations/iteracie-dec_2-delta_4.txt",
        "iterations/iteracie-dec_2-delta_5.txt"
    ]
    decoder_names = [
        "SF",
        "$\delta$=0",
        "$\delta$=1",
        "$\delta$=2",
        "$\delta$=3",
        "$\delta$=4",
        "$\delta$=5"
    ]
    for filename in filenames:
        with open(filename, "r") as f:
            line = f.readline()
        line = line.replace("\n", "").replace(";", " ").strip().split(" ")
        line = [int(x) for x in line]
        iterations_per_decoder.append(line)
    plt.figure(figsize=(10, 15))
    plt.yticks([i * 2 for i in range(51)])
    plt.ylim(ymin=0, ymax=100)
    plt.grid(axis="y", alpha=0.4)
    plt.xlabel('Dekodér')
    plt.ylabel('Počet iterácií')
    plt.boxplot(iterations_per_decoder, labels=decoder_names, notch=True)
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname + ".png", bbox_inches="tight")


def plot_iterations_threshold(fname=None):
    iterations_per_decoder = []
    filenames = [
        "iterations/iteracie-dec_0-delta_0.txt",
        "iterations/iteracie-dec_3-threshold_1.txt",
        "iterations/iteracie-dec_3-threshold_2.txt",
        "iterations/iteracie-dec_3-threshold_3.txt",
        "iterations/iteracie-dec_3-threshold_4.txt",
    ]
    decoder_names = [
        "SF",
        "$T_1(s)$",
        "$T_2(s)$",
        "$T_3(s)$",
        "$T_4(s)$"
    ]
    for filename in filenames:
        with open(filename, "r") as f:
            line = f.readline()
        line = line.replace("\n", "").replace(";", " ").strip().split(" ")
        line = [int(x) for x in line]
        iterations_per_decoder.append(line)
    plt.figure(figsize=(10, 15))
    plt.yticks([i * 2 for i in range(51)])
    plt.ylim(ymin=0, ymax=100)
    plt.grid(axis="y", alpha=0.4)
    plt.xlabel('Dekodér')
    plt.ylabel('Počet iterácií')
    plt.boxplot(iterations_per_decoder, labels=decoder_names, notch=True)
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname + ".png", bbox_inches="tight")


if __name__ == "__main__":
    x = [88, 90, 92, 94, 96]
    """
    y1 = [3, 38, 232, 760, 1943]  # 100 100 2293 37 88-96 0
    y2 = [2, 15, 91, 476, 1582]   # 200 100 2339 37 88-96 0
    y1 = [float(val)/10000 for val in y1]
    y2 = [float(val)/20000 for val in y2]
    ys = [y1, y2]
    labels = ["k=2293", "k=2339"]
    """
    """
    y = [3 / 10000.0, 38 / 10000.0, 232 / 10000.0, 760 / 10000.0, 1943 / 10000.0]  # 100 100 2293 37 88-96 0
    y0 = [10 / 20000.0, 47 / 10000.0, 228 / 10000.0, 732 / 10000.0, 2028 / 10000.0]  # delta=0
    y1 = [17 / 20000.0, 50 / 10000.0, 316 / 10000.0, 1020 / 10000.0, 2538 / 10000.0]  # delta=1
    y2 = [24 / 20000.0, 78 / 10000.0, 408 / 10000.0, 1354 / 10000.0, 3308 / 10000.0]  # delta=2
    y3 = [27 / 20000.0, 119 / 10000.0, 637 / 10000.0, 1995 / 10000.0, 4175 / 10000.0]  # delta=3
    y4 = [95 / 20000.0, 261 / 10000.0, 1102 / 10000.0, 3183 / 10000.0, 5693 / 10000.0]  # delta=4
    y5 = [300 / 20000.0, 757 / 10000.0, 2254 / 10000.0, 4793 / 10000.0, 7293 / 10000.0]  # delta=5
    ys = [y, y0, y1, y2, y3, y4, y5]
    labels = ["SF", "$\delta=0$", "$\delta=1$", "$\delta=2$", "$\delta=3$", "$\delta=4$", "$\delta=5$"]
    plot_multi_dfr(x, ys, labels, "sf-delta-dfr")
    """

    """
    y0 = [3 / 10000.0, 38 / 10000.0, 232 / 10000.0, 760 / 10000.0, 1943 / 10000.0]  # 100 100 2293 37 88-96 0
    y1 = [2538 / 20_000.0, 2930 / 10_000.0, 5385 / 10_000.0, 7693 / 10_000.0, 9068 / 10_000.0]
    y2 = [178 / 20_000.0, 493 / 10_000.0, 1554 / 10_000.0, 3496 / 10_000.0, 5947 / 10_000.0]
    y3 = [45 / 20_000.0, 154 / 10_000.0, 601 / 10_000.0, 1712 / 10_000.0, 3533 / 10_000.0]
    y4 = [39 / 20_000.0, 142 / 10_000.0, 506 / 10_000.0, 1819 / 10_000.0, 3881 / 10_000.0]
    ys = [y0, y1, y2, y3, y4]
    labels = ["SF", "$T_1(s)$", "$T_2(s)$", "$T_3(s)$", "$T_4(s)$"]
    plot_multi_dfr(x, ys, labels, "sf-threshold-dfr")
    """

    # plot_avg_syndrome_weight()
    # plot_same_as_max("sf-count-same-as-sigma-max")
    find_threshold_per_syndrome_weight_all("sf-threshold-approx")
    # plot_histogram_of_sigmas_per_syndrome_weight("sf-histogram-sigmas")
    # plot_iterations_delta("sf-delta-iterations")
    # plot_iterations_threshold("sf-threshold-iterations")
