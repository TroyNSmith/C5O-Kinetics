import matplotlib.pyplot as plt


def visualize_scan_energy(steps: list, energies: list):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(steps, energies, marker="o", linestyle="None")
    for x, y in zip(steps, energies):
        ax.text(x, y, f"{x}", fontsize=8)
    ax.set_xlabel("Step")
    ax.set_ylabel("Energy")
    ax.set_xticks([2 * i for i in range(0, 15)])
    plt.show()
