import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from matplotlib.backends.backend_pdf import PdfPages

def save_pdf(pdf_fn):
    with PdfPages(pdf_fn+'.pdf') as pdf:
        pdf.savefig()
        
def sobel_cpu_chart():
    files = glob.glob("logs/cpu/sobel_cpu*.txt")
    nThreadsArr = []
    el_times = []
    for f in files:
        nThreads = f.split("_")[2].split(".")[0]
        nThreads = int(nThreads)
        nThreadsArr.append(nThreads)        
        with open(f, 'r') as fn:
            lines = fn.readlines()
            line = lines[1]
            el_time = float(line.split(" : ")[1])
            el_times.append(el_time)
    sorted_t_indices = sorted(range(len(nThreadsArr)), key=lambda k: nThreadsArr[k])
    plt.title("Sobel CPU Performance - Concurrency vs Elapsed Time")
    plt.xticks([i for i in range(len(nThreadsArr))],
               [nThreadsArr[t] for t in sorted_t_indices])
    plt.yscale("linear")
    plt.grid(True)
    plt.xlabel("Number of threads")
    plt.ylabel("Elapsed time")
    plt.plot([el_times[t] for t in sorted_t_indices],linestyle='solid', marker='x')
    # nThreadsArr.sort()
    print(nThreadsArr)
    print(el_times)
    save_pdf("sobel_cpu")
    plt.clf()

thread_blocks = ["1", "4", "16", "64", "256", "1024", "4096"] # x axis, 7 of them
threads_per_block = ['32', '64', '128', '256', '512', '1024'] # y axis, 6 of them

b_dict = {}
for idx, k in enumerate(thread_blocks):
    b_dict[k]=idx

t_dict = {}
for idx, k in enumerate(threads_per_block):
    t_dict[k]=idx

def get_runtime():
    files = glob.glob("logs/gpu/exec/sobel_gpu_*.txt")
    bs = []
    ts = []
    runtime = np.zeros((6,7))
    print(b_dict)
    print(t_dict)
    for f in files:
        xs = f.split("_")
        b = xs[2]
        t = xs[3].split(".")[0]
        bs.append(b)
        ts.append(t)
        with open(f, 'r') as fn:
            lines = fn.readlines()[3:]
            rt = float(lines[2].split(",")[2])
            runtime[t_dict[t]][b_dict[b]] = round(rt,1)
    return runtime

def get_efficiency():
    files = glob.glob("logs/gpu/eff/sobel_gpu_*.txt")
    bs = []
    ts = []
    runtime = np.zeros((6,7))
    print(b_dict)
    print(t_dict)
    for f in files:
        xs = f.split("_")
        b = xs[3]
        t = xs[4].split(".")[0]
        bs.append(b)
        ts.append(t)
        with open(f, 'r') as fn:
            lines = fn.readlines()
            line = lines[5]
            eff = line.split(",")[12][:-1]
            eff = float(eff)
            eff = round(eff, 1)
            runtime[t_dict[t]][b_dict[b]] = eff
    return runtime

color_schemes = ["coolwarm", "coolwarm_r"]

def sobel_gpu_chart():
    idx = 0
    for runtime in [get_runtime(), get_efficiency()]:
        fig, ax = plt.subplots()
        im = ax.imshow(runtime, cmap=color_schemes[idx])
        ax.set_xticks(np.arange(len(thread_blocks)))
        ax.set_yticks(np.arange(len(threads_per_block)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(thread_blocks)
        ax.set_yticklabels(threads_per_block)
        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(threads_per_block)): # y axis
            for j in range(len(thread_blocks)): # x axis
                text = ax.text(j, i, runtime[i, j],
                               ha="center", va="center", color="k")

        ax.set_title("Runtime on Cori at Varying Block Size and Number of Blocks")
        ax.set_ylabel('Threads per block')
        ax.set_xlabel('Block Sizes')
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        save_pdf(f"sobel_gpu_{idx}")
        plt.clf()
        idx+=1
    
if __name__ == "__main__":
    # sobel_cpu_chart()
    sobel_gpu_chart()
