import csv
import pandas as pd
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

N = [128, 512, 2048]
B = [4, 16]
T = [1,4,16,64]

MODES=["basic"]
       # "blas", "blocked"]

def get_runtime(df):
    return float(df.iloc[[5]][1])

def get_key(ps, bs):
    return f"{ps}*{bs}"

def save_pdf(pdf_fn):
    with PdfPages(pdf_fn+'.pdf') as pdf:
        pdf.savefig()

def plot(mode, xy):
    plt.title(mode)
    ps = {}
    for k in xy.keys():
        ps[k.split('*')[0]] = True
    plt.xticks([i for i in range(len(ps))],[int(i) for i in ps.keys()])
    plt.plot(xy.values(),linestyle='solid', marker='x')
    plt.yscale("linear")
    plt.grid(True)
    plt.xlabel("Problem size")
    plt.ylabel("Speedups")
    plt.legend(xy.keys())
    save_pdf(mode)
    plt.clf()
# S/T(n,p)
def speedup(mode_data):
    # Find serial time
    serial_ts = {}
    for i in range(len(mode_data[0])):
        ps,bs,cs,dfs = mode_data[0],mode_data[1],mode_data[2],mode_data[3]
        if cs[i] == "0":
            serial_ts[get_key(ps[i],bs[i])] = get_runtime(dfs[i])

    speedups = {}         
    for i in range(len(mode_data[0])):
        ps,bs,cs,dfs = mode_data[0],mode_data[1],mode_data[2],mode_data[3]
        if cs[i] == "0":
            continue
        speedups[ps[i]] = serial_ts[get_key(ps[i], bs[i])]/get_runtime(dfs[i])
    return speedups
if __name__ == "__main__":
    for mode in MODES:
        files = glob.glob("{}-*.csv".format(mode))
        mode_data = [[],[],[],[]]
        for f in files:
            print(f"Parsing-{f}")
            with open(f, 'r') as fh:
                df = pd.read_csv(fh, header=None, skiprows=[4,56])
                splits = f.split(".csv")[0].split("-")
                mode_data[0].append(splits[1])
                mode_data[1].append(splits[2])
                mode_data[2].append(splits[3])
                mode_data[3].append(df)
            print(f, df.shape)
        speedups = speedup(mode_data)
        plot(mode, speedups)
