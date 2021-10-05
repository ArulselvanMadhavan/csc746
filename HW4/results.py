import csv
import pandas as pd
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

N = [128, 512, 2048]
B = [4, 16]
T = [1,4,16,64]
LINE_NO = {"0":59, "default":126}
MODES=["blas", "basic", "blocked"]

def get_runtime(cs, df):
    ln = -1
    cn = -1
    if LINE_NO.get(cs, None) is None:
        ln,cn = LINE_NO["default"], 3
    else:
        ln,cn = LINE_NO[cs],1
    return float(df.iloc[[ln - 3]][cn])

def get_key(ps, bs, cs):
    return f"{ps}*{bs}*{cs}"

def save_pdf(pdf_fn):
    with PdfPages(pdf_fn+'.pdf') as pdf:
        pdf.savefig()

def plot(mode, stimes, ptimes, xy):
    plt.title(mode)
    ps,bs = {},{}
    for kv in sorted(xy.items(),key=lambda kv:[int(k) for k in kv[0].split('*')]):
        print("Key", kv)
        k = kv[0]
        p,b,c = k.split('*')
        ps[p] = True
        if p == b:
            b = 'all'
        cs = bs.get(b, {})
        cps = cs.get(c,{})
        cps[p] = xy[k]
        cs[c] = cps
        bs[b] = cs
    print("PS", ps)
    print("BS-CS", bs)
    plt.xticks([i for i in range(len(ps))],[int(i) for i in ps.keys()])
    plt.yscale("linear")
    plt.grid(True)
    plt.xlabel("Problem size")
    plt.ylabel("Speedups")

    legend = []
    for b in bs.keys():
        for ck, cv in sorted(bs[b].items(), key=lambda kv:int(kv[0])):
            ck_i = int(ck) + 1
            if b != "all":
                lbl = f"{b}*{ck_i}"
            else:
                lbl = f"{ck_i}"
            legend.append(lbl)
            plotpts_t = sorted(cv.items(), key=lambda kv:int(kv[0]))
            print("Sorted Plot pts", plotpts_t)
            plotpts = [rt for p,rt in plotpts_t]
            # ptims = []
            # for p, rt in plotpts_t:
            #     bkey = p if b == 'all' else b
            #     ptims.append(ptimes[f"{p}*{bkey}*{ck}"])
            # print("Ptimes", ptims)
            plt.plot(plotpts,linestyle='solid', marker='x')
    plt.legend(legend)
    save_pdf(mode)
    plt.clf()

ROW_NAMES = ["CBLAS", "BasicMM-omp","BMMCO-omp"]

def get_ddr_read(cs, df):
    ln = -1
    cn = -1
    if cs == "0":
        ln, cn = 67, 1
    else:
        ln, cn = 134, 1
    return float(df.iloc[[ln]][cn])

def get_mcdram_read(cs, df):
    # 128,1
    ln, cn = -1, -1
    if cs == "0":
        ln, cn = 60,1
    else:
        ln, cn = 127,1
    return float(df.iloc[[ln]][cn])

def tables(idx, mode_data):
    
    for i in range(len(mode_data[0])):
        ps,bs,cs,dfs = mode_data[0],mode_data[1],mode_data[2],mode_data[3]
        if cs[i] == "0" and ps[i] == "2048":
            if ps[i] == bs[i]:
                rname = f"{ROW_NAMES[idx]}"
            else:
                rname = f"{ROW_NAMES[idx]}(b={bs[i]})"
            runtime = get_runtime(cs[i], dfs[i])
            mcdram = get_mcdram_read(cs[i], dfs[i])
            ddr_read = get_ddr_read(cs[i],dfs[i])            
            print(f"{rname} & {runtime} & {mcdram} & {ddr_read} & \\\\")
            
def table2(idx, mode_data):
    for i in range(len(mode_data[0])):
        ps,bs,cs,dfs = mode_data[0],mode_data[1],mode_data[2],mode_data[3]    
        if cs[i] == "63" and ps[i] == "2048":
            if ps[i] == bs[i]:
                rname = f"{ROW_NAMES[idx]}"
            else:
                rname = f"{ROW_NAMES[idx]}(b={bs[i]})"
            runtime = get_runtime(cs[i], dfs[i])
            mcdram = get_mcdram_read(cs[i], dfs[i])
            ddr_read = get_ddr_read(cs[i],dfs[i])            
            print(f"{rname} & {runtime} & {mcdram} & {ddr_read} & \\\\")        
        
            
def speedup(mode_data):
    # Find serial time
    serial_ts = {}
    for i in range(len(mode_data[0])):
        ps,bs,cs,dfs = mode_data[0],mode_data[1],mode_data[2],mode_data[3]
        if cs[i] == "0":
            serial_ts[get_key(ps[i],bs[i],cs[i])] = get_runtime(cs[i], dfs[i])
    print("Serial", serial_ts)
    speedups = {}
    ptimes = {}
    stimes = {}
    for i in range(len(mode_data[0])):
        ps,bs,cs,dfs = mode_data[0],mode_data[1],mode_data[2],mode_data[3]
        if cs[i] == "0":
            continue
        ptime = get_runtime(cs[i], dfs[i])
        stime = serial_ts[get_key(ps[i], bs[i], "0")]
        key = get_key(ps[i], bs[i], cs[i])
        speedups[key] = stime / ptime
        stimes[key] = stime
        ptimes[key] = ptime
    return stimes, ptimes, speedups

if __name__ == "__main__":
    for idx in range(0, len(MODES)):
        mode = MODES[idx]
        files = glob.glob("{}-*.csv".format(mode))
        mode_data = [[],[],[],[]]
        for f in files:
            # print(f"Parsing-{f}")
            with open(f, 'r') as fh:
                df = pd.read_csv(fh, header=None, skiprows=[4,56])
                splits = f.split(".csv")[0].split("-")
                mode_data[0].append(splits[1])
                mode_data[1].append(splits[2])
                mode_data[2].append(splits[3])
                mode_data[3].append(df)
        # stimes, ptimes, speedups = speedup(mode_data)
        # plot(mode, stimes, ptimes, speedups)
        tables(idx, mode_data)
    
        # table2(idx, mode_data)
