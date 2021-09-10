import pandas as pd
import csv
import matplotlib.pyplot as plt

def read_lines(filename):
    with open(filename, 'r') as f:
        content = f.readlines()
    lines = [line.rstrip() for line in content]
    return lines

block_sizes = [2,4,16,32,64]

def infer_ps(lines):
    line_no = 0
    pss = []
    els = []
    refs = []
    while line_no < len(lines):
        if "N=" in lines[line_no]:
            ps = lines[line_no].split("N=")[1]
            pss.append(ps)
            line_no+=1
            if "Blocked" in lines[line_no]: # Jackpot
                el_b = []
                re_b = []
                line_no+=1
                for b in block_sizes:
                    el_time = lines[line_no].split(":")[1].strip()
                    line_no+=1
                    ref_time= lines[line_no].split(":")[1].strip()
                    line_no+=1
                    line_no+=1
                    el_b.append(el_time)
                    re_b.append(ref_time)
                els.append(el_b)
                refs.append(re_b)
                line_no-=1
            else:
                el_time = lines[line_no].split(":")[1].strip()
                line_no+=1
                ref_time= lines[line_no].split(":")[1].strip()
                line_no+=1
                els.append(el_time)
                refs.append(ref_time)
        line_no+=1
    print(pss)
    print(els)
    print(refs)
    return pss, els, refs

def build_csv(fn, pss, els, refs):
    with open(fn, 'w+') as f:
        for i in range(len(pss)):
            writer = csv.writer(f, delimiter=',')
            if isinstance(els[i],list):
                row = []
                row.append(pss[i])
                for idx in range(len(block_sizes)):
                    row.append(els[i][idx])
                row.append(refs[i][0])
                    # b = block_sizes[idx]
                writer.writerow(row)
            else:
                writer.writerow((pss[i], els[i], refs[i]))      
        
PLOT_FILENAMES = ["basic.csv", "basic-copy.csv", "blocked.csv"]
PLOT_COLORS = ['k','c','m','y','g','b','r']

def produce_plots(fn):
    df = pd.read_csv(fn, header=None)
    plt.title(fn)
    plt.xticks([i for i in range(len(df))], df[0])
    # plt.plot(df[1], "r-+")
    # plt.plot(df[2], "b-x")
    rows, cols = df.shape
    for col in range(1,cols):
        plt.plot(df[col], PLOT_COLORS[col])
    plt.yscale("linear")
    plt.grid(True)
    plt.xlabel("Problem Size")
    plt.ylabel("Runtime (in secs)")
    if cols > 3:
        varNames = block_sizes + ["CBLAS"]
    else:
        varNames = ["Custom Impl", "CBLAS"]
    plt.legend(varNames, loc='best')
    plt.show()

def produce_tables(fn):
    df = pd.read_csv(fn, header=None)
    rows, cols = df.shape
    result = []
    for idx, row in df.iterrows():
        rs = " & ".join([str(el) for el in row])
        rs = rs + " \\\\"
        result.append(rs)
    print("\n".join(result))
    
if __name__ == "__main__":
    filenames = ["basic.txt", "basic-copy.txt","blocked.txt"]
    for idx in range(len(filenames)):
        filename = filenames[idx]
        lines = read_lines(filename)
        lines = lines[2:]
        pss,els, refs = infer_ps(lines)
        build_csv(PLOT_FILENAMES[idx], pss, els, refs)
        produce_plots(PLOT_FILENAMES[idx])
        produce_tables(PLOT_FILENAMES[idx])
