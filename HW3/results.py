import pandas as pd
import csv
import matplotlib.pyplot as plt
import math
from matplotlib.backends.backend_pdf import PdfPages

def read_lines(filename):
    with open(filename, 'r') as f:
        content = f.readlines()
    lines = [line.rstrip() for line in content]
    return lines

def infer_threads(line):
    ps,threads,schedule = line.split("N=")[1].split("\t")
    ts = threads.split("=")[1]
    ss = schedule.split("=")[1]
    return ps, ts, ss

def write_csv(idx, schedules):
    plot_fn = PLOT_FILENAMES[idx]
    for s_k,s_v in schedules.items():        
        for t_k in sorted(s_v):
            fn = plot_fn + '_' + s_k + '_' + t_k +'.csv'
            with open(fn, 'w+') as f:
                writer = csv.writer(f, delimiter=',')
                ps, es, rs = s_v[t_k]
                for p_i in range(len(ps)):
                    ps_float = float(ps[p_i])
                    es_float = float(es[p_i])
                    rs_float = float(rs[p_i])
                    writer.writerow((ps_float, es_float, rs_float))

def total_flops(n):
    return 2 * n * n

def mflops(ps_col, val_col):
    no_zeros = [1.0 if val_col[ix] == 0.0 else val_col[ix] for ix in range(len(val_col))]
    return [int(total_flops(ps_col[ix]) / (1000000*no_zeros[ix])) for ix in range(len(val_col))]
def save_pdf(pdf_fn):
    with PdfPages(pdf_fn+'.pdf') as pdf:
        pdf.savefig()

PLOT_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
PROBLEM_SIZES = [1024,2048,4096,8192,16384]

def speedup(thread_count, serial_col, parallel_col):
    def calc(idx,denom):
        num = serial_col[1][idx]
        return num/thread_count if denom == 0 else num/denom
    return [calc(idx, parallel_col[idx]) for idx in range(len(serial_col))]

def produce_flops(idx):
    plot_fn = PLOT_FILENAMES[idx]
    if plot_fn == 'basic':
        df = pd.read_csv(plot_fn + '_(null)_(null).csv', header=None)
        plt.title(plot_fn)
        plt.xticks([i for i in range(len(df))], [int(i) for i in df[0]])
        plt.plot(mflops(df[0],df[1]), linestyle='solid',color='b',marker='x')
        plt.plot(mflops(df[0],df[2]), linestyle='dashed',color='r',marker='+')    
        plt.yscale("linear")
        plt.grid(True)
        plt.xlabel("Problem size")
        plt.ylabel("MFLOPS")
        plt.legend(["Custom Impl", "CBLAS"], loc='best')
        save_pdf(plot_fn)
        plt.clf()
    else:
        schedules = ['static','dynamic']
        serial_df = pd.read_csv("basic" + '_(null)_(null).csv', header=None)
        for sch in schedules:
            plot_s_fn = plot_fn + '_' + sch
            plt.title(plot_s_fn)
            plt.xticks([i for i in range(len(PROBLEM_SIZES))], PROBLEM_SIZES)
            for t_idx in range(len(THREADS)):
                t = THREADS[t_idx]
                color = PLOT_COLORS[t_idx]
                plot_s_t_fn =  plot_s_fn + '_' + str(t)
                df = pd.read_csv(plot_s_t_fn + '.csv', header=None)
                speedups = speedup(t, serial_df,df[1])
                print("Thread\t speedup", t, speedups)
                plt.plot(speedups, linestyle='solid', marker='x', color=color)
                # plt.plot(mflops(df[0],df[2]), linestyle='dashed',color=color,marker='+')
            plt.yscale("linear")
            plt.grid(True)
            plt.xlabel("Problem size")
            plt.ylabel("Speedup")
            plt.legend(THREADS, loc='best')
            save_pdf(plot_s_fn)
            plt.clf()

# THREADS = [1,2,4,8,16,32,64,272]
THREADS=[1,2,4,8,16,32,64,272]

def produce_best_plot():
    best_df = pd.read_csv("openmp_static_272.csv", header=None)
    serial_df = pd.read_csv("basic_(null)_(null).csv", header=None)
    plt.title("Best plot - omp-static-272")
    plt.xticks([i for i in range(len(PROBLEM_SIZES))], PROBLEM_SIZES)
    plt.plot(mflops(best_df[0], best_df[1]), linestyle='solid', marker='x', color='b')
    plt.plot(mflops(serial_df[0], serial_df[2]), linestyle='dashed',marker='+', color='r')
    plt.yscale("linear")
    plt.grid(True)
    plt.xlabel("Problem size")
    plt.ylabel("MFLOPS")
    plt.legend(["omp-static-272", "serial CBLAS"], loc='best')
    save_pdf("omp-static-272")
    plt.clf()    
    
def infer_ps(lines):
    line_no = 0
    schedules = {}
    while line_no < len(lines):
        if "N=" in lines[line_no]:
            ps,ts, ss = infer_threads(lines[line_no])
            if schedules.get(ss) is None:
                schedules[ss] = {}
            if schedules.get(ss).get(ts) is None:
                schedules[ss][ts] = ([],[],[])
            line_no+=1
            el_time = lines[line_no].split(":")[1].strip()
            line_no+=1
            ref_time= lines[line_no].split(":")[1].strip()

            pss, els, refs = schedules[ss][ts]
            pss.append(ps)
            els.append(el_time)
            refs.append(ref_time)
            schedules[ss][ts] = (pss, els, refs)
        line_no+=1
    return schedules

PLOT_FILENAMES = ["basic", "openmp"]

if __name__ == "__main__":
    filenames = ["basic.txt","openmp.txt"]
    for idx in range(len(filenames)):
        fname = filenames[idx]
        lines = read_lines(fname)
        lines = lines[2:]
        schedules = infer_ps(lines)
        write_csv(idx, schedules)
        produce_flops(idx)
        produce_best_plot()
