import numpy
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math

OUTPUT_FILE = "output.txt"
P = [4,9,16,25,36,49,64,81]
# P = [4]
# serial_time = [168.24, 0.044, 154.41,16.55]

def save_pdf(pdf_fn):
    with PdfPages(pdf_fn+'.pdf') as pdf:
        pdf.savefig()

def find_times(lines, lineno):
    while(lineno < len(lines)):
        line = lines[lineno]
        if "Timing results" in line:
            st = parse_time(lines[lineno+1])
            ht = parse_time(lines[lineno+2])
            ot = parse_time(lines[lineno+3])
            gt = parse_time(lines[lineno+4])
            return lineno+4, st, ht, ot, gt
        lineno+=1
    return -1, None, None, None, None

def parse_time(line):
    return float(line.split(":")[1].split("(")[0].strip())

phases = ["scatter", "ghost", "sobel", "gather"]
modenames = ["row", "column", "tiled"]
modeOpt = ["-g 1", "-g 2", "-g 3"]

def read_file(fname, mode):
    times_array = [[],[],[],[]]
    with open(fname, 'r') as fn:
        lines = fn.readlines()
        lineno = 0
        while(lineno < len(lines)):
            line = lines[lineno]
            if "srun" in line and mode in line:
                lineno, st,ht,ot,gt = find_times(lines, lineno)
                times_array[0].append(st)
                times_array[1].append(ht)
                times_array[2].append(ot)
                times_array[3].append(gt)
            lineno+=1
    print(f"Total lines:{lineno}")
    return times_array

def get_serial_time(mode):
    with open(f"serial_{mode}.txt", 'r') as fn:
        lines = fn.readlines()
        lno, st, ht, ot, gt = find_times(lines, 0)
    return st,ht,ot,gt
        
def plot_speedup(mode, times_array):
    plt.title(mode.capitalize()+" decomposition")
    plt.xticks([i for i in range(len(P))],P)
    plt.yscale("linear")
    plt.grid(True)
    plt.xlabel("Parallel Processes Count")
    plt.ylabel("Speedups")
    pid = 0
    for phase in phases:
        t = times_array[pid]
        stime = get_serial_time(mode)
        # print(stime)
        sp = [serial_time[pid]/t[idx] for idx in range(len(P))]
        plt.plot(sp, linestyle='solid', marker='x')
        pid+=1
    plt.legend(phases)
    save_pdf(f"{mode}_speedup")
    plt.clf()

WIDTH = 7112
HEIGHT = 5146
TOTAL_PIXELS = 7112*5146

def plot_data_mov():
    plt.title("Data Movement")
    plt.xticks([i for i in range(len(P))],P)
    plt.yscale("linear")
    plt.grid(True)
    plt.xlabel("Parallel Processes Count")
    plt.ylabel("Data Moved")
    
    
def calculate_tiles(mode, n):
    if mode == "row":
        xd = WIDTH
        yd = (HEIGHT)/n
    elif mode == "column":
        xd = (WIDTH)/n
        yd = HEIGHT
    elif mode == "tiled":
        procs = int(math.sqrt(n))
        xd = WIDTH/procs
        yd = HEIGHT/procs
    return int(xd), int(yd)
        
def scatter_data(xd, yd):
    return TOTAL_PIXELS - (xd * yd)

def gather_data(xd, yd):
    return TOTAL_PIXELS - (xd * yd)

nhalo = 1
MPI_PROC_NULL = -1

def calc_gdata(rank, nx, ny, nX, nY):
    nleft = rank - 1 if (nx > 0) else MPI_PROC_NULL
    nrght = rank + 1 if (nx < nX - 1) else  MPI_PROC_NULL
    ntop = rank - nX if (ny > 0) else MPI_PROC_NULL
    nbot = rank + nX if (ny < nY - 1)  else MPI_PROC_NULL
    gy = [yd for np in [nleft, nrght] if np != MPI_PROC_NULL]
    gx = [xd for np in [ntop, nbot] if np != MPI_PROC_NULL]
    return sum(gx) + sum(gy)
    
    
def ghost_per_tile(mode, p, xd, yd):
    msgCount = 0
    gd = 0
    if mode == "row":
        for i in range(p):
            gd += calc_gdata(i, 0, i, 1, p)
    elif mode == "column":
        for i in range(p):
            gd += calc_gdata(i, i, 0, p, 1)
    elif mode == "tiled":
        px = int(math.sqrt(p))
        py = px
        rank = 0
        for j in range(py):
            for i in range(px):
                gd += calc_gdata(rank, i, j, px, py)
                rank += 1
    return gd

def calc_gmsg(rank, nx, ny, nX, nY):
    nleft = rank - 1 if (nx > 0) else MPI_PROC_NULL
    nrght = rank + 1 if (nx < nX - 1) else  MPI_PROC_NULL
    ntop = rank - nX if (ny > 0) else MPI_PROC_NULL
    nbot = rank + nX if (ny < nY - 1)  else MPI_PROC_NULL
    gy = [1 for np in [nleft, nrght] if np != MPI_PROC_NULL]
    gx = [1 for np in [ntop, nbot] if np != MPI_PROC_NULL]
    return sum(gx) + sum(gy)

def msgCount(mode, p):
    tmsgs = 0
    if mode == "row":
        for i in range(p):
            tmsgs += calc_gmsg(i, 0, i, 1, p)
    if mode == "column":
        for i in range(p):
            tmsgs += calc_gmsg(i, i, 0, p, 1)
    if mode == "tiled":
        # Ghost
        px = int(math.sqrt(p))
        py = px
        rank = 0
        for j in range(py):
            for i in range(px):
                tmsgs += calc_gmsg(rank, i, j, px, py)
                rank += 1
    tmsgs += p - 1              # Scatter
    tmsgs += p - 1              # Gather
    return tmsgs

def toMB(pixels):
    pbytes = 8 * pixels
    return pbytes/1000000

if __name__ == "__main__":
    mid = 0
    # for mode in modenames:
    #     t_array = read_file(OUTPUT_FILE, modeOpt[mid])
    #     plot_speedup(mode, t_array)
    #     mid += 1
    for p in P:
        data_row = []
        msg_row = []
        for mode in modenames:
            xd,yd = calculate_tiles(mode, p)
            sdata = toMB(scatter_data(xd, yd))
            tdata = toMB(gather_data(xd, yd))
            gdata = toMB(ghost_per_tile(mode, p, xd, yd))
            # print(f"Pixels:{8*TOTAL_PIXELS}\ttotalMB:{toMB(8*TOTAL_PIXELS)}SD:{sdata}\tTD:{tdata}\tGD:{gdata}")
            total = sdata + gdata + tdata
            data_row.append(total)
            msgs = msgCount(mode, p)
            msg_row.append(msgs)
        print(f"{p} & {msg_row[0]:.2f} & {data_row[0]:.2f} & {msg_row[1]:.2f} & {data_row[1]:.2f} & {msg_row[2]:.2f} & {data_row[2]:.2f}")
