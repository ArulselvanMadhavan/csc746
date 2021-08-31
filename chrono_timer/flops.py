import csv
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def calculate_flop(line_no, lines):
    count = 3
    found = 0
    result = []
    while (line_no < len(lines)):
        l = lines[line_no]
        if "Elapsed time" in l:
            el = l.split(":")[1]
            result.append(float(el.strip()))
            found += 1
        if count == found:
            return line_no+1, result
        line_no+=1

def infer_problem_size(line):
    ps = line.split("=")[1]
    ps = ps.strip()
    return int(ps)

def calculate_flops(lines):
    lines = lines[1:]
    line_no = 0
    final = []
    while (line_no < len(lines)):
        line = lines[line_no]
        if "Problem size" in line:
            ps = infer_problem_size(line)
            line_no, result = calculate_flop(line_no, lines)
            # print(ps, result)
            final.append([ps, result])
        line_no += 1
    return final

def read_lines(filename):
    with open(filename, 'r') as f:
        content = f.readlines()
    lines = [line.rstrip() for line in content]
    return lines

THEORETICAL_PEAK_FLOP = 44.8*math.pow(10,9)
THEORETICAL_MAX_BANDWIDTH = 102*math.pow(10,9)
SCIENTIFIC_NOTATION_FORMAT = "{:.1e}"

def calculate_peak_flop(flops):
    max_flop = max(flops)
    max_index = flops.index(max_flop)
    return max_index, max_flop/THEORETICAL_PEAK_FLOP

def number_of_accesses(problem_size):
    return [problem_size, 2*problem_size, 4*problem_size]
# Check if this is right
def bytes_accessed(problem_size):
    return [12*problem_size, (16*problem_size),(8+8+4+8)*problem_size]

def calculate_peak_bandwidth(bws):
    max_bw = max(bws)
    max_index = bws.index(max_bw)
    return max_index, max_bw/THEORETICAL_MAX_BANDWIDTH

def print_row_sub(xs):
    return " & ".join([SCIENTIFIC_NOTATION_FORMAT.format(x) for x in xs]) + " \\\\"

def interleave_entries(left, right):
    result = []
    for j in range(len(left)):
        result.append(left[j])
        result.append(right[j])
    return result

def run_through(left, right):
    for row in range(len(left)):
        l_row = left[row]
        r_row = right[row]
        ps, l_entries = l_row
        interleaved_entries = interleave_entries(l_entries, r_row[1])
        print("{}".format(int(ps/math.pow(10,5))) + " & " + print_row_sub(interleaved_entries))
        
def print_latex_row(result):
    for i in range(len(result)-1):
        f_0,b_0,l_0 = result[i]
        f_1,b_1,l_1 = result[i+1]
        run_through(f_0, f_1)   # flops
        print("**************")
        run_through(b_0, b_1)   # bw
        print("**************")
        run_through(l_0, l_1)   # latencies

def get_stats(filename):
    lines = read_lines(filename)
    result = (calculate_flops(lines))
    max_flop_percentage = []
    max_bw_percentage = []
    flops_rows = []
    bw_rows = []
    latency_rows = []
    elapsed_time = []
    for res in result:
        ps, ms = res[0],res[1]
        flops = []
        elapsed = []
        for el_time in ms:
            elapsed.append(el_time)
            flops.append(ps/el_time)
        print("FLOPS:", ps, ",".join([SCIENTIFIC_NOTATION_FORMAT.format(f) for f in flops]))
        flops_rows.append((ps, flops))
        elapsed_time.append((ps, elapsed))
        
        max_idx, max_flop = calculate_peak_flop(flops)
        max_flop_percentage.append(max_flop)

        bytes_acc = bytes_accessed(ps)
        bandwidth_used = [bytes_acc[i]/ms[i] for i in range(len(ms))]
        bw_rows.append((ps, bandwidth_used))
        print("Bandwidth used", ",".join([SCIENTIFIC_NOTATION_FORMAT.format(f) for f in bandwidth_used]))
        
        max_idx, max_bw = calculate_peak_bandwidth(bandwidth_used)
        max_bw_percentage.append(max_bw)
        accesses_count = number_of_accesses(ps)
        latency = [ms[i]/accesses_count[i] for i in range(len(ms))]
        latency_rows.append((ps, latency))
        print("Latency", ",".join([SCIENTIFIC_NOTATION_FORMAT.format(f) for f in latency]))

    return elapsed_time, flops_rows, bw_rows, latency_rows

def extract_bw(result):
    _,b,_ = result
    final = []
    for b_row in b:
        ps,entries = b_row
        final.append(entries)
    return final

def extract_l(result):
    _,_,l = result
    return [(number_of_accesses(ps), entries)for ps,entries in l]

def merge_rows(acc, l_r, r_r):
    row = []
    for i in range(1,len(acc)):
        row.append(acc[i])
        row.append(l_r[i])
        row.append(r_r[i])
    return row

def print_latency_rows(l_rows, r_rows):
    print(r_rows)
    for i in range(len(l_rows)):
        acc, l_row = l_rows[i]
        _, r_row = r_rows[i]
        row = merge_rows(acc, l_row, r_row)
        print(print_row_sub(row))
    
def print_acc_vs_latency(results):
    l_entries = [extract_l(result) for result in results]
    for i in range(len(l_entries)-1):
        l_rows = l_entries[i]   # O0
        r_rows = l_entries[i+1]   # O3
        print_latency_rows(l_rows, r_rows)

PLOT_FILENAMES = ["nm.csv", "sm.csv", "um.csv"]

def produce_plot_data(el_entries):
    csv_entries = [[],[],[]]
    for i in range(len(el_entries) - 1):
        left = el_entries[i]
        right = el_entries[i + 1]
        for row_i in range(len(left)):
            ps, l_r = left[row_i]
            ps, r_r = right[row_i]
            for m_i in range(len(l_r)):
                l_e = l_r[m_i]
                r_e = r_r[m_i]
                csv_entries[m_i].append((ps, l_e, r_e))
                
    for i in range(len(csv_entries)):
        with open(PLOT_FILENAMES[i], 'w+') as f:
            c = csv_entries[i]
            for row in c:
                writer = csv.writer(f, delimiter=',')
                writer.writerow(row)
    return csv_entries

PLOT_NAMES = ["No Memory accesses", "Structured Memory accesses", "Unstructured memory accesses"]

def produce_plots():
    for i in range(len(PLOT_FILENAMES)):
        pf = PLOT_FILENAMES[i]
        pn = PLOT_NAMES[i]
        df = pd.read_csv(pf, header=None)
        plt.title(pn)
        plt.xticks([i for i in range(len(df))], df[0])
        plt.plot(df[1], "r-+")
        plt.plot(df[2], "b-x")
        plt.yscale("linear")
        plt.grid(True)
        plt.xlabel("Problem Size")
        plt.ylabel("Runtime (in secs)")
        varNames = ["O0", "O3"]
        plt.legend(varNames, loc='best')
        plt.show()
        # print(xlocs)

        
def print_bw_stats(results):
    bw_entries = [extract_bw(result) for result in results]
    arr = np.asarray(bw_entries)
    print(arr.shape)
    print("Min")
    min_arr = np.min(arr, axis=1)    
    print(min_arr)
    print(print_row_sub(np.reshape(min_arr, (-1,))))
    print("Max")
    max_arr = np.max(arr, axis=1)
    print(max_arr)
    print(print_row_sub(np.reshape(max_arr, (-1,))))
    print("Average")
    avg_arr = np.average(arr, axis=1)    
    print(avg_arr)
    print(print_row_sub(np.reshape(avg_arr, (-1,))))
            
if __name__ == "__main__":
    filenames = ["no_opt.log", "opt.log"]
    elapsed = []
    result = []
    for filename in filenames:
        (e, f, b, l) = get_stats(filename)
        result.append((f,b,l))
        elapsed.append(e)

    # print_latex_row(result)
    # print_bw_stats(result)
    # print_acc_vs_latency(result)
    produce_plot_data(elapsed)
    produce_plots()
        # bytes accessed / elapsed time = bytes/sec = bandwidth used
        # What % of peak bandwidth are you using?
        # Elapsed time / # accesses = average latency
