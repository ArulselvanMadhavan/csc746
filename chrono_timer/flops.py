import csv
import math

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

def calculate_peak_flop(flops):
    max_flop = max(flops)
    max_index = flops.index(max_flop)
    return max_index, max_flop/THEORETICAL_PEAK_FLOP

def bytes_accessed(problem_size):
    return [16, (8)+(8*problem_size),8+(2*8*problem_size)]

def calculate_peak_bandwidth(bws):
    max_bw = max(bws)
    max_index = bws.index(max_bw)
    return max_index, max_bw/THEORETICAL_MAX_BANDWIDTH

if __name__ == "__main__":
    filename = "no_opt.log"
    lines = read_lines(filename)
    result = (calculate_flops(lines))
    max_flop_percentage = []
    max_bw_percentage = []
    for res in result:
        ps, ms = res[0],res[1]
        flops = []
        for el_time in ms:
            flops.append(ps/el_time)
        print("FLOPS:", ps, ",".join(["{:2e}".format(f) for f in flops]))
        max_idx, max_flop = calculate_peak_flop(flops)
        max_flop_percentage.append(max_flop)

        bytes_acc = bytes_accessed(ps)
        bandwidth_used = [bytes_acc[i]/ms[i] for i in range(len(ms))]
        print("Bandwidth used", ",".join(["{:2e}".format(f) for f in bandwidth_used]))
        max_idx, max_bw = calculate_peak_bandwidth(bandwidth_used)
        max_bw_percentage.append(max_bw)
        latency = [ms[i]/bytes_acc[i] for i in range(len(ms))]
        print("Latency", ",".join(["{:2e}".format(f) for f in latency]))
        # bytes accessed / elapsed time = bytes/sec = bandwidth used
        # What % of peak bandwidth are you using?
        # Elapsed time / # accesses = average latency
