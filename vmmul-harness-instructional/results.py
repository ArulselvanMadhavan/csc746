import pandas as pd
import csv
import matplotlib.pyplot as plt
import math

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

if __name__ == "__main__":
    filenames = ["basic.txt","openmp.txt"]
    for idx in range(len(filenames)):
        fname = filenames[idx]
        lines = read_lines(fname)
        lines = lines[2:]
        schedules = infer_ps(lines)
        print(schedules)
