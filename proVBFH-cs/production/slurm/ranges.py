#!/usr/bin/env python3
"""Compress a comma list of integers (stdin) into Slurm array ranges a-b,c,..."""
import sys
v = sorted({int(x) for x in sys.stdin.read().replace(',', ' ').split()})
out, i = [], 0
while i < len(v):
    j = i
    while j + 1 < len(v) and v[j + 1] == v[j] + 1:
        j += 1
    out.append(str(v[i]) if i == j else f"{v[i]}-{v[j]}")
    i = j + 1
print(','.join(out))
