#!/usr/bin/env python3
"""Resubmit the lines of a jobs.list that started but did not finish (no
"done" marker; a time.log present or fed before by feed_missing.sh) and are not pending or running in any
array of the given job name (tasks stuck COMPLETING on dead nodes count as
failed: their run never finished, so it never wrote "done"). Unlike resubmit.sh this works while the arrays
are still active (nodes failing during the production).

  resubmit_failed.py <name> <binary> <jobs.list> <time limit> <cap> [nice]

Array task i of an array runs line i + OFFSET of jobs.list. OFFSET is 0 for
every array except the ones listed in OFFSETS below (submitted by hand with
OFFSET=10000).
"""
import os
import subprocess
import sys

OFFSETS = {"48603439": 10000, "48603568": 10000}
HERE = os.path.dirname(os.path.abspath(__file__))


def sh(cmd):
    return subprocess.run(cmd, check=True, text=True, capture_output=True).stdout


def main():
    name, binary, jobs_list, tlim, cap = sys.argv[1:6]
    nice = sys.argv[6] if len(sys.argv) > 6 else "0"
    cap = int(cap)
    user = os.environ["USER"]
    active = set()
    for line in sh(["squeue", "-u", user, "-h", "-r", "-n", name,
                    "-t", "PD,R,CF", "-o", "%F %K"]).split("\n"):
        f = line.split()
        if len(f) == 2 and f[1].isdigit():
            active.add(int(f[1]) + OFFSETS.get(f[0], 0))
    dirs = [l.strip() for l in open(jobs_list)]
    # lines handed to Slurm before by feed_missing.sh (some die before
    # their time.log is written, e.g. NODE_FAIL at launch)
    fed = set()
    if os.path.exists(jobs_list + ".fed"):
        fed = {int(l) for l in open(jobs_list + ".fed") if l.strip().isdigit()}
    failed = [i + 1 for i, d in enumerate(dirs)
              if i + 1 not in active
              and not os.path.exists(os.path.join(d, "done"))
              and (os.path.exists(os.path.join(d, "time.log")) or i + 1 in fed)]
    queued = len(sh(["squeue", "-u", user, "-h", "-r"]).splitlines())
    room = cap - queued
    if not failed:
        print(f"{name}: no failed lines")
        return
    if room <= 0:
        print(f"{name}: {len(failed)} failed lines, queue at cap ({cap})")
        return
    failed = failed[:room]
    p = subprocess.run([os.path.join(HERE, "ranges.py")], input=",".join(map(str, failed)),
                       text=True, capture_output=True, check=True)
    ranges = p.stdout.strip()
    root = os.path.dirname(os.path.dirname(os.path.dirname(jobs_list)))
    logs = os.path.join(root, "..", "logs", os.path.basename(root))
    bad = ""
    try:
        bad = ",".join(l.strip() for l in open("/ptmp/mpp/akarlber/cs-production/bad_nodes") if l.strip())
    except OSError:
        pass
    cmd = ["sbatch", "--parsable", f"--nice={nice}", f"--array={ranges}", f"--time={tlim}",
           "-J", name, f"--export=ALL,BIN={binary},LIST={jobs_list},OFFSET=0",
           "-o", os.path.join(logs, "%x.%A_%a.out")]
    if bad:
        cmd.append(f"--exclude={bad}")
    cmd.append(os.path.join(HERE, "run_array.sh"))
    jid = sh(cmd).strip()
    msg = f"{name}: resubmitted {len(failed)} failed lines as array {jid}"
    print(msg)
    with open(jobs_list + ".submitted", "a") as f:
        f.write(msg + "\n")


if __name__ == "__main__":
    main()
