#!/usr/bin/env python3
import subprocess
import statistics
import argparse
import time
import os
import tempfile

def format_time(seconds):
    """Format seconds as H:MM:SS or M:SS"""
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    if h > 0:
        return f"{h}h {m:02d}m {s:02d}s"
    else:
        return f"{m}m {s:02d}s"

def parse_instances(inst_path, orig_path):
    instances = []
    originals = []
    with open(inst_path) as f:
        inst_lines = [l.rstrip() for l in f]
    with open(orig_path) as f:
        orig_lines = [l.rstrip() for l in f]

    def split_chunks(lines):
        chunks = []
        current = []
        for l in lines:
            if l.startswith('# Instance'):
                if current:
                    chunks.append(current)
                current = []
                continue
            if l.strip():
                current.append(l)
        if current:
            chunks.append(current)
        return chunks

    inst_chunks = split_chunks(inst_lines)
    orig_chunks = split_chunks(orig_lines)
    if len(inst_chunks) != len(orig_chunks):
        raise ValueError("Mismatch between instance and original counts")

    for inst, orig in zip(inst_chunks, orig_chunks):
        instances.append(inst)
        originals.append(orig[0])
    return instances, originals

def run_experiments(inst_file, orig_file, ants_list, times_list, betas_list, sbh_script, common_args):
    instances, originals = parse_instances(inst_file, orig_file)
    total_tasks = len(instances) * len(ants_list) * len(times_list) * len(betas_list)
    completed_tasks = 0
    start_time = time.time()
    estimated_total_printed = False
    results = {}

    for ants in ants_list:
        for t_budget in times_list:
            for beta in betas_list:
                dists = []
                for inst, orig in zip(instances, originals):
                    # Write to unique temporary files
                    with tempfile.NamedTemporaryFile('w+', delete=False) as tf_inst, \
                         tempfile.NamedTemporaryFile('w+', delete=False) as tf_orig:
                        tf_inst.write("\n".join(inst) + "\n")
                        tf_orig.write(orig + "\n")
                        tf_inst.flush()
                        tf_orig.flush()
                        tmp_inst = tf_inst.name
                        tmp_orig = tf_orig.name

                    cmd = [
                        'python', sbh_script,
                        tmp_inst, '-o', tmp_orig,
                        '-a', str(ants), '-t', str(t_budget),
                        '--beta', str(beta)
                    ] + common_args
                    proc = subprocess.run(cmd, capture_output=True, text=True)

                    # Clean up temp files
                    os.remove(tmp_inst)
                    os.remove(tmp_orig)

                    # Parse Levenshtein from output
                    last_d = None
                    for line in reversed(proc.stdout.splitlines()):
                        if 'Levenshtein distance to original:' in line:
                            last_d = int(line.split(':')[-1].strip())
                            dists.append(last_d)
                            break

                    # Progress reporting
                    completed_tasks += 1
                    elapsed = time.time() - start_time
                    avg_time = elapsed / completed_tasks
                    remaining = total_tasks - completed_tasks
                    eta = remaining * avg_time
                    if not estimated_total_printed:
                        estimated_total = avg_time * total_tasks
                        print(f"Estimated total time: {format_time(estimated_total)} "
                              f"for {total_tasks} tasks")
                        estimated_total_printed = True

                    print(f"Progress [ants={ants}, time={t_budget}s, beta={beta}]: "
                          f"{completed_tasks}/{total_tasks} | Elapsed: {format_time(elapsed)} "
                          f"| ETA: {format_time(eta)} | d={last_d}")

                results[(ants, t_budget, beta)] = statistics.mean(dists)

    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--instances', default='instance-50-700.txt')
    parser.add_argument('--originals', default='original-50-700.txt')
    parser.add_argument('--sbh', default='SBH.py')
    parser.add_argument('--ants', nargs='+', type=int, default=[100])
    parser.add_argument('--time', '-t', nargs='+', type=int, default=[10])
    parser.add_argument('--betas', nargs='+', type=float, default=[8.0])
    parser.add_argument('--alpha', type=float, default=0.7)
    parser.add_argument('--rho', type=float, default=0.8)
    parser.add_argument('-q','--Q', type=float, default=20.0)
    parser.add_argument('-i','--iter', type=int, default=0)
    args = parser.parse_args()

    common = [
        '--alpha', str(args.alpha),
        '--rho', str(args.rho),
        '-q', str(args.Q),
        '-i', str(args.iter)
    ]

    print("Starting batch run...")
    results = run_experiments(
        args.instances, args.originals,
        args.ants, args.time, args.betas,
        args.sbh, common
    )

    print("\n=== Average Levenshtein by (ants, time, beta) ===")
    for (ants, t_budget, beta), avg in sorted(results.items()):
        print(f"Ants={ants}, Time={t_budget}s, Beta={beta}: avg Levenshtein = {avg:.2f}")

if __name__ == '__main__':
    main()
