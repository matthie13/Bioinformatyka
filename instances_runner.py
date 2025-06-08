#!/usr/bin/env python3
import subprocess
import statistics
import argparse
import time

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


def run_experiments(inst_file, orig_file, ants_list, sbh_script, common_args):
    instances, originals = parse_instances(inst_file, orig_file)
    total_tasks = len(instances) * len(ants_list)
    completed_tasks = 0
    start_time = time.time()
    estimated_total_printed = False
    results = {}

    for ants in ants_list:
        dists = []
        for idx, (inst, orig) in enumerate(zip(instances, originals), start=1):
            with open('tmp_instance.txt', 'w') as f:
                f.write("\n".join(inst) + "\n")
            with open('tmp_original.txt', 'w') as f:
                f.write(orig + "\n")

            cmd = ['python', sbh_script, 'tmp_instance.txt', '-o', 'tmp_original.txt',
                   '-a', str(ants)] + common_args
            proc = subprocess.run(cmd, capture_output=True, text=True)
            # Parse Levenshtein
            for line in reversed(proc.stdout.splitlines()):
                if 'Levenshtein distance to original:' in line:
                    d = int(line.split(':')[-1].strip())
                    dists.append(d)
                    break

            # Progress update
            completed_tasks += 1
            elapsed = time.time() - start_time
            avg_time = elapsed / completed_tasks
            remaining = total_tasks - completed_tasks
            eta = remaining * avg_time

            # Print estimated total once after first task
            if not estimated_total_printed and completed_tasks > 0:
                estimated_total = avg_time * total_tasks
                print(f"Estimated total time: {format_time(estimated_total)} for {total_tasks} tasks")
                estimated_total_printed = True

            print(f"Progress: {completed_tasks}/{total_tasks} | Elapsed: {format_time(elapsed)} | ETA: {format_time(eta)}")

        results[ants] = statistics.mean(dists)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--instances', default='instance-50-500.txt')
    parser.add_argument('--originals', default='original-50-500.txt')
    parser.add_argument('--sbh', default='SBH.py')
    parser.add_argument('--ants', nargs='+', type=int, default=[50,100,150,200,250])
    parser.add_argument('--alpha', type=float, default=0.7)
    parser.add_argument('--beta', type=float, default=6.0)
    parser.add_argument('--rho', type=float, default=0.8)
    parser.add_argument('-q', '--Q', type=float, default=20.0)
    parser.add_argument('--tau0', type=float, default=0.7)
    parser.add_argument('-t', '--time', type=int, default=5)
    parser.add_argument('-i', '--iter', type=int, default=0)
    args = parser.parse_args()

    common = [
        '--alpha', str(args.alpha),
        '--beta', str(args.beta),
        '--rho', str(args.rho),
        '-q', str(args.Q),
        '--tau0', str(args.tau0),
        '-t', str(args.time),
        '-i', str(args.iter)
    ]
    print("Starting batch run...")
    results = run_experiments(
        args.instances, args.originals,
        args.ants, args.sbh, common
    )
    print("\n=== Average Levenshtein by number of ants ===")
    for ants, avg in results.items():
        print(f"Ants={ants}: avg Levenshtein={avg:.2f}")

if __name__ == '__main__':
    main()
