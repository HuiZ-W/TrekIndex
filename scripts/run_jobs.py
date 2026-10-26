#!/usr/bin/env python3
"""Run a list of single-threaded experiment jobs, one job per dedicated core.

Each worker owns one CPU and runs its jobs one after another under
`nice -n 19 taskset -c CPU /usr/bin/time -v`, so no two timed processes ever
share a core. Jobs are taken longest-estimate first. Where taskset or GNU time
is not available (e.g. macOS), jobs run unpinned and without resource report;
use one worker there (--cores 1).

Jobs file: JSON lines {"name", "cmd", "est_min"}; cmd is run by bash in --cwd.
For every job, <out>/<name>.stdout / .stderr hold its output and
<out>/<name>.status its core, wall time and exit code. A job whose status says
exit 0 is skipped on restart, so the runner can be relaunched after a crash.
"""
import argparse, json, os, pathlib, queue, shutil, subprocess, sys, threading, time

def pin(cpu):
    """taskset -c CPU /usr/bin/time -v where available (Linux), nothing otherwise."""
    out = ["taskset", "-c", cpu] if shutil.which("taskset") else []
    if sys.platform.startswith("linux") and os.path.exists("/usr/bin/time"):
        out += ["/usr/bin/time", "-v"]
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", required=True, type=pathlib.Path)
    ap.add_argument("--cores", required=True, help="comma-separated CPU ids, one worker each")
    ap.add_argument("--cwd", required=True)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    jobs = [json.loads(l) for l in a.jobs.read_text().splitlines() if l.strip()]
    names = [j["name"] for j in jobs]
    assert len(names) == len(set(names)), "duplicate job names"
    todo = queue.Queue()
    for j in sorted(jobs, key=lambda j: -j.get("est_min", 0)):
        st = a.out / f"{j['name']}.status"
        if st.exists() and json.loads(st.read_text()).get("exit") == 0:
            continue
        todo.put(j)
    log = open(a.out / "runner.log", "a")
    lock = threading.Lock()
    def say(msg):
        with lock:
            log.write(time.strftime("%F %T ") + msg + "\n"); log.flush()

    def worker(cpu):
        while True:
            try: j = todo.get_nowait()
            except queue.Empty: return
            say(f"start {j['name']} on cpu {cpu}")
            t = time.time()
            with open(a.out / f"{j['name']}.stdout", "w") as o, open(a.out / f"{j['name']}.stderr", "w") as e:
                rc = subprocess.run(["nice", "-n", "19", *pin(cpu), "bash", "-c", j["cmd"]], cwd=a.cwd, stdout=o, stderr=e,
                                    stdin=subprocess.DEVNULL).returncode
            status = {"name": j["name"], "cpu": int(cpu), "wall_s": round(time.time() - t, 1),
                      "exit": rc, "cmd": j["cmd"], "end": time.strftime("%FT%T")}
            (a.out / f"{j['name']}.status").write_text(json.dumps(status) + "\n")
            say(f"done  {j['name']} exit {rc} in {status['wall_s']:.0f}s")

    say(f"runner start: {todo.qsize()} jobs on cpus {a.cores}")
    threads = [threading.Thread(target=worker, args=(c,)) for c in a.cores.split(",")]
    for th in threads: th.start()
    for th in threads: th.join()
    say("runner finished")

if __name__ == "__main__":
    main()
