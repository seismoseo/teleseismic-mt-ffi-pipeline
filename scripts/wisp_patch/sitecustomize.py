# Runtime patch (no edit of vendored WISP): widen the strong-motion station distance
# cutoff so regional stations beyond WISP's default 2 deg are retained. Eagerly imports
# and patches ffm.data_processing; the later `import ffm.data_processing` in the CLI gets
# the already-patched module from sys.modules.
import os
CUT = float(os.environ.get("WISP_STRONG_MAXDIST_DEG", "4.5"))
try:
    import inspect
    import ffm.data_processing as dp
    d = dp.__dict__
    src = inspect.getsource(d["__select_str_files"])
    src = src.replace("distance = 2 if time_shift < 50 else 4", f"distance = {CUT}")
    src = src.replace("def __select_str_files", "def _patched_select_str_files")
    ns = dict(d)
    exec(src, ns)
    d["__select_str_files"] = ns["_patched_select_str_files"]
    print(f"WISP PATCH active: strong-motion distance cutoff -> {CUT} deg", flush=True)
except Exception as e:
    print("WISP PATCH FAILED:", e, flush=True)
