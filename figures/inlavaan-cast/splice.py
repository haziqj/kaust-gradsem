# Turn a raw `asciinema rec` of blavaan-*-rec.R into a presentation cast,
# matching the style of inlavaan-kidney.cast:
#   - library() and the model fit are shown in full, then held before "Enter"
#     (1s and 1.5s);
#   - later commands (e.g. print()) are typed out after an empty prompt line;
#   - all output keeps its real timings.
#
# Usage (from the repo root): python3 splice.py raw.cast out.cast COLS ROWS
import json
import sys

raw_path, out_path, cols, rows = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
raw = [json.loads(l) for l in open(raw_path, encoding="utf-8").read().splitlines()[1:]]
txt = [e[2] if e[1] == "o" else "" for e in raw]

i_lib = txt.index("> library(blavaan)\r\n")
i_cmds = [i for i, t in enumerate(txt) if i > i_lib and t.startswith("> ") and t.strip() != ">"]
i_end = next(i for i, t in enumerate(txt) if i > i_cmds[-1] and t.endswith("> "))

# Keystroke gaps borrowed from `timing(fit)` in inlavaan-kidney.cast
KEYS = [0.141, 0.285, 0.133, 0.143, 0.282, 0.203, 0.337, 0.091, 0.167, 0.230]

# Theme etc. from the INLAvaan cast so the GIFs match
hdr = json.loads(open("figures/inlavaan-cast/inlavaan-kidney.cast", encoding="utf-8").readline())
hdr["term"]["cols"], hdr["term"]["rows"] = cols, rows

ev = [
    [0.000, "o", "> library(blavaan)"],
    [1.000, "o", "\r\n"],
]
prev = i_lib
for k, i in enumerate(i_cmds):
    cmd = txt[i][2:].removesuffix("\r\n")
    ev += raw[prev + 1 : i]                         # output of previous command
    ev += [[raw[i][0], "o", "> \r\n"],              # empty prompt line (keeps real delay)
           [0.001, "o", "> "]]
    if k == 0:                                      # model fit: shown in full
        ev += [[1.000, "o", cmd], [1.500, "o", "\r\n"]]
    else:                                           # follow-ups: typed out
        ev += [[1.000 if j == 0 else KEYS[(j - 1) % len(KEYS)], "o", ch] for j, ch in enumerate(cmd)]
        ev += [[0.500, "o", "\r\n"]]
    prev = i
ev += raw[prev + 1 : i_end + 1]                     # last output and final prompt
ev += [[0.005, "x", "0"]]


def fmt(e):
    return f"[{e[0]:.3f}, {json.dumps(e[1])}, {json.dumps(e[2], ensure_ascii=False)}]"


with open(out_path, "w", encoding="utf-8") as f:
    f.write(json.dumps(hdr, separators=(",", ":")) + "\n")
    f.write("\n".join(fmt(e) for e in ev) + "\n")
print(f"{out_path}: {len(ev)} events, {sum(e[0] for e in ev):.2f} s")
