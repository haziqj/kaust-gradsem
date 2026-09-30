# Turn a raw `asciinema rec` of blavaan-*-rec.R into a presentation cast:
# each command is shown in full, then held before "Enter" (1s for library(),
# 1.5s for the fit, as in inlavaan-kidney.cast); all output keeps real timings.
#
# Usage: python3 splice.py raw.cast out.cast COLS ROWS
import json
import sys

raw_path, out_path, cols, rows = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
raw = [json.loads(l) for l in open(raw_path, encoding="utf-8").read().splitlines()[1:]]
txt = [e[2] if e[1] == "o" else "" for e in raw]

i_lib = txt.index("> library(blavaan)\r\n")
i_cmd = next(i for i, t in enumerate(txt) if t.startswith("> fit_blav"))
i_end = next(i for i, t in enumerate(txt) if t.endswith("> ") and i > i_cmd)

# Theme etc. from the INLAvaan cast so both GIFs match
hdr = json.loads(open("figures/inlavaan-cast/inlavaan-kidney.cast", encoding="utf-8").readline())
hdr["term"]["cols"], hdr["term"]["rows"] = cols, rows

ev = [
    [0.000, "o", "> library(blavaan)"],
    [1.000, "o", "\r\n"],
    *raw[i_lib + 1 : i_cmd],                      # startup messages
    [0.003, "o", "> \r\n"],                        # empty prompt line
    [0.001, "o", "> "],
    [1.000, "o", txt[i_cmd][2:].removesuffix("\r\n")],
    [1.500, "o", "\r\n"],
    *raw[i_cmd + 1 : i_end + 1],                  # sampler output, final prompt
    [0.005, "x", "0"],
]


def fmt(e):
    return f"[{e[0]:.3f}, {json.dumps(e[1])}, {json.dumps(e[2], ensure_ascii=False)}]"


with open(out_path, "w", encoding="utf-8") as f:
    f.write(json.dumps(hdr, separators=(",", ":")) + "\n")
    f.write("\n".join(fmt(e) for e in ev) + "\n")
print(f"{out_path}: {len(ev)} events, {sum(e[0] for e in ev):.2f} s")
