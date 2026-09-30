QMD  := index.qmd
OUT  := docs
TMP  := _pdf
SIZE := 1250x760

.PHONY: all pdf html clean

all: pdf html

# Static (no-GIF) render into a scratch dir, then print it with decktape
pdf:
	quarto render $(QMD) -M nogif:true --output-dir $(TMP)
	decktape reveal --size $(SIZE) $(TMP)/index.html $(OUT)/slides.pdf
	rm -rf $(TMP)

# Normal render into docs/. Rendering a single file (not the whole project)
# leaves everything else in docs/ in place.
html:
	quarto render $(QMD) -M nogif:false

clean:
	rm -rf $(TMP)
