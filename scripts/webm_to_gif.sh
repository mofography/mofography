#!/usr/bin/env bash
# Convert a screen-recorded .webm (or any ffmpeg-readable video) to a README-sized GIF.
# Pipeline: ffmpeg (frames) -> gifski (high-quality palette) -> gifsicle (lossless -O3, optional --lossy).
# Requires: ffmpeg, gifski, gifsicle   (brew install ffmpeg gifski gifsicle)
set -euo pipefail

usage() {
  echo "usage: $0 in.webm [out.gif] [--width 900] [--fps 12] [--quality 90] [--lossy 0] [--start S] [--duration S] [--crop W:H:X:Y]"
  echo "  --lossy N   gifsicle lossy level (0 = lossless only; 30-80 = smaller, mild artefacts)"
  exit 1
}
[ $# -ge 1 ] || usage
in="$1"; shift
out="${in%.*}.gif"
[ $# -gt 0 ] && [[ "$1" != --* ]] && { out="$1"; shift; }
width=900; fps=12; quality=90; lossy=0; start=""; dur=""; crop=""
while [ $# -gt 0 ]; do
  case "$1" in
    --width) width="$2";; --fps) fps="$2";; --quality) quality="$2";;
    --lossy) lossy="$2";; --start) start="$2";; --duration) dur="$2";; --crop) crop="$2";;
    *) usage;;
  esac; shift 2
done

tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
trim=(); [ -n "$start" ] && trim+=(-ss "$start"); [ -n "$dur" ] && trim+=(-t "$dur")
# libvpx-vp9 decoder keeps the alpha channel of transparent webm; harmless otherwise
ffmpeg -v error -y -c:v libvpx-vp9 ${trim[@]+"${trim[@]}"} -i "$in" \
  -vf "${crop:+crop=$crop,}fps=$fps,scale=$width:-2:flags=lanczos" "$tmp/f%05d.png" 2>/dev/null \
  || ffmpeg -v error -y ${trim[@]+"${trim[@]}"} -i "$in" -vf "${crop:+crop=$crop,}fps=$fps,scale=$width:-2:flags=lanczos" "$tmp/f%05d.png"

gifski --quiet --fps "$fps" --quality "$quality" -o "$tmp/raw.gif" "$tmp"/f*.png
opt=(-O3); [ "$lossy" -gt 0 ] && opt+=(--lossy="$lossy")
gifsicle "${opt[@]}" "$tmp/raw.gif" -o "$out"

sz() { du -h "$1" | cut -f1; }
echo "raw $(sz "$tmp/raw.gif") -> $out $(sz "$out")  (${width}px, ${fps}fps, lossy=$lossy)"
