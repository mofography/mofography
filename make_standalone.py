"""Inline all external CDN scripts into the generated HTML for offline use."""
import os
from pathlib import Path

HTML_PATH = Path(os.environ.get("HTML_PATH", "/Users/arungopalan/repos/mofography/P1-IRMOF-1_padded_threejs.html"))
OUT_PATH = Path(os.environ.get("OUT_PATH", "/Users/arungopalan/repos/mofography/P1-IRMOF-1_padded_threejs_standalone.html"))
DEPS_DIR = Path("/tmp/threejs_deps")

# Map CDN URLs -> local downloaded file names
CDN_MAP = {
    "https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js": "three.min.js",
    "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js": "OrbitControls.js",
    "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/shaders/CopyShader.js": "CopyShader.js",
    "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/shaders/SSAOShader.js": "SSAOShader.js",
    "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/math/SimplexNoise.js": "SimplexNoise.js",
    "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/postprocessing/EffectComposer.js": "EffectComposer.js",
    "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/postprocessing/ShaderPass.js": "ShaderPass.js",
    "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/postprocessing/MaskPass.js": "MaskPass.js",
    "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/postprocessing/RenderPass.js": "RenderPass.js",
    "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/postprocessing/SSAOPass.js": "SSAOPass.js",
    "https://unpkg.com/gif.js@0.2.0/dist/gif.js": "gif.js",
}

def main():
    html = HTML_PATH.read_text(encoding="utf-8")

    # 1. Inline each <script src="CDN"> tag
    for cdn_url, filename in CDN_MAP.items():
        local_path = DEPS_DIR / filename
        js_content = local_path.read_text(encoding="utf-8")
        # Escape any </script> inside the JS to prevent breaking the HTML
        js_content = js_content.replace("</script>", "<\\/script>")
        js_content = js_content.replace("</SCRIPT>", "<\\/SCRIPT>")

        old_tag = f'<script src="{cdn_url}"></script>'
        new_tag = f"<script>\n{js_content}\n</script>"

        if old_tag not in html:
            print(f"Warning: {cdn_url} not found in HTML")
            continue

        html = html.replace(old_tag, new_tag, 1)
        print(f"Inlined {filename} ({len(js_content):,} bytes)")

    # 2. Handle gif.worker.js dynamic fetch -> replace the whole function body
    worker_js = (DEPS_DIR / "gif.worker.js").read_text(encoding="utf-8")
    worker_js_escaped = worker_js.replace("\\", "\\\\").replace("`", "\\`").replace("$", "\\$")

    old_func = '''    async function getGifWorkerUrl() {
      if (gifWorkerBlobUrl) return gifWorkerBlobUrl;
      try {
        const resp = await fetch("https://unpkg.com/gif.js@0.2.0/dist/gif.worker.js");
        if (!resp.ok) throw new Error("fetch failed");
        const text = await resp.text();
        const blob = new Blob([text], { type: "application/javascript" });
        gifWorkerBlobUrl = URL.createObjectURL(blob);
        return gifWorkerBlobUrl;
      } catch(e) {
        return "https://unpkg.com/gif.js@0.2.0/dist/gif.worker.js";
      }
    }'''

    new_func = f'''    async function getGifWorkerUrl() {{
      if (gifWorkerBlobUrl) return gifWorkerBlobUrl;
      const text = `{worker_js_escaped}`;
      const blob = new Blob([text], {{ type: "application/javascript" }});
      gifWorkerBlobUrl = URL.createObjectURL(blob);
      return gifWorkerBlobUrl;
    }}'''

    if old_func in html:
        html = html.replace(old_func, new_func, 1)
        print("Inlined gif.worker.js via Blob (inline function replaced)")
    else:
        print("Warning: getGifWorkerUrl() pattern not found; trying string replacement")
        html = html.replace(
            '"https://unpkg.com/gif.js@0.2.0/dist/gif.worker.js"',
            '"data:text/javascript;base64,"',  # fallback won't work well
        )

    # 3. Write output
    OUT_PATH.write_text(html, encoding="utf-8")
    original_size = HTML_PATH.stat().st_size
    new_size = OUT_PATH.stat().st_size
    print(f"\nDone.")
    print(f"  Original: {original_size:,} bytes")
    print(f"  Standalone: {new_size:,} bytes")
    print(f"  Saved to: {OUT_PATH}")


if __name__ == "__main__":
    main()
