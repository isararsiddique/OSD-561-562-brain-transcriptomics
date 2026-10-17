"""Inspect scientific PDF embedding, page boundaries and integrated headings."""
import json
import hashlib
from pypdf import PdfReader
import pdfplumber
from pipeline.paths import stage_root, results_root

records = []
for group in ["core", "extended", "deep"]:
    for path in sorted((stage_root(group) / "figures").glob("*.pdf")):
        reader = PdfReader(path)
        assert len(reader.pages) == 1, path
        fonts = []
        for ref in reader.pages[0]["/Resources"]["/Font"].values():
            font = ref.get_object()
            base = font["/DescendantFonts"][0].get_object() if "/DescendantFonts" in font else font
            descriptor = base.get("/FontDescriptor")
            assert descriptor is not None, (path, font)
            descriptor = descriptor.get_object()
            key = next(k for k in ["/FontFile", "/FontFile2", "/FontFile3"] if k in descriptor)
            raw = descriptor[key].get_object().get_data()
            assert raw, path
            fonts.append({"name": str(font["/BaseFont"]), "embedded_bytes": len(raw)})
        with pdfplumber.open(path) as doc:
            page = doc.pages[0]
            chars = [c for c in page.chars if c["text"].strip()]
            assert chars, path
            if group != "core":
                top = min(c["top"] for c in chars)
                heading = [c for c in chars if c["top"] < top + 1]
                def black(c):
                    color = c["non_stroking_color"]
                    return abs(color) < 1e-8 if isinstance(color, (int, float)) else all(abs(float(v)) < 1e-8 for v in color)
                assert heading and all(black(c) for c in heading), path
            assert all(c["x0"] >= -1 and c["x1"] <= page.width + 1 and c["top"] >= -1 and c["bottom"] <= page.height + 1 for c in chars), path
        records.append({"path": path.relative_to(results_root()).as_posix(), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "pages": 1, "fonts": fonts, "text_within_page": True, "black_integrated_heading": group != "core"})
assert len([r for r in records if r["path"].startswith("extended/")]) == 8
assert len([r for r in records if r["path"].startswith("deep/")]) == 3
(results_root() / "metadata/figure_file_checks.json").write_text(json.dumps(records, indent=2) + "\n")
print(f"Inspected {len(records)} PDF figures: embedded fonts, single pages and page bounds.")
