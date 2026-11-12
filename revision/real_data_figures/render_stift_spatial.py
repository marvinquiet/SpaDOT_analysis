from collections import defaultdict
from csv import DictReader
from hashlib import sha256
import json
from math import ceil, isfinite
from pathlib import Path
import os
import sys

from reportlab.lib.colors import HexColor
from reportlab.lib.units import cm
from reportlab.pdfbase.pdfmetrics import registerFont, stringWidth
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas


DATASETS = {
    "Axolotl": ("axolotl_figure5_stift_20261109", ("St.44", "St.54", "St.57"),
                (1477, 2929, 4410), (6, 7, 8)),
    "ChickenHeart": ("chicken_figure4_stift_20261109", ("Day 4", "Day 7", "Day 10", "Day 14"),
                     (747, 1966, 1916, 1967), (5, 7, 7, 6)),
    "MouseEmbryo": ("mouse_figure3b_stift_20261109", ("E12.5", "E13.5", "E15.5"),
                    (2397, 2312, 3546), (5, 5, 5)),
}
# Default ggplot2 discrete hue palettes for the displayed label counts.
PALETTES = {
    5: ("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"),
    7: ("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7"),
    8: ("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"),
}
WIDTH, HEIGHT = 8.59 * cm, 3.43 * cm


def read_stages(root, dataset, chicken_day4):
    folder, titles, spot_counts, cluster_counts = DATASETS[dataset]
    directory = root / folder
    path = directory / "matched_stift_reference.csv"
    audit_path = directory / "prediction_audit.json"
    if dataset == "MouseEmbryo":
        audit_path = directory / "computed" / "prediction_audit.json"
    audit = json.loads(audit_path.read_text())
    assert sha256(path.read_bytes()).hexdigest() == audit["matched_stift_table_sha256"]

    historical = {}
    if dataset == "ChickenHeart":
        with chicken_day4.open() as handle:
            history_rows = list(DictReader(handle))
        historical = {row["spot"]: (float(row["spatial_x"]), float(row["spatial_y"])) for row in history_rows}
        assert len(historical) == len(history_rows) == spot_counts[0]

    seen, replaced = set(), set()
    stages = defaultdict(list)
    with path.open() as handle:
        for row in DictReader(handle):
            spot, stage = row["spot"], row["timepoint"]
            label_stage, raw_label = row["STIFT_KMeans_per_timepoint"].split("_")
            assert spot not in seen and stage == label_stage
            seen.add(spot)
            x, y = float(row["spatial_x"]), float(row["spatial_y"])
            if dataset == "ChickenHeart":
                if stage == "t0":
                    x, y = historical[spot]
                    replaced.add(spot)
                x, y = y, x
            elif dataset == "MouseEmbryo":
                x, y = (y if stage == "t1" else -y), -x
            assert isfinite(x) and isfinite(y), spot
            stages[stage].append((x, y, int(raw_label)))
    assert replaced == set(historical)
    assert set(stages) == {f"t{i}" for i in range(len(titles))}
    for i, (nspots, k) in enumerate(zip(spot_counts, cluster_counts)):
        stage = stages[f"t{i}"]
        assert len(stage) == nspots and {row[2] for row in stage} == set(range(k))
    return [stages[f"t{i}"] for i in range(len(titles))]


def add_text(pdf, text, x, y, size, *, bold=False, align="center"):
    pdf.setFillColor(HexColor("#171717"))
    pdf.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    if align == "center":
        pdf.drawCentredString(x, y, text)
    else:
        pdf.drawString(x, y, text)


def draw_stage(pdf, rows, center, area_width, colors, dataset):
    xmin, xmax = min(r[0] for r in rows), max(r[0] for r in rows)
    ymin, ymax = min(r[1] for r in rows), max(r[1] for r in rows)
    xmid, ymid = (xmin + xmax) / 2, (ymin + ymax) / 2
    area_height, center_y = (63.5, 52) if dataset == "MouseEmbryo" else (50, 50)
    sx, sy = area_width / (xmax - xmin), area_height / (ymax - ymin)
    if dataset == "ChickenHeart":
        sx = sy = min(sx, sy)
    radius = {"Axolotl": 0.36, "ChickenHeart": 0.34, "MouseEmbryo": 0.40}[dataset]
    for x, y, label in rows:
        pdf.setFillColor(colors[label])
        pdf.circle(center + (x - xmid) * sx, center_y + (y - ymid) * sy,
                   radius, stroke=0, fill=1)


def draw_legend(pdf, colors):
    columns = ceil(len(colors) / 2)
    gap = 13.5
    left = (WIDTH - (columns - 1) * gap) / 2 - 3
    for label, color in enumerate(colors):
        x = left + (label // 2) * gap
        y = 15 - (label % 2) * 7.6
        pdf.setFillColor(color)
        pdf.circle(x, y, 1.32, stroke=0, fill=1)
        add_text(pdf, str(label), x + 4.2, y - 2.1, 6, align="left")


def draw_mouse_legend(pdf, colors):
    title = "STIFT Predicted Labels"
    title_size, label_size = 6.1, 5.4
    font_path = os.environ.get("SPADOT_BOLD_FONT")
    font = "Helvetica-Bold"
    if font_path:
        font = "STIFTLegendBold"
        registerFont(TTFont(font, font_path))
    title_width = stringWidth(title, font, title_size)
    first_key, key_spacing, text_offset = title_width + 5.3, 13.55, 5.5
    group_width = first_key + (len(colors) - 1) * key_spacing + text_offset + stringWidth("4", "Helvetica", label_size)
    left = (WIDTH - group_width) / 2
    pdf.setFillColor(HexColor("#171717"))
    pdf.setFont(font, title_size)
    pdf.drawString(left, 5.4, title)
    for label, color in enumerate(colors):
        x = left + first_key + label * key_spacing
        pdf.setFillColor(color)
        pdf.circle(x, 7.6, 1.35, stroke=0, fill=1)
        add_text(pdf, str(label), x + text_offset, 5.4, label_size, align="left")


def render(root, out, dataset, chicken_day4):
    _, titles, counts, cluster_counts = DATASETS[dataset]
    stages = read_stages(root, dataset, chicken_day4)
    colors = [HexColor(c) for c in PALETTES[max(cluster_counts)]]
    path = out / f"{dataset}_STIFT_prediction.pdf"
    pdf = Canvas(str(path), pagesize=(WIDTH, HEIGHT), invariant=True)
    pdf.setTitle(f"{dataset} STIFT Prediction")
    is_mouse = dataset == "MouseEmbryo"
    if is_mouse:
        pdf.setFillAlpha(1.0)
        pdf.setStrokeAlpha(1.0)
    if not is_mouse:
        add_text(pdf, "STIFT Prediction", WIDTH / 2, 88.0, 7.0)
    outer_margin = 1.2 if is_mouse else 8
    cell_width = (WIDTH - 2 * outer_margin) / len(stages)
    heading_y, heading_size, gap = (88.8, 4.75, 8.8) if is_mouse else (80.3, 5.7, 8)
    for i, (title, rows) in enumerate(zip(titles, stages)):
        center = outer_margin + (i + 0.5) * cell_width
        add_text(pdf, title, center, heading_y, heading_size, bold=True)
        draw_stage(pdf, rows, center, cell_width - gap, colors, dataset)
    if is_mouse:
        draw_mouse_legend(pdf, colors)
    else:
        draw_legend(pdf, colors)
    pdf.showPage()
    pdf.save()
    print(dataset, sum(counts), "spots;", cluster_counts, "clusters per stage")


if __name__ == "__main__":
    root, output, day4 = map(Path, sys.argv[1:])
    for name in DATASETS:
        render(root, output, name, day4)
