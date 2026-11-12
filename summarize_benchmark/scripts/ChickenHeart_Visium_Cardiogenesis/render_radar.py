from csv import DictReader
from math import cos, isclose, isfinite, pi, sin
from pathlib import Path
from statistics import mean
import sys

from reportlab.lib.colors import HexColor
from reportlab.lib.units import cm
from reportlab.pdfgen.canvas import Canvas

inputs = Path(sys.argv[1])
metrics = ("ARI", "AMI", "CHAOS", "median-LISI", "PAS", "H")
methods = ("SpaDOT", "spatiAlign", "spatialPCA", "spaVAE", "SPIRAL", "SEDR",
           "STADIA", "PRECAST", "spaGCN", "STAligner", "GraphST", "STAMP", "STIFT")
palette = ("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0",
           "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae", "#ABDDDE")
with (inputs / "palette_reference.csv").open() as file:
    archive_order = tuple(dict.fromkeys(row["method"] for row in DictReader(file)))
assert set(archive_order) == set(methods[:-1])
colors = {method: HexColor(color) for method, color in zip(archive_order, palette)}
colors["STIFT"] = HexColor("#6A4C93")

with (inputs / "ChickenHeart.csv").open() as file:
    rows = list(DictReader(file))
scores = {(row["method"], row["index"], row["timepoint"]): float(row["value"]) for row in rows}
stages = ("t0", "t1", "t2", "t3")
assert len(scores) == len(rows) == 390 and all(map(isfinite, scores.values()))
for metric in metrics:
    for method in methods:
        value = mean(scores[method, metric, stage] for stage in stages)
        assert isclose(value, scores[method, metric, "avg"], abs_tol=1e-12)

lower_is_better = {"CHAOS", "median-LISI", "PAS"}
labels = ("ARI", "AMI", "1 - CHAOS", "1 - median-LISI", "1 - PAS", "H")
vectors = [(cos(pi / 2 + i * pi / 3), sin(pi / 2 + i * pi / 3)) for i in range(6)]
grid_values = (0, 0.25, 0.5, 0.75, 1)

average_width = 7.04 * cm
day_width = 7.55 * cm
width, height = average_width + day_width, 9 * cm
pdf = Canvas("ChickenHeart_benchmark_radar.pdf", pagesize=(width, height), invariant=True)
pdf.setTitle("Chicken Heart benchmark radar")
pdf.setLineJoin(1)


def normalized_scores(stage):
    # Follow the original radar: scale each metric within the plotted timepoint.
    result = {}
    for metric in metrics:
        values = [scores[method, metric, stage] for method in methods]
        low, high = min(values), max(values)
        assert high > low, (stage, metric)
        fractions = [(value - low) / (high - low) for value in values]
        result[metric] = [1 - value for value in fractions] if metric in lower_is_better else fractions
    return result


def polygon(points):
    path = pdf.beginPath()
    path.moveTo(*points[0])
    for xy in points[1:]:
        path.lineTo(*xy)
    path.close()
    pdf.drawPath(path)


def text(x, y, label, size=9, color="#202020", align="left", bold=False):
    pdf.setFillColor(HexColor(color))
    pdf.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    {"left": pdf.drawString, "center": pdf.drawCentredString, "right": pdf.drawRightString}[align](x, y, label)


def draw_radar(stage, title, cx, cy, radius, title_y, *, large=False):
    normalized = normalized_scores(stage)
    label_size = 8.0 if large else 5.9

    def point(index, distance):
        x, y = vectors[index]
        return cx + distance * x, cy + distance * y

    def radial(score):
        # Retain the original radar's inner zero ring.
        return radius * (0.2 + 0.8 * score)

    pdf.setStrokeColor(HexColor("#c7c7c7"))
    pdf.setLineWidth(0.5 if large else 0.35)
    for score in grid_values:
        polygon([point(i, radial(score)) for i in range(6)])
    for i in range(6):
        pdf.line(*point(i, radial(0)), *point(i, radius))

    for method in (*methods[1:], methods[0]):
        idx = methods.index(method)
        points = [point(i, radial(normalized[metric][idx])) for i, metric in enumerate(metrics)]
        is_spadot = method == "SpaDOT"
        pdf.setStrokeColor(colors[method])
        pdf.setFillColor(colors[method])
        pdf.setStrokeAlpha(0.95 if is_spadot else 0.80)
        pdf.setFillAlpha(1 if is_spadot else 0.85)
        stroke, dot = (1.8, 1.35) if is_spadot else (0.95, 0.88)
        scale = 1 if large else 0.66
        pdf.setLineWidth(stroke * scale)
        polygon(points)
        for xy in points:
            pdf.circle(*xy, dot * scale, stroke=0, fill=1)
    pdf.setStrokeAlpha(1)
    pdf.setFillAlpha(1)

    # Move diagonal axis labels above and below their vertices to fit the
    # two narrow columns without shrinking the day polygons.
    shoulder = radius * 0.76
    upper = cy + radius * 0.5 + (11 if large else 8)
    lower = cy - radius * 0.5 - (23 if large else 17.5)
    label_positions = (
        (cx, cy + radius + 4.5),
        (cx - shoulder, upper),
        (cx - shoulder, lower),
        (cx, cy - radius - (11 if large else 7.5)),
        (cx + shoulder, lower),
        (cx + shoulder, upper),
    )
    for label, (x, y) in zip(labels, label_positions):
        text(x, y, label, size=label_size, align="center")
    for score in grid_values if large else (0, 0.5, 1):
        x, y = point(0, radial(score))
        text(x - 3, y - 1.7, f"{score:g}", size=6.4 if large else 5.0,
             color="#6b6b6b", align="right")
    text(cx, title_y, title,
         size=12 if large else 8.4, align="center", bold=True)


draw_radar("avg", "Average", average_width / 2, 153, 72, 243, large=True)

day_centers = (average_width + day_width / 4, average_width + 3 * day_width / 4)
for index, (stage, day) in enumerate(zip(stages, (4, 7, 10, 14))):
    row, column = divmod(index, 2)
    cy = (186.5, 60.5)[row]
    draw_radar(stage, f"Day {day}", day_centers[column], cy, 48.5, cy + 59.7)

# Center the legend in the space below the Average radar.
legend_center_y = 34
pdf.setStrokeColor(HexColor("#c9c9c9"))
pdf.setLineWidth(0.45)
pdf.line(43.5, legend_center_y - 18, 43.5, legend_center_y + 18)
for idx, method in enumerate(methods):
    if idx == 0:
        x, y = 9, legend_center_y
        dot, label_offset = 2.8, 5.1
    else:
        col, row = divmod(idx - 1, 3)
        x, y = (50, 89.5, 128, 166)[col], legend_center_y + (1 - row) * 11.5
        dot, label_offset = 1.8, 4.0
    pdf.setFillColor(colors[method])
    pdf.circle(x, y, dot, stroke=0, fill=1)
    text(x + label_offset, y - 2.0, method, size=6.0)
pdf.showPage()
pdf.save()
