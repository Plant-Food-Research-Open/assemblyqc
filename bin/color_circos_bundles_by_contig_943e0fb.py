#!/usr/bin/env python

import sys
import re

bundled_links_file_name = sys.argv[1]


def natural_key(string):
    """Return a list of keys that sort naturally."""
    return [int(s) if s.isdigit() else s for s in re.split(r"(\d+)", string)]


def hsv2rgb(h, s, v):
    """Convert HSV color to RGB color."""
    h = float(h)
    s = float(s)
    v = float(v)
    h60 = h / 60.0
    h60f = int(h60)
    hi = int(h60f) % 6
    f = h60 - h60f
    p = v * (1 - s)
    q = v * (1 - f * s)
    t = v * (1 - (1 - f) * s)
    r, g, b = 0, 0, 0
    if hi == 0:
        r, g, b = v, t, p
    elif hi == 1:
        r, g, b = q, v, p
    elif hi == 2:
        r, g, b = p, v, t
    elif hi == 3:
        r, g, b = p, q, v
    elif hi == 4:
        r, g, b = t, p, v
    elif hi == 5:
        r, g, b = v, p, q
    return int(r * 255), int(g * 255), int(b * 255)


def generate_colors(num_colors):
    """Generate a list of colors"""
    hue_step = int(360 / num_colors)
    hue = 0
    colors = []
    for i in range(num_colors):
        red, green, blue = hsv2rgb(hue, 0.8, 0.8)
        colors.append(f"{red},{green},{blue},0.5")
        hue += hue_step
    return colors


def read_file_lines(file_path):
    with open(file_path, "r") as f:
        return f.readlines()


def generate_colors_by_ids(bundle_file_lines):
    """Create a dictionary to map unique target ids to colors"""
    unique_ids = set(
        line.split()[3] for line in bundle_file_lines
    )  # index 3: Target ids
    num_unique_ids = len(unique_ids)
    colors = generate_colors(num_unique_ids)
    return dict(zip(sorted(unique_ids, key=natural_key), colors))


if __name__ == "__main__":
    bundle_file_lines = read_file_lines(bundled_links_file_name)
    id_to_color = generate_colors_by_ids(bundle_file_lines)

    for line in bundle_file_lines:
        parts = line.strip().split()
        unique_id = parts[3]  # index 3: Target ids
        color = id_to_color[unique_id]
        print(" ".join(parts[0:6] + [f"color=({color})", parts[6]]))
