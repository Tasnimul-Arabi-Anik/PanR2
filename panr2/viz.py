import textwrap


PLOT_STYLE_CONFIGS = {
    "publication": {
        "dpi": 300,
        "label_max_length": 55,
        "wrap_width": 18,
        "heatmap_max_features": 50,
        "network_max_edges": 75,
    },
    "dashboard": {
        "dpi": 180,
        "label_max_length": 75,
        "wrap_width": 24,
        "heatmap_max_features": 80,
        "network_max_edges": 120,
    },
    "compact": {
        "dpi": 180,
        "label_max_length": 38,
        "wrap_width": 14,
        "heatmap_max_features": 30,
        "network_max_edges": 45,
    },
}


def plot_style_config(style="publication", label_max_length=None):
    """Return deterministic plotting limits for a named PanR2 plot style."""
    config = dict(PLOT_STYLE_CONFIGS.get(style or "publication", PLOT_STYLE_CONFIGS["publication"]))
    config["style"] = style or "publication"
    if label_max_length is not None:
        config["label_max_length"] = int(label_max_length)
    return config


def shorten_label(value, max_length=55, wrap_width=18):
    """Shorten and wrap labels for crowded static/HTML figures."""
    text = str(value)
    if max_length and len(text) > max_length:
        text = text[: max(1, max_length - 1)] + "..."
    if wrap_width and len(text) > wrap_width:
        return "\n".join(textwrap.wrap(text, width=wrap_width, break_long_words=False) or [text])
    return text


def label_map(labels, max_length=55, wrap_width=18):
    return {label: shorten_label(label, max_length=max_length, wrap_width=wrap_width) for label in labels}


def label_warning_rows(labels, max_length=55, figure=""):
    rows = []
    long_count = sum(1 for label in labels if len(str(label)) > max_length)
    if long_count:
        rows.append({
            "figure": figure,
            "warning_type": "label_truncation",
            "detail": f"{long_count} label(s) exceeded {max_length} characters and were shortened for readability.",
            "value": long_count,
            "limit": max_length,
        })
    return rows


def crowding_warning_row(figure, count, limit, item="features"):
    if limit and count > limit:
        return {
            "figure": figure,
            "warning_type": "crowded_plot",
            "detail": f"{count} {item} were available; the plot was limited to {limit} for readability.",
            "value": count,
            "limit": limit,
        }
    return None
