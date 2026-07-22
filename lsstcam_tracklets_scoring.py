from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from sklearn.metrics import precision_recall_curve, roc_auc_score, roc_curve
from sorcha.modules.PPFootprintFilter import Footprint, radec_to_focal_plane, radec_to_tangent_plane


POP_COLORS = {
    "NEO": "#d32f2f",
    "MBA": "#1976d2",
    "TNO": "#e6a700",
    "Trojan": "#2e8b57",
    "other": "#6b7280",
    "unmatched": "#a0a0a0",
}


@dataclass
class AnalysisData:
    night: int
    root: Path
    figure_dir: Path
    footprint: Footprint
    visits: pd.DataFrame
    sorcha: pd.DataFrame
    benchmark: pd.DataFrame
    matched: pd.DataFrame
    active_area_deg2: float
    camera_width_deg: float


def _nearest_visit_indices(times: np.ndarray, visits: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    midpoints = visits["observationStartMJD"].to_numpy() + visits["visitTime"].to_numpy() / 172800.0
    insertion = np.clip(np.searchsorted(midpoints, times), 1, len(midpoints) - 1)
    candidates = np.column_stack((insertion - 1, insertion))
    distances = np.abs(times[:, None] - midpoints[candidates])
    best = candidates[np.arange(len(times)), np.argmin(distances, axis=1)]
    return best, np.abs(times - midpoints[best]) * 86400.0


def _attach_visits(tracklets: pd.DataFrame, visits: pd.DataFrame) -> pd.DataFrame:
    out = tracklets.copy().reset_index(drop=True)
    for detection in (0, 1):
        best, error_s = _nearest_visit_indices(out[f"mjd{detection}_tai"].to_numpy(), visits)
        for column in visits.columns:
            out[f"{column}{detection}"] = visits[column].to_numpy()[best]
        out[f"visit_match_error_s{detection}"] = error_s

        db_band = out[f"band{detection}"].astype(str).to_numpy()
        tracklet_band = out[f"filter{detection}"].astype(str).to_numpy()
        if not np.array_equal(db_band, tracklet_band):
            raise ValueError(f"Detection {detection} visit matches do not preserve the observing band.")
        if error_s.max() > 0.01:
            raise ValueError(
                f"Detection {detection} has a visit midpoint mismatch of {error_s.max():.6f} s."
            )

    out["field_key"] = (
        out["observationId0"].astype(int).astype(str)
        + "-"
        + out["observationId1"].astype(int).astype(str)
    )
    return out


def _validate_detectors(tracklets: pd.DataFrame, footprint: Footprint) -> pd.DataFrame:
    out = tracklets.copy()
    for detection in (0, 1):
        frame = pd.DataFrame(
            {
                "RA_deg": out[f"ra{detection}"],
                "Dec_deg": out[f"dec{detection}"],
                "FieldID": out[f"observationId{detection}"],
                "fieldRA_deg": out[f"fieldRA{detection}"],
                "fieldDec_deg": out[f"fieldDec{detection}"],
                "fieldRotSkyPos_deg": out[f"rotSkyPos{detection}"],
            }
        )
        indices, detector_ids = footprint.applyFootprint(frame, edge_thresh=2.0)
        if len(indices) != len(np.unique(indices)):
            raise ValueError("A detection was assigned to more than one LSSTCam detector.")
        on_sensor = np.zeros(len(out), dtype=bool)
        detector = np.full(len(out), -1, dtype=int)
        on_sensor[indices] = True
        detector[indices] = detector_ids
        out[f"on_sensor{detection}"] = on_sensor
        out[f"detector{detection}"] = detector
    return out


def _camera_metrics(footprint: Footprint) -> tuple[float, float]:
    area = 0.0
    all_x = []
    for detector in footprint.detectors:
        x = np.degrees(detector.x)
        y = np.degrees(detector.y)
        area += 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
        all_x.extend(x)
    width = float(np.max(all_x) - np.min(all_x))
    return float(area), width


def prepare_analysis(root: str | Path, night: int = 61642) -> AnalysisData:
    root = Path(root).resolve()
    figure_dir = root / "Figures"
    figure_dir.mkdir(exist_ok=True)

    pointing_db = root.parent / "baseline_v5.0.0_2yrs.db"
    if not pointing_db.exists():
        raise FileNotFoundError(f"Pointing database not found: {pointing_db}")

    sorcha_path = root / "sorcha_comparison_case1_nbody_Vband.parquet"
    benchmark_path = (
        root
        / "outputs"
        / "phase2_benchmark_s3m_nbody_mag245"
        / "benchmark_comparison_s3m_nbody_mag245.parquet"
    )
    sorcha = pd.read_parquet(sorcha_path)
    sorcha = sorcha.loc[sorcha["night"].eq(night)].copy().reset_index(drop=True)
    benchmark = pd.read_parquet(benchmark_path)

    query = """
        SELECT observationId, observationStartMJD, visitExposureTime, visitTime,
               fieldRA, fieldDec, rotSkyPos, band
        FROM observations
        WHERE observationStartMJD >= ? AND observationStartMJD < ?
        ORDER BY observationStartMJD
    """
    with sqlite3.connect(pointing_db) as connection:
        visits = pd.read_sql_query(query, connection, params=(night - 0.1, night + 1.1))

    sorcha = _attach_visits(sorcha, visits)
    footprint = Footprint()
    sorcha = _validate_detectors(sorcha, footprint)
    if not (sorcha["on_sensor0"].all() and sorcha["on_sensor1"].all()):
        raise ValueError("At least one retained Sorcha detection failed the exact LSSTCam footprint.")

    truth = benchmark.set_index("s3m_objid")
    matched = sorcha.join(
        truth[["population", "P_NEO_vdp", "P_NEO_d2", "mean_mag", "vlam", "vbeta"]].rename(
            columns={
                "population": "population_bench",
                "P_NEO_vdp": "P_bench_vdp",
                "P_NEO_d2": "P_bench_d2",
                "mean_mag": "mean_mag_bench",
                "vlam": "vlam_bench",
                "vbeta": "vbeta_bench",
            }
        ),
        on="ObjID",
        how="inner",
    )
    matched = matched.rename(columns={"population": "population_sorcha"})
    matched["dP_vdp"] = matched["P_NEO_vdp_Vband"] - matched["P_bench_vdp"]
    area, width = _camera_metrics(footprint)

    return AnalysisData(
        night=night,
        root=root,
        figure_dir=figure_dir,
        footprint=footprint,
        visits=visits,
        sorcha=sorcha,
        benchmark=benchmark,
        matched=matched,
        active_area_deg2=area,
        camera_width_deg=width,
    )


def analysis_audit(data: AnalysisData) -> pd.Series:
    return pd.Series(
        {
            "night": data.night,
            "pointing visits loaded": len(data.visits),
            "Sorcha tracklets": len(data.sorcha),
            "identity matches": len(data.matched),
            "matched NEOs": int(data.matched["population_bench"].eq("NEO").sum()),
            "unique visit pairs": data.sorcha["field_key"].nunique(),
            "maximum visit-match error (s)": max(
                data.sorcha["visit_match_error_s0"].max(),
                data.sorcha["visit_match_error_s1"].max(),
            ),
            "detections on active silicon": int(
                data.sorcha[["on_sensor0", "on_sensor1"]].to_numpy().sum()
            ),
            "active detectors": data.footprint.N,
            "active focal-plane area (deg^2, tangent plane)": data.active_area_deg2,
            "camera width (deg, tangent plane)": data.camera_width_deg,
        },
        name="value",
    )


def plot_camera_geometry(data: AnalysisData):
    fig, ax = plt.subplots(figsize=(7.5, 7.2))
    for detector in data.footprint.detectors:
        points = np.column_stack((np.degrees(detector.x), np.degrees(detector.y)))
        ax.add_patch(Polygon(points, facecolor="#d8e8f5", edgecolor="#315f7d", lw=0.32))
    ax.axhline(0, color="0.4", lw=0.5)
    ax.axvline(0, color="0.4", lw=0.5)
    ax.set_aspect("equal")
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.9, 1.9)
    ax.set_xlabel("Focal-plane x (deg)")
    ax.set_ylabel("Focal-plane y (deg)")
    ax.set_title(
        f"Sorcha built-in LSSTCam footprint: {data.footprint.N} detectors\n"
        f"active area={data.active_area_deg2:.2f} deg$^2$, width={data.camera_width_deg:.2f} deg"
    )
    ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(data.figure_dir / "LSSTcam_01_geometry.png", dpi=180, bbox_inches="tight")
    return fig


def _angular_separation(ra, dec, ras, decs):
    ra, dec, ras, decs = map(np.radians, (ra, dec, ras, decs))
    cosine = np.sin(dec) * np.sin(decs) + np.cos(dec) * np.cos(decs) * np.cos(ra - ras)
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def _field_label(index: int) -> str:
    label = ""
    value = index + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        label = chr(65 + remainder) + label
    return label


def select_dense_fields(
    data: AnalysisData,
    min_center_separation_deg: float = 3.5,
    target_neo: int = 75,
    target_nonneo: int = 75,
    max_fields: int = 40,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fields = (
        data.sorcha.groupby(["field_key", "observationId0", "observationId1"], as_index=False)
        .agg(
            n_tracklets=("tracklet_id", "size"),
            n_objects=("ObjID", "nunique"),
            fieldRA=("fieldRA0", "first"),
            fieldDec=("fieldDec0", "first"),
            rotSkyPos0=("rotSkyPos0", "first"),
            rotSkyPos1=("rotSkyPos1", "first"),
            mjd0=("observationStartMJD0", "first"),
            mjd1=("observationStartMJD1", "first"),
            band0=("band0", "first"),
            band1=("band1", "first"),
        )
        .sort_values(["n_tracklets", "mjd0", "observationId0"], ascending=[False, True, True])
        .reset_index(drop=True)
    )
    matched_counts = data.matched.groupby("field_key")["population_bench"].agg(
        n_matched="size", n_neo=lambda values: int(values.eq("NEO").sum())
    )
    fields = fields.join(matched_counts, on="field_key")
    fields[["n_matched", "n_neo"]] = fields[["n_matched", "n_neo"]].fillna(0).astype(int)
    fields["n_nonneo"] = fields["n_matched"] - fields["n_neo"]

    selected_indices = []
    cumulative_neo = 0
    cumulative_nonneo = 0
    for row in fields.itertuples():
        if selected_indices:
            prior = fields.loc[selected_indices]
            separation = _angular_separation(
                row.fieldRA,
                row.fieldDec,
                prior["fieldRA"].to_numpy(),
                prior["fieldDec"].to_numpy(),
            )
            if np.any(separation < min_center_separation_deg):
                continue
        selected_indices.append(row.Index)
        cumulative_neo += row.n_neo
        cumulative_nonneo += row.n_nonneo
        if cumulative_neo >= target_neo and cumulative_nonneo >= target_nonneo:
            break
        if len(selected_indices) >= max_fields:
            break

    selected = fields.loc[selected_indices].copy().reset_index(drop=True)
    selected["field_label"] = [_field_label(i) for i in range(len(selected))]
    selected["cum_neo"] = selected["n_neo"].cumsum()
    selected["cum_nonneo"] = selected["n_nonneo"].cumsum()
    selected["cum_matched"] = selected["n_matched"].cumsum()
    selected["cum_tracklets"] = selected["n_tracklets"].cumsum()

    label_map = selected.set_index("field_key")["field_label"]
    selected_tracklets = data.sorcha.loc[data.sorcha["field_key"].isin(label_map.index)].copy()
    selected_tracklets["field_label"] = selected_tracklets["field_key"].map(label_map)
    selected_matched = data.matched.loc[data.matched["field_key"].isin(label_map.index)].copy()
    selected_matched["field_label"] = selected_matched["field_key"].map(label_map)

    fields.head(30).to_csv(data.root / "outputs" / "LSSTcam_top30_field_candidates.csv", index=False)
    selected.to_csv(data.root / "outputs" / "LSSTcam_selected_fields.csv", index=False)
    selected_tracklets.to_csv(data.root / "outputs" / "LSSTcam_selected_tracklets.csv", index=False)
    return selected, selected_tracklets, selected_matched


def _focal_to_radec(x, y, field_ra_deg, field_dec_deg, rotation_deg):
    local = (np.asarray(x) + 1j * np.asarray(y)) * np.exp(-1j * np.radians(rotation_deg))
    tx = np.real(local)
    ty = np.imag(local)
    ra0 = np.radians(field_ra_deg)
    dec0 = np.radians(field_dec_deg)
    center = np.array([np.cos(ra0) * np.cos(dec0), np.sin(ra0) * np.cos(dec0), np.sin(dec0)])
    east = np.array([-np.sin(ra0), np.cos(ra0), 0.0])
    north = np.cross(center, east)
    vectors = center[:, None] + east[:, None] * tx + north[:, None] * ty
    vectors /= np.linalg.norm(vectors, axis=0)
    ra = np.mod(np.degrees(np.arctan2(vectors[1], vectors[0])), 360.0)
    dec = np.degrees(np.arcsin(vectors[2]))
    return ra, dec


def _camera_polygons_in_reference(footprint, visit, reference_ra, reference_dec):
    polygons = []
    for detector in footprint.detectors:
        ra, dec = _focal_to_radec(
            detector.x,
            detector.y,
            visit.fieldRA,
            visit.fieldDec,
            visit.rotSkyPos,
        )
        x, y = radec_to_tangent_plane(
            np.radians(ra),
            np.radians(dec),
            np.full(len(ra), np.radians(reference_ra)),
            np.full(len(dec), np.radians(reference_dec)),
        )
        polygons.append(np.column_stack((np.degrees(x), np.degrees(y))))
    return polygons


def _tracklet_segments(frame, reference_ra, reference_dec):
    reference_ra_array = np.full(len(frame), np.radians(reference_ra))
    reference_dec_array = np.full(len(frame), np.radians(reference_dec))
    x0, y0 = radec_to_tangent_plane(
        np.radians(frame["ra0"].to_numpy()),
        np.radians(frame["dec0"].to_numpy()),
        reference_ra_array,
        reference_dec_array,
    )
    x1, y1 = radec_to_tangent_plane(
        np.radians(frame["ra1"].to_numpy()),
        np.radians(frame["dec1"].to_numpy()),
        reference_ra_array,
        reference_dec_array,
    )
    return np.stack(
        [np.column_stack((np.degrees(x0), np.degrees(y0))), np.column_stack((np.degrees(x1), np.degrees(y1)))],
        axis=1,
    )


def plot_field_a(data: AnalysisData, selected: pd.DataFrame):
    field_a = selected.iloc[0]
    tracklets = data.sorcha.loc[data.sorcha["field_key"].eq(field_a.field_key)].copy()
    matched_population = data.matched.set_index("ObjID")["population_bench"]
    tracklets["population_plot"] = tracklets["ObjID"].map(matched_population).fillna("unmatched")
    segments = _tracklet_segments(tracklets, field_a.fieldRA, field_a.fieldDec)

    visit0 = pd.Series(
        {"fieldRA": tracklets["fieldRA0"].iloc[0], "fieldDec": tracklets["fieldDec0"].iloc[0], "rotSkyPos": tracklets["rotSkyPos0"].iloc[0]}
    )
    visit1 = pd.Series(
        {"fieldRA": tracklets["fieldRA1"].iloc[0], "fieldDec": tracklets["fieldDec1"].iloc[0], "rotSkyPos": tracklets["rotSkyPos1"].iloc[0]}
    )
    polygons0 = _camera_polygons_in_reference(data.footprint, visit0, field_a.fieldRA, field_a.fieldDec)
    polygons1 = _camera_polygons_in_reference(data.footprint, visit1, field_a.fieldRA, field_a.fieldDec)

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.4), sharex=True, sharey=True)
    for ax in axes:
        for points in polygons0:
            ax.add_patch(Polygon(points, fill=False, edgecolor="#2171b5", lw=0.22, alpha=0.55))
        for points in polygons1:
            ax.add_patch(Polygon(points, fill=False, edgecolor="#d95f0e", lw=0.22, alpha=0.45, ls="--"))
        ax.set_aspect("equal")
        ax.set_xlim(-1.9, 1.9)
        ax.set_ylim(-1.9, 1.9)
        ax.set_xlabel("East offset from field center (deg)")
        ax.grid(alpha=0.12)
    axes[0].set_ylabel("North offset from field center (deg)")

    for population in ["NEO", "MBA", "TNO", "other", "unmatched"]:
        mask = tracklets["population_plot"].eq(population).to_numpy()
        if mask.any():
            axes[0].add_collection(
                LineCollection(segments[mask], colors=POP_COLORS[population], linewidths=1.6, alpha=0.9)
            )
            axes[0].scatter(
                segments[mask, 0, 0], segments[mask, 0, 1], s=18,
                color=POP_COLORS[population], label=f"{population} (n={mask.sum()})", zorder=3,
            )
    axes[0].legend(fontsize=8, loc="lower left")
    axes[0].set_title("Tracklets by benchmark truth; gray = unmatched")

    score = tracklets["P_NEO_vdp_Vband"].to_numpy()
    collection = LineCollection(segments, cmap="viridis", linewidths=2.0, alpha=0.95)
    collection.set_array(score)
    collection.set_clim(0, 1)
    axes[1].add_collection(collection)
    axes[1].scatter(segments[:, 0, 0], segments[:, 0, 1], c=score, cmap="viridis", vmin=0, vmax=1, s=18)
    fig.colorbar(collection, ax=axes[1], label="Sorcha V-band VDP score")
    axes[1].set_title("Same tracklets colored by VDP score")

    fig.suptitle(
        f"Field A: LSSTCam visit pair {int(field_a.observationId0)} / {int(field_a.observationId1)}; "
        f"n={int(field_a.n_tracklets)}\n"
        f"center=({field_a.fieldRA:.3f} deg, {field_a.fieldDec:.3f} deg), "
        f"rotations={field_a.rotSkyPos0:.2f} deg / {field_a.rotSkyPos1:.2f} deg"
    )
    fig.tight_layout()
    fig.savefig(data.figure_dir / "LSSTcam_02_field_A_tracklets.png", dpi=180, bbox_inches="tight")
    return fig


def plot_selected_field_sky(data: AnalysisData, selected: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(13, 6.5))
    for row in selected.itertuples():
        visit = pd.Series({"fieldRA": row.fieldRA, "fieldDec": row.fieldDec, "rotSkyPos": row.rotSkyPos0})
        for detector in data.footprint.detectors:
            ra, dec = _focal_to_radec(
                detector.x, detector.y, visit.fieldRA, visit.fieldDec, visit.rotSkyPos
            )
            wrapped_ra = (ra + 180.0) % 360.0 - 180.0
            ax.plot(wrapped_ra, dec, color="#3b82a0", lw=0.12, alpha=0.35)
        center_ra = (row.fieldRA + 180.0) % 360.0 - 180.0
        ax.text(center_ra, row.fieldDec, row.field_label, ha="center", va="center", fontsize=7, weight="bold")
    all_ra = (data.sorcha["mean_ra"].to_numpy() + 180.0) % 360.0 - 180.0
    ax.scatter(all_ra, data.sorcha["mean_dec"], s=2, color="0.75", alpha=0.18, label="all night-61642 tracklets")
    ax.set_xlim(180, -180)
    ax.set_ylim(-65, 20)
    ax.set_xlabel("Right ascension (deg; wrapped)")
    ax.set_ylabel("Declination (deg)")
    ax.set_title(
        f"Density-selected, spatially distinct LSSTCam fields (n={len(selected)})\n"
        f"target sample: {selected.n_neo.sum()} matched NEOs, {selected.n_nonneo.sum()} matched non-NEOs"
    )
    ax.grid(alpha=0.16)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(data.figure_dir / "LSSTcam_03_selected_fields.png", dpi=180, bbox_inches="tight")
    return fig


def _pr_values(y, score):
    precision, recall, thresholds = precision_recall_curve(y, score)
    f1 = np.divide(2 * precision * recall, precision + recall, out=np.zeros_like(precision), where=(precision + recall) > 0)
    best = int(np.argmax(f1[:-1]))
    return recall * 100.0, (1.0 - precision) * 100.0, float(f1[best]), float(thresholds[best])


def _cluster_pr_band(frame, score_column, field_column="field_key", repeats=400, seed=61642):
    rng = np.random.default_rng(seed)
    keys = frame[field_column].unique()
    groups = {key: frame.loc[frame[field_column].eq(key)] for key in keys}
    recall_grid = np.linspace(0, 100, 101)
    curves = []
    for _ in range(repeats):
        sample = pd.concat([groups[key] for key in rng.choice(keys, size=len(keys), replace=True)], ignore_index=True)
        y = sample["population_bench"].eq("NEO").astype(int).to_numpy()
        if y.min() == y.max():
            continue
        recall, contamination, _, _ = _pr_values(y, sample[score_column].to_numpy())
        order = np.argsort(recall)
        curves.append(np.interp(recall_grid, recall[order], contamination[order]))
    low, high = np.nanpercentile(np.asarray(curves), [16, 84], axis=0)
    return recall_grid, low, high


def plot_dense_precision_recall(data: AnalysisData, selected_matched: pd.DataFrame):
    y = selected_matched["population_bench"].eq("NEO").astype(int).to_numpy()
    panels = [
        ("Benchmark", [("P_bench_vdp", "VDP", "#1976d2"), ("P_bench_d2", "digest2", "#2e8b57")]),
        ("Sorcha", [("P_NEO_vdp_Vband", "VDP", "#d32f2f"), ("P_NEO_d2", "digest2", "#e07a1f")]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), sharex=True, sharey=True)
    for ax, (title, series) in zip(axes, panels):
        for column, label, color in series:
            completeness, contamination, f1, threshold = _pr_values(y, selected_matched[column].to_numpy())
            ax.plot(completeness, contamination, color=color, lw=2.2, label=f"{label}: F1={f1:.3f}, t={threshold:.3f}")
            grid, low, high = _cluster_pr_band(selected_matched, column)
            ax.fill_between(grid, low, high, color=color, alpha=0.13)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.set_xlabel("NEO completeness (%)")
        ax.set_title(title)
        ax.grid(alpha=0.22)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Contamination (%)")
    fig.suptitle(
        f"Dense LSSTCam fields: VDP vs digest2, same {len(selected_matched)} identity-matched objects\n"
        "bands are 16-84% field-cluster bootstrap intervals"
    )
    fig.tight_layout()
    fig.savefig(data.figure_dir / "LSSTcam_04_dense_fields_precision_recall.png", dpi=180, bbox_inches="tight")
    return fig


def plot_dense_roc(data: AnalysisData, selected_matched: pd.DataFrame):
    y = selected_matched["population_bench"].eq("NEO").astype(int).to_numpy()
    panels = [
        ("Benchmark", [("P_bench_vdp", "VDP", "#1976d2"), ("P_bench_d2", "digest2", "#2e8b57")]),
        ("Sorcha", [("P_NEO_vdp_Vband", "VDP", "#d32f2f"), ("P_NEO_d2", "digest2", "#e07a1f")]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), sharex=True, sharey=True)
    for ax, (title, series) in zip(axes, panels):
        for column, label, color in series:
            fpr, tpr, _ = roc_curve(y, selected_matched[column].to_numpy())
            auc = roc_auc_score(y, selected_matched[column].to_numpy())
            ax.plot(fpr, tpr, lw=2.2, color=color, label=f"{label}: AUC={auc:.3f}")
        ax.plot([0, 1], [0, 1], color="0.5", lw=1, ls="--")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("False-positive rate")
        ax.set_title(title)
        ax.grid(alpha=0.22)
        ax.legend(fontsize=9, loc="lower right")
    axes[0].set_ylabel("True-positive rate")
    fig.suptitle("Conventional ROC for the selected dense LSSTCam fields")
    fig.tight_layout()
    fig.savefig(data.figure_dir / "LSSTcam_05_dense_fields_roc.png", dpi=180, bbox_inches="tight")
    return fig


def plot_zero_diagnostics(data: AnalysisData, selected_matched: pd.DataFrame):
    zero = selected_matched["P_NEO_vdp_Vband"].eq(0)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
    for population, frame in selected_matched.groupby("population_bench"):
        axes[0].scatter(
            frame["P_bench_vdp"], frame["P_NEO_vdp_Vband"], s=22,
            alpha=0.7, color=POP_COLORS.get(population, POP_COLORS["other"]), label=f"{population} (n={len(frame)})",
        )
    axes[0].plot([0, 1], [0, 1], color="0.25", ls="--", lw=1)
    axes[0].set_xlim(-0.02, 1.02)
    axes[0].set_ylim(-0.02, 1.02)
    axes[0].set_xlabel("Benchmark VDP score")
    axes[0].set_ylabel("Sorcha V-band VDP score")
    axes[0].set_title("Paired VDP scores")
    axes[0].grid(alpha=0.2)
    axes[0].legend(fontsize=8)

    zero_frame = selected_matched.loc[zero]
    if len(zero_frame):
        axes[1].hist(
            [
                zero_frame.loc[zero_frame["population_bench"].eq("NEO"), "P_NEO_d2"],
                zero_frame.loc[~zero_frame["population_bench"].eq("NEO"), "P_NEO_d2"],
            ],
            bins=np.linspace(0, 1, 21), stacked=True,
            color=[POP_COLORS["NEO"], POP_COLORS["MBA"]], label=["NEO", "non-NEO"], alpha=0.85,
        )
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("Sorcha digest2 score")
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"Objects assigned Sorcha VDP = 0 (n={zero.sum()})")
    axes[1].grid(alpha=0.2)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(data.figure_dir / "LSSTcam_06_zero_diagnostics.png", dpi=180, bbox_inches="tight")
    return fig
