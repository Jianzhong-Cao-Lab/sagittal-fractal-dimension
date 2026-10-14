"""Extract 2D and 3D box-counting fractal-dimension features from NIfTI masks."""

from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from tqdm.auto import tqdm


DEFAULT_SCALES = (1, 2, 4, 8, 16)


def box_count(coords: np.ndarray, box_size: int, rng: np.random.Generator, n_shifts: int) -> float:
    counts = []
    for _ in range(n_shifts):
        shift = rng.uniform(-box_size, 0, size=coords.shape[1])
        boxes = np.floor((coords + shift) / box_size).astype(np.int32)
        counts.append(np.unique(boxes, axis=0).shape[0])
    return float(np.mean(counts))


def fractal_dimension(
    coords: np.ndarray,
    shape: tuple[int, ...],
    rng: np.random.Generator,
    n_shifts: int = 20,
    scales: tuple[int, ...] = DEFAULT_SCALES,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    if coords.size == 0:
        return np.nan, np.nan, np.array([]), np.array([])

    valid_scales = np.asarray([scale for scale in scales if scale < min(shape) // 2])
    if valid_scales.size < 2:
        return np.nan, np.nan, valid_scales, np.array([])

    counts = np.asarray(
        [box_count(coords, int(scale), rng, n_shifts) for scale in valid_scales]
    )
    positive = counts > 0
    valid_scales = valid_scales[positive]
    counts = counts[positive]
    if counts.size < 2:
        return np.nan, np.nan, valid_scales, counts

    log_inverse_scale = np.log(1.0 / valid_scales)
    log_count = np.log(counts)
    slope, intercept = np.polyfit(log_inverse_scale, log_count, 1)
    fitted = slope * log_inverse_scale + intercept
    residual_sum = np.sum((log_count - fitted) ** 2)
    total_sum = np.sum((log_count - np.mean(log_count)) ** 2)
    r_squared = 1 - residual_sum / total_sum if total_sum > 0 else np.nan
    return float(slope), float(r_squared), valid_scales, counts


def plane_features(
    volume: np.ndarray,
    axis: int,
    rng: np.random.Generator,
    n_repeats: int = 5,
    n_shifts: int = 20,
) -> tuple[float, float, float, float, int]:
    values = []
    for index in range(volume.shape[axis]):
        image = np.take(volume, index, axis=axis)
        coords = np.argwhere(image != 0)
        candidates = [
            fractal_dimension(coords, image.shape, rng, n_shifts)[:2]
            for _ in range(n_repeats)
        ]
        candidates = [(fd, r2) for fd, r2 in candidates if np.isfinite(fd) and np.isfinite(r2)]
        if candidates:
            values.append(max(candidates, key=lambda item: item[1])[0])

    if not values:
        return np.nan, np.nan, np.nan, np.nan, 0
    values_array = np.asarray(values)
    return (
        float(values_array.mean()),
        float(values_array.max()),
        float(values_array.min()),
        float(np.median(values_array)),
        len(values),
    )


def extract_mask_features(
    nii_path: Path,
    rng: np.random.Generator,
    n_repeats: int = 5,
    n_shifts: int = 20,
) -> dict[str, float | int | str]:
    image = nib.as_closest_canonical(nib.load(str(nii_path)))
    volume = np.asarray(image.dataobj)
    patient_id = nii_path.name.removesuffix(".nii.gz")

    features: dict[str, float | int | str] = {"patient_id": patient_id}
    for plane, axis in {"sagittal": 0, "coronal": 1, "axial": 2}.items():
        mean_fd, max_fd, min_fd, median_fd, valid_slices = plane_features(
            volume, axis, rng, n_repeats, n_shifts
        )
        features.update(
            {
                f"{plane}_mean_2d_fd": mean_fd,
                f"{plane}_max_2d_fd": max_fd,
                f"{plane}_min_2d_fd": min_fd,
                f"{plane}_median_2d_fd": median_fd,
                f"{plane}_valid_slices": valid_slices,
            }
        )

    coords_3d = np.argwhere(volume != 0)
    candidates_3d = [
        fractal_dimension(coords_3d, volume.shape, rng, n_shifts)[:2]
        for _ in range(n_repeats)
    ]
    candidates_3d = [(fd, r2) for fd, r2 in candidates_3d if np.isfinite(fd) and np.isfinite(r2)]
    features["global_3d_fd"] = (
        max(candidates_3d, key=lambda item: item[1])[0] if candidates_3d else np.nan
    )
    return features


def batch_process(
    input_dir: Path,
    output_csv: Path,
    n_repeats: int = 5,
    n_shifts: int = 20,
    seed: int = 123,
) -> pd.DataFrame:
    files = sorted(input_dir.glob("*.nii.gz"))
    if not files:
        raise FileNotFoundError(f"No .nii.gz files found in {input_dir}")

    rng = np.random.default_rng(seed)
    rows = []
    for path in tqdm(files, desc="Extracting fractal-dimension features"):
        try:
            rows.append(extract_mask_features(path, rng, n_repeats, n_shifts))
        except Exception as error:
            tqdm.write(f"Skipped {path.name}: {error}")
            rows.append({"patient_id": path.name.removesuffix(".nii.gz")})

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(output_csv, index=False, encoding="utf-8")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True, help="Directory containing NIfTI masks")
    parser.add_argument("--output-csv", type=Path, required=True, help="Destination feature table")
    parser.add_argument("--repeats", type=int, default=5, help="Repeated estimates per slice or volume")
    parser.add_argument("--shifts", type=int, default=20, help="Random grid shifts per box size")
    parser.add_argument("--seed", type=int, default=123, help="Random seed")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    batch_process(args.input_dir, args.output_csv, args.repeats, args.shifts, args.seed)


if __name__ == "__main__":
    main()
