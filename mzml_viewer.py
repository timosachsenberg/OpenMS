#!/usr/bin/env python3
"""
Fast mzML Peak Map Viewer using NiceGUI + Datashader + pyOpenMS

Designed to handle 50+ million peaks with smooth zooming and panning.
Uses datashader for server-side rendering of massive datasets.
Supports FeatureMap overlay with centroids, bounding boxes, and convex hulls.
Supports idXML overlay showing peptide identification precursor positions.
Includes annotated MS2 spectrum viewer for peptide identifications.
Displays Total Ion Chromatogram (TIC) with clickable MS1 spectrum viewer.

Usage:
    python mzml_viewer.py                           # Start with empty viewer
    python mzml_viewer.py sample.mzML               # Load mzML file
    python mzml_viewer.py sample.mzML features.featureXML  # Load mzML + features
    python mzml_viewer.py sample.mzML ids.idXML     # Load mzML + identifications
    python mzml_viewer.py sample.mzML features.featureXML ids.idXML  # All three
"""

import io
import sys
import base64
import math
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any

import click
import plotly.graph_objects as go

# Datashader for fast rendering
import datashader as ds
import datashader.transfer_functions as tf
from colorcet import fire

# PIL for drawing overlays and axes
from PIL import Image, ImageDraw, ImageFont

# pyOpenMS for file loading and spectrum annotation
from pyopenms import (
    MSExperiment, MzMLFile,
    FeatureMap, FeatureXMLFile,
    IdXMLFile, PeptideIdentification, ProteinIdentification,
    TheoreticalSpectrumGenerator, AASequence, Param,
    MSSpectrum
)

# NiceGUI for the web interface
from nicegui import ui, app


# Global viewer instance for CLI file loading
_viewer_instance = None
_cli_files = {'mzml': None, 'featurexml': None, 'idxml': None}

# Ion type colors for spectrum annotation
ION_COLORS = {
    'b': '#1f77b4',  # Blue
    'y': '#d62728',  # Red
    'a': '#2ca02c',  # Green
    'c': '#9467bd',  # Purple
    'x': '#8c564b',  # Brown
    'z': '#e377c2',  # Pink
    'precursor': '#ff7f0e',  # Orange
    'unknown': '#7f7f7f',  # Gray
}


def calculate_nice_ticks(vmin: float, vmax: float, num_ticks: int = 6) -> List[float]:
    """Calculate nice round tick values for an axis."""
    if vmin >= vmax:
        return [vmin]

    range_val = vmax - vmin
    rough_step = range_val / (num_ticks - 1)

    mag = math.floor(math.log10(rough_step))
    pow10 = 10 ** mag
    norm_step = rough_step / pow10

    if norm_step < 1.5:
        nice_step = 1
    elif norm_step < 3:
        nice_step = 2
    elif norm_step < 7:
        nice_step = 5
    else:
        nice_step = 10

    step = nice_step * pow10
    first_tick = math.ceil(vmin / step) * step
    ticks = []
    tick = first_tick
    while tick <= vmax + step * 0.001:
        ticks.append(tick)
        tick += step

    return ticks


def format_tick_label(value: float, range_val: float) -> str:
    """Format a tick label based on the value and range."""
    if range_val >= 1000:
        if abs(value) >= 1000:
            return f"{value:.0f}"
        return f"{value:.1f}"
    elif range_val >= 10:
        return f"{value:.1f}"
    elif range_val >= 1:
        return f"{value:.2f}"
    else:
        return f"{value:.3f}"


def generate_theoretical_spectrum(sequence: AASequence, charge: int) -> Dict[str, List[Tuple[float, str]]]:
    """Generate theoretical b/y ion spectrum for annotation."""
    tsg = TheoreticalSpectrumGenerator()
    spec = MSSpectrum()

    # Configure for b and y ions
    params = tsg.getParameters()
    params.setValue("add_b_ions", "true")
    params.setValue("add_y_ions", "true")
    params.setValue("add_a_ions", "false")
    params.setValue("add_c_ions", "false")
    params.setValue("add_x_ions", "false")
    params.setValue("add_z_ions", "false")
    params.setValue("add_metainfo", "true")
    tsg.setParameters(params)

    tsg.getSpectrum(spec, sequence, 1, min(charge, 2))

    ions = {'b': [], 'y': [], 'other': []}

    for i in range(spec.size()):
        mz = spec[i].getMZ()
        intensity = spec[i].getIntensity()

        # Get ion annotation from metadata
        ion_name = ""
        if spec[i].metaValueExists("IonName"):
            ion_name = spec[i].getMetaValue("IonName")

        if ion_name.startswith('b'):
            ions['b'].append((mz, ion_name))
        elif ion_name.startswith('y'):
            ions['y'].append((mz, ion_name))
        else:
            ions['other'].append((mz, ion_name))

    return ions


def create_annotated_spectrum_plot(
    exp_mz: np.ndarray,
    exp_int: np.ndarray,
    sequence_str: str,
    charge: int,
    precursor_mz: float,
    tolerance_da: float = 0.5
) -> go.Figure:
    """Create an annotated spectrum plot using Plotly."""

    # Normalize intensities to percentage
    max_int = exp_int.max() if len(exp_int) > 0 else 1
    exp_int_norm = (exp_int / max_int) * 100

    # Create figure
    fig = go.Figure()

    # Add experimental spectrum as gray bars
    fig.add_trace(go.Bar(
        x=exp_mz,
        y=exp_int_norm,
        marker_color='gray',
        name='Experimental',
        width=0.5,
        opacity=0.6,
        hovertemplate='m/z: %{x:.4f}<br>Intensity: %{y:.1f}%<extra></extra>'
    ))

    # Try to generate theoretical spectrum for annotation
    try:
        seq = AASequence.fromString(sequence_str)
        theo_ions = generate_theoretical_spectrum(seq, charge)

        # Match theoretical to experimental and annotate
        annotations = []
        matched_mz = []
        matched_int = []
        matched_labels = []
        matched_colors = []

        for ion_type, ions in [('b', theo_ions['b']), ('y', theo_ions['y'])]:
            color = ION_COLORS[ion_type]
            for theo_mz, ion_name in ions:
                # Find closest experimental peak
                if len(exp_mz) > 0:
                    diffs = np.abs(exp_mz - theo_mz)
                    min_idx = np.argmin(diffs)
                    if diffs[min_idx] <= tolerance_da:
                        matched_mz.append(exp_mz[min_idx])
                        matched_int.append(exp_int_norm[min_idx])
                        matched_labels.append(ion_name)
                        matched_colors.append(color)

        # Add matched peaks as colored bars
        if matched_mz:
            for i, (mz, intensity, label, color) in enumerate(zip(matched_mz, matched_int, matched_labels, matched_colors)):
                fig.add_trace(go.Bar(
                    x=[mz],
                    y=[intensity],
                    marker_color=color,
                    name=label if i < 10 else None,  # Only show first 10 in legend
                    showlegend=(i < 10),
                    width=1.0,
                    hovertemplate=f'{label}<br>m/z: {mz:.4f}<br>Intensity: {intensity:.1f}%<extra></extra>'
                ))

                # Add text annotation
                fig.add_annotation(
                    x=mz,
                    y=intensity + 3,
                    text=label,
                    showarrow=False,
                    font=dict(size=9, color=color),
                    textangle=-45
                )

    except Exception as e:
        # If annotation fails, just show the raw spectrum
        pass

    # Add precursor marker
    fig.add_vline(x=precursor_mz, line_dash="dash", line_color="orange",
                  annotation_text=f"Precursor ({precursor_mz:.2f})")

    # Update layout
    fig.update_layout(
        title=dict(
            text=f"MS2 Spectrum: {sequence_str} (z={charge}+)",
            font=dict(size=14)
        ),
        xaxis_title="m/z",
        yaxis_title="Relative Intensity (%)",
        template="plotly_dark",
        height=400,
        margin=dict(l=60, r=20, t=50, b=50),
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
            font=dict(size=10)
        ),
        barmode='overlay'
    )

    fig.update_xaxes(range=[0, max(exp_mz) * 1.05] if len(exp_mz) > 0 else [0, 2000])
    fig.update_yaxes(range=[0, 110])

    return fig


class MzMLViewer:
    """High-performance mzML peak map viewer using datashader with feature and ID overlay."""

    def __init__(self):
        self.exp = None
        self.df = None
        self.current_file = None

        # FeatureMap data
        self.feature_map = None
        self.features_file = None
        self.feature_data = []

        # Identification data
        self.peptide_ids = []
        self.protein_ids = []
        self.id_file = None
        self.id_data = []

        # TIC data
        self.tic_rt = None
        self.tic_intensity = None

        # FAIMS data
        self.faims_cvs = []  # List of unique CV values
        self.faims_data = {}  # Dict: CV -> DataFrame of peaks
        self.faims_tic = {}  # Dict: CV -> (rt_array, intensity_array)
        self.has_faims = False
        self.show_faims_view = False  # Toggle for FAIMS multi-panel view

        # Spectrum browser data
        self.spectrum_data = []  # List of spectrum metadata for table
        self.selected_spectrum_idx = None

        # View bounds
        self.rt_min = 0
        self.rt_max = 1
        self.mz_min = 0
        self.mz_max = 1

        # Current view
        self.view_rt_min = None
        self.view_rt_max = None
        self.view_mz_min = None
        self.view_mz_max = None

        # Selected indices
        self.selected_feature_idx = None
        self.selected_id_idx = None

        # Image dimensions
        self.plot_width = 1100
        self.plot_height = 550

        # Margins
        self.margin_left = 80
        self.margin_right = 20
        self.margin_top = 20
        self.margin_bottom = 50

        self.canvas_width = self.plot_width + self.margin_left + self.margin_right
        self.canvas_height = self.plot_height + self.margin_top + self.margin_bottom

        # Display options
        self.show_centroids = True
        self.show_bounding_boxes = True
        self.show_convex_hulls = True
        self.show_ids = True
        self.show_spectrum_marker = True  # Show RT marker for selected spectrum

        # Colors
        self.centroid_color = (0, 255, 100, 255)
        self.bbox_color = (255, 255, 0, 200)
        self.hull_color = (0, 200, 255, 150)
        self.selected_color = (255, 100, 255, 255)
        self.id_color = (255, 150, 50, 255)
        self.id_selected_color = (255, 50, 50, 255)

        self.axis_color = (200, 200, 200, 255)
        self.tick_color = (180, 180, 180, 255)
        self.label_color = (220, 220, 220, 255)
        self.grid_color = (60, 60, 60, 255)

        # UI elements
        self.image_element = None
        self.status_label = None
        self.info_label = None
        self.feature_info_label = None
        self.id_info_label = None
        self.rt_range_label = None
        self.mz_range_label = None
        self.feature_table = None
        self.id_table = None
        self.spectrum_plot = None
        self.spectrum_info_label = None
        self.tic_plot = None
        self.ms1_spectrum_plot = None
        self.ms1_spectrum_info_label = None

        # Spectrum browser UI elements
        self.spectrum_table = None
        self.spectrum_browser_plot = None
        self.spectrum_browser_info = None
        self.spectrum_nav_label = None

        # FAIMS UI elements
        self.faims_container = None  # Container for multiple peak maps
        self.faims_images = {}  # Dict: CV -> image element
        self.faims_toggle = None
        self.faims_info_label = None

        # UI update flags
        self._updating_from_tic = False  # Prevent circular TIC updates

    def _get_cv_from_spectrum(self, spec) -> Optional[float]:
        """Extract FAIMS compensation voltage from spectrum metadata."""
        # Try common CV metadata names
        cv_names = [
            "FAIMS compensation voltage",
            "ion mobility drift time",  # Alternative
            "MS:1001581",  # CV accession for FAIMS CV
        ]
        for name in cv_names:
            if spec.metaValueExists(name):
                try:
                    return float(spec.getMetaValue(name))
                except (ValueError, TypeError):
                    pass

        # Check in acquisition info / scan windows
        try:
            # Try to get from instrument settings or other metadata
            acq = spec.getAcquisitionInfo()
            if acq:
                for a in acq:
                    for name in cv_names:
                        if a.metaValueExists(name):
                            return float(a.getMetaValue(name))
        except Exception:
            pass

        return None

    def load_mzml(self, filepath: str) -> bool:
        """Load mzML file and extract peak data."""
        try:
            if self.status_label:
                self.status_label.set_text(f"Loading {Path(filepath).name}...")
            ui.notify(f"Loading {filepath}...", type="info")

            self.exp = MSExperiment()
            MzMLFile().load(filepath, self.exp)

            if self.status_label:
                self.status_label.set_text("Extracting peaks...")

            total_peaks = sum(spec.size() for spec in self.exp)

            if total_peaks == 0:
                ui.notify("No peaks found in file!", type="warning")
                return False

            # First pass: detect FAIMS CVs
            cv_set = set()
            for spec in self.exp:
                if spec.getMSLevel() == 1:
                    cv = self._get_cv_from_spectrum(spec)
                    if cv is not None:
                        cv_set.add(cv)

            self.has_faims = len(cv_set) > 1
            self.faims_cvs = sorted(cv_set) if self.has_faims else []

            # Data structures for peak extraction
            rts = np.empty(total_peaks, dtype=np.float32)
            mzs = np.empty(total_peaks, dtype=np.float32)
            intensities = np.empty(total_peaks, dtype=np.float32)
            cvs = np.empty(total_peaks, dtype=np.float32) if self.has_faims else None

            # Also compute TIC (overall and per-CV)
            tic_rts = []
            tic_intensities = []
            faims_tic_data = {cv: {'rt': [], 'int': []} for cv in self.faims_cvs} if self.has_faims else {}

            idx = 0
            for spec in self.exp:
                if spec.getMSLevel() != 1:
                    continue
                rt = spec.getRT()
                mz_array, int_array = spec.get_peaks()
                n = len(mz_array)

                cv = self._get_cv_from_spectrum(spec) if self.has_faims else None

                if n > 0:
                    rts[idx:idx+n] = rt
                    mzs[idx:idx+n] = mz_array
                    intensities[idx:idx+n] = int_array
                    if self.has_faims and cv is not None:
                        cvs[idx:idx+n] = cv
                    idx += n

                    # TIC: sum of all intensities for this spectrum
                    tic_sum = float(np.sum(int_array))
                    tic_rts.append(rt)
                    tic_intensities.append(tic_sum)

                    # Per-CV TIC
                    if self.has_faims and cv is not None:
                        faims_tic_data[cv]['rt'].append(rt)
                        faims_tic_data[cv]['int'].append(tic_sum)

            rts = rts[:idx]
            mzs = mzs[:idx]
            intensities = intensities[:idx]
            if self.has_faims:
                cvs = cvs[:idx]

            # Store TIC data
            self.tic_rt = np.array(tic_rts, dtype=np.float32)
            self.tic_intensity = np.array(tic_intensities, dtype=np.float32)

            # Store per-CV TIC data
            self.faims_tic = {}
            for cv in self.faims_cvs:
                self.faims_tic[cv] = (
                    np.array(faims_tic_data[cv]['rt'], dtype=np.float32),
                    np.array(faims_tic_data[cv]['int'], dtype=np.float32)
                )

            # Extract spectrum metadata for browser
            self.spectrum_data = self._extract_spectrum_data()

            # Create main DataFrame
            self.df = pd.DataFrame({
                'rt': rts,
                'mz': mzs,
                'intensity': intensities
            })
            if self.has_faims:
                self.df['cv'] = cvs
            self.df['log_intensity'] = np.log1p(self.df['intensity'])

            # Create per-CV DataFrames for FAIMS view
            self.faims_data = {}
            if self.has_faims:
                for cv in self.faims_cvs:
                    cv_df = self.df[self.df['cv'] == cv].copy()
                    self.faims_data[cv] = cv_df

            self.rt_min = float(self.df['rt'].min())
            self.rt_max = float(self.df['rt'].max())
            self.mz_min = float(self.df['mz'].min())
            self.mz_max = float(self.df['mz'].max())

            self.view_rt_min = self.rt_min
            self.view_rt_max = self.rt_max
            self.view_mz_min = self.mz_min
            self.view_mz_max = self.mz_max

            self.current_file = filepath

            # Build info text
            info_text = (
                f"Loaded: {Path(filepath).name} | "
                f"Spectra: {self.exp.size():,} | "
                f"Peaks: {len(self.df):,}"
            )
            if self.has_faims:
                info_text += f" | FAIMS: {len(self.faims_cvs)} CVs"

            if self.info_label:
                self.info_label.set_text(info_text)
            if self.status_label:
                self.status_label.set_text("Ready")

            # Update FAIMS UI
            if self.has_faims:
                if self.faims_info_label:
                    cv_str = ", ".join([f"{cv:.1f}V" for cv in self.faims_cvs])
                    self.faims_info_label.set_text(f"FAIMS CVs detected: {cv_str}")
                    self.faims_info_label.set_visibility(True)
                if self.faims_toggle:
                    self.faims_toggle.set_visibility(True)
                # Create FAIMS image elements
                if hasattr(self, '_create_faims_images') and self._create_faims_images:
                    self._create_faims_images()
            else:
                if self.faims_info_label:
                    self.faims_info_label.set_visibility(False)
                if self.faims_toggle:
                    self.faims_toggle.set_visibility(False)
                    self.show_faims_view = False
                if self.faims_container:
                    self.faims_container.set_visibility(False)

            # Update spectrum browser table
            if self.spectrum_table is not None:
                self.spectrum_table.update_rows(self.spectrum_data)

            ui.notify(f"Loaded {len(self.df):,} peaks", type="positive")
            if self.has_faims:
                ui.notify(f"FAIMS data detected: {len(self.faims_cvs)} compensation voltages", type="info")

            return True

        except Exception as e:
            if self.status_label:
                self.status_label.set_text(f"Error: {e}")
            ui.notify(f"Error loading file: {e}", type="negative")
            return False

    def _extract_feature_data(self) -> List[Dict[str, Any]]:
        """Extract feature data for table display."""
        if self.feature_map is None:
            return []

        data = []
        for idx, feature in enumerate(self.feature_map):
            rt = feature.getRT()
            mz = feature.getMZ()
            intensity = feature.getIntensity()
            charge = feature.getCharge()
            quality = feature.getOverallQuality()

            hulls = feature.getConvexHulls()
            rt_width = 0
            mz_width = 0
            if hulls:
                all_points = []
                for hull in hulls:
                    points = hull.getHullPoints()
                    all_points.extend([(p[0], p[1]) for p in points])
                if all_points:
                    rt_coords = [p[0] for p in all_points]
                    mz_coords = [p[1] for p in all_points]
                    rt_width = max(rt_coords) - min(rt_coords)
                    mz_width = max(mz_coords) - min(mz_coords)

            data.append({
                'idx': idx,
                'rt': round(rt, 2),
                'mz': round(mz, 4),
                'intensity': f"{intensity:.2e}",
                'charge': charge if charge != 0 else '-',
                'quality': round(quality, 3) if quality > 0 else '-',
                'rt_width': round(rt_width, 2) if rt_width > 0 else '-',
                'mz_width': round(mz_width, 4) if mz_width > 0 else '-',
            })

        return data

    def _extract_spectrum_data(self) -> List[Dict[str, Any]]:
        """Extract spectrum metadata for the spectrum browser table."""
        if self.exp is None:
            return []

        data = []
        for idx in range(self.exp.size()):
            spec = self.exp[idx]
            rt = spec.getRT()
            ms_level = spec.getMSLevel()
            n_peaks = spec.size()

            # Get peaks for TIC calculation
            mz_array, int_array = spec.get_peaks()
            tic = float(np.sum(int_array)) if len(int_array) > 0 else 0

            # Get m/z range
            mz_min = float(mz_array.min()) if len(mz_array) > 0 else 0
            mz_max = float(mz_array.max()) if len(mz_array) > 0 else 0

            # Get precursor info for MS2+
            precursor_mz = '-'
            precursor_charge = '-'
            if ms_level > 1:
                precursors = spec.getPrecursors()
                if precursors:
                    precursor_mz = round(precursors[0].getMZ(), 4)
                    charge = precursors[0].getCharge()
                    precursor_charge = charge if charge > 0 else '-'

            data.append({
                'idx': idx,
                'rt': round(rt, 2),
                'ms_level': ms_level,
                'n_peaks': n_peaks,
                'tic': f"{tic:.2e}",
                'mz_range': f"{mz_min:.1f}-{mz_max:.1f}" if n_peaks > 0 else '-',
                'precursor_mz': precursor_mz,
                'precursor_z': precursor_charge,
            })

        return data

    def show_spectrum_in_browser(self, spectrum_idx: int):
        """Display a spectrum in the 1D browser view."""
        if self.exp is None or spectrum_idx < 0 or spectrum_idx >= self.exp.size():
            return

        self.selected_spectrum_idx = spectrum_idx
        spec = self.exp[spectrum_idx]

        mz_array, int_array = spec.get_peaks()
        rt = spec.getRT()
        ms_level = spec.getMSLevel()

        if len(mz_array) == 0:
            ui.notify("Spectrum is empty", type="warning")
            return

        # Normalize intensities
        max_int = int_array.max() if len(int_array) > 0 else 1
        int_norm = (int_array / max_int) * 100

        # Create figure
        fig = go.Figure()

        # Color based on MS level
        color = '#00d4ff' if ms_level == 1 else '#ff6b6b'

        # Add spectrum as bars
        fig.add_trace(go.Bar(
            x=mz_array,
            y=int_norm,
            marker_color=color,
            width=0.5,
            opacity=0.8,
            hovertemplate='m/z: %{x:.4f}<br>Intensity: %{y:.1f}%<extra></extra>'
        ))

        # Title with spectrum info
        title = f"Spectrum #{spectrum_idx} | MS{ms_level} | RT={rt:.2f}s | {len(mz_array):,} peaks"

        # Add precursor line for MS2+
        if ms_level > 1:
            precursors = spec.getPrecursors()
            if precursors:
                prec_mz = precursors[0].getMZ()
                fig.add_vline(x=prec_mz, line_dash="dash", line_color="orange",
                              annotation_text=f"Precursor ({prec_mz:.2f})")
                title += f" | Precursor: {prec_mz:.4f}"

        fig.update_layout(
            title=dict(text=title, font=dict(size=14)),
            xaxis_title="m/z",
            yaxis_title="Relative Intensity (%)",
            template="plotly_dark",
            height=350,
            margin=dict(l=60, r=20, t=50, b=50),
            showlegend=False
        )

        fig.update_yaxes(range=[0, 105])

        # Update plot
        if self.spectrum_browser_plot is not None:
            self.spectrum_browser_plot.update_figure(fig)

        # Update navigation label
        if self.spectrum_nav_label is not None:
            self.spectrum_nav_label.set_text(f"Spectrum {spectrum_idx + 1} of {self.exp.size()}")

        # Update info label
        if self.spectrum_browser_info is not None:
            tic = float(np.sum(int_array))
            mz_range = f"{mz_array.min():.2f} - {mz_array.max():.2f}" if len(mz_array) > 0 else "-"
            self.spectrum_browser_info.set_text(
                f"RT: {rt:.2f}s | MS Level: {ms_level} | Peaks: {len(mz_array):,} | TIC: {tic:.2e} | m/z: {mz_range}"
            )

        # Update peak map to show the spectrum marker
        if self.show_spectrum_marker and self.df is not None:
            self.update_plot()

    def navigate_spectrum(self, direction: int):
        """Navigate to prev/next spectrum."""
        if self.exp is None or self.exp.size() == 0:
            return

        if self.selected_spectrum_idx is None:
            new_idx = 0
        else:
            new_idx = self.selected_spectrum_idx + direction

        # Clamp to valid range
        new_idx = max(0, min(self.exp.size() - 1, new_idx))
        self.show_spectrum_in_browser(new_idx)

    def navigate_spectrum_by_ms_level(self, direction: int, ms_level: int):
        """Navigate to prev/next spectrum of specific MS level."""
        if self.exp is None or self.exp.size() == 0:
            return

        start_idx = self.selected_spectrum_idx if self.selected_spectrum_idx is not None else 0

        if direction > 0:
            # Search forward
            for i in range(start_idx + 1, self.exp.size()):
                if self.exp[i].getMSLevel() == ms_level:
                    self.show_spectrum_in_browser(i)
                    return
        else:
            # Search backward
            for i in range(start_idx - 1, -1, -1):
                if self.exp[i].getMSLevel() == ms_level:
                    self.show_spectrum_in_browser(i)
                    return

        ui.notify(f"No more MS{ms_level} spectra in that direction", type="info")

    def load_featuremap(self, filepath: str) -> bool:
        """Load featureXML file."""
        try:
            if self.status_label:
                self.status_label.set_text(f"Loading features from {Path(filepath).name}...")
            ui.notify(f"Loading {filepath}...", type="info")

            self.feature_map = FeatureMap()
            FeatureXMLFile().load(filepath, self.feature_map)

            self.features_file = filepath
            self.selected_feature_idx = None
            self.feature_data = self._extract_feature_data()

            n_features = self.feature_map.size()
            if self.feature_info_label:
                self.feature_info_label.set_text(f"Features: {n_features:,}")
            if self.status_label:
                self.status_label.set_text("Ready")
            ui.notify(f"Loaded {n_features:,} features", type="positive")

            if self.feature_table is not None:
                self.feature_table.update_rows(self.feature_data)

            return True

        except Exception as e:
            if self.status_label:
                self.status_label.set_text(f"Error: {e}")
            ui.notify(f"Error loading features: {e}", type="negative")
            return False

    def clear_features(self):
        """Clear loaded feature map."""
        self.feature_map = None
        self.features_file = None
        self.feature_data = []
        self.selected_feature_idx = None
        if self.feature_info_label:
            self.feature_info_label.set_text("Features: None")
        if self.feature_table is not None:
            self.feature_table.update_rows([])
        ui.notify("Features cleared", type="info")

    def _extract_id_data(self) -> List[Dict[str, Any]]:
        """Extract peptide ID data for table display."""
        if not self.peptide_ids:
            return []

        data = []
        idx = 0
        for pep_id in self.peptide_ids:
            rt = pep_id.getRT()
            mz = pep_id.getMZ()
            hits = pep_id.getHits()

            if hits:
                best_hit = hits[0]
                sequence = best_hit.getSequence().toString()
                score = best_hit.getScore()
                charge = best_hit.getCharge()
            else:
                sequence = "-"
                score = 0
                charge = 0

            data.append({
                'idx': idx,
                'rt': round(rt, 2),
                'mz': round(mz, 4),
                'sequence': sequence[:30] + "..." if len(sequence) > 30 else sequence,
                'full_sequence': sequence,
                'charge': charge if charge != 0 else '-',
                'score': round(score, 4) if score != 0 else '-',
            })
            idx += 1

        return data

    def load_idxml(self, filepath: str) -> bool:
        """Load idXML file with peptide identifications."""
        try:
            if self.status_label:
                self.status_label.set_text(f"Loading IDs from {Path(filepath).name}...")
            ui.notify(f"Loading {filepath}...", type="info")

            self.protein_ids = []
            self.peptide_ids = []
            IdXMLFile().load(filepath, self.protein_ids, self.peptide_ids)

            self.id_file = filepath
            self.selected_id_idx = None
            self.id_data = self._extract_id_data()

            n_ids = len(self.peptide_ids)
            if self.id_info_label:
                self.id_info_label.set_text(f"IDs: {n_ids:,}")
            if self.status_label:
                self.status_label.set_text("Ready")
            ui.notify(f"Loaded {n_ids:,} peptide IDs", type="positive")

            if self.id_table is not None:
                self.id_table.update_rows(self.id_data)

            return True

        except Exception as e:
            if self.status_label:
                self.status_label.set_text(f"Error: {e}")
            ui.notify(f"Error loading IDs: {e}", type="negative")
            return False

    def clear_ids(self):
        """Clear loaded identifications."""
        self.peptide_ids = []
        self.protein_ids = []
        self.id_file = None
        self.id_data = []
        self.selected_id_idx = None
        if self.id_info_label:
            self.id_info_label.set_text("IDs: None")
        if self.id_table is not None:
            self.id_table.update_rows([])
        if self.spectrum_plot is not None:
            self.spectrum_plot.update_figure(go.Figure())
        if self.spectrum_info_label is not None:
            self.spectrum_info_label.set_text("Click an identification to view its annotated MS2 spectrum")
        ui.notify("Identifications cleared", type="info")

    def find_ms2_spectrum(self, rt: float, precursor_mz: float, rt_tolerance: float = 5.0, mz_tolerance: float = 0.5) -> Optional[MSSpectrum]:
        """Find the MS2 spectrum matching the given RT and precursor m/z."""
        if self.exp is None:
            return None

        best_spec = None
        best_rt_diff = float('inf')

        for spec in self.exp:
            if spec.getMSLevel() != 2:
                continue

            spec_rt = spec.getRT()
            if abs(spec_rt - rt) > rt_tolerance:
                continue

            # Check precursor m/z
            precursors = spec.getPrecursors()
            if precursors:
                prec_mz = precursors[0].getMZ()
                if abs(prec_mz - precursor_mz) <= mz_tolerance:
                    rt_diff = abs(spec_rt - rt)
                    if rt_diff < best_rt_diff:
                        best_rt_diff = rt_diff
                        best_spec = spec

        return best_spec

    def show_annotated_spectrum(self, id_idx: int):
        """Show annotated MS2 spectrum for the selected peptide ID."""
        if not self.peptide_ids or id_idx >= len(self.peptide_ids):
            return

        if self.exp is None:
            ui.notify("Load mzML file first to view spectra", type="warning")
            return

        pep_id = self.peptide_ids[id_idx]
        rt = pep_id.getRT()
        mz = pep_id.getMZ()

        hits = pep_id.getHits()
        if not hits:
            ui.notify("No peptide hits for this identification", type="warning")
            return

        best_hit = hits[0]
        sequence_str = best_hit.getSequence().toString()
        charge = best_hit.getCharge()

        # Find matching MS2 spectrum
        ms2_spec = self.find_ms2_spectrum(rt, mz)

        if ms2_spec is None:
            ui.notify(f"No MS2 spectrum found near RT={rt:.1f}s, m/z={mz:.2f}", type="warning")
            return

        # Get spectrum data
        mz_array, int_array = ms2_spec.get_peaks()

        if len(mz_array) == 0:
            ui.notify("MS2 spectrum is empty", type="warning")
            return

        # Create annotated spectrum plot
        fig = create_annotated_spectrum_plot(
            mz_array, int_array,
            sequence_str, charge, mz
        )

        # Update the plot
        if self.spectrum_plot is not None:
            self.spectrum_plot.update_figure(fig)

        if self.spectrum_info_label is not None:
            self.spectrum_info_label.set_text(
                f"Spectrum: {sequence_str} | RT: {rt:.2f}s | Precursor m/z: {mz:.4f} | Charge: {charge}+"
            )

    def create_tic_plot(self) -> go.Figure:
        """Create TIC (Total Ion Chromatogram) plot."""
        fig = go.Figure()

        if self.tic_rt is None or len(self.tic_rt) == 0:
            fig.update_layout(
                title="TIC - No data loaded",
                template="plotly_dark",
                height=200
            )
            return fig

        # Create TIC trace
        fig.add_trace(go.Scatter(
            x=self.tic_rt,
            y=self.tic_intensity,
            mode='lines',
            name='TIC',
            line=dict(color='#00d4ff', width=1),
            fill='tozeroy',
            fillcolor='rgba(0, 212, 255, 0.2)',
            hovertemplate='RT: %{x:.2f}s<br>Intensity: %{y:.2e}<extra></extra>'
        ))

        # Add view range indicator
        if self.view_rt_min is not None and self.view_rt_max is not None:
            fig.add_vrect(
                x0=self.view_rt_min,
                x1=self.view_rt_max,
                fillcolor="rgba(255, 255, 0, 0.15)",
                layer="below",
                line_width=1,
                line_color="rgba(255, 255, 0, 0.5)"
            )

        fig.update_layout(
            title=dict(text="Total Ion Chromatogram (TIC) - Click to view MS1 spectrum", font=dict(size=14)),
            xaxis_title="RT (s)",
            yaxis_title="Total Intensity",
            template="plotly_dark",
            height=200,
            margin=dict(l=60, r=20, t=40, b=40),
            showlegend=False,
            hovermode='x unified'
        )

        # Set x-axis range to match data
        if len(self.tic_rt) > 0:
            fig.update_xaxes(range=[self.rt_min, self.rt_max])

        return fig

    def update_tic_plot(self):
        """Update the TIC plot display."""
        # Skip if we're updating from a TIC interaction to prevent circular updates
        if self._updating_from_tic:
            return
        if self.tic_plot is not None:
            fig = self.create_tic_plot()
            self.tic_plot.update_figure(fig)

    def find_ms1_spectrum_at_rt(self, target_rt: float) -> Optional[MSSpectrum]:
        """Find the MS1 spectrum closest to the given RT."""
        if self.exp is None:
            return None

        best_spec = None
        best_rt_diff = float('inf')

        for spec in self.exp:
            if spec.getMSLevel() != 1:
                continue

            spec_rt = spec.getRT()
            rt_diff = abs(spec_rt - target_rt)
            if rt_diff < best_rt_diff:
                best_rt_diff = rt_diff
                best_spec = spec

        return best_spec

    def show_ms1_spectrum(self, rt: float):
        """Display MS1 spectrum at the given retention time."""
        if self.exp is None:
            ui.notify("Load mzML file first", type="warning")
            return

        spec = self.find_ms1_spectrum_at_rt(rt)
        if spec is None:
            ui.notify(f"No MS1 spectrum found near RT={rt:.1f}s", type="warning")
            return

        mz_array, int_array = spec.get_peaks()
        actual_rt = spec.getRT()

        if len(mz_array) == 0:
            ui.notify("Spectrum is empty", type="warning")
            return

        # Normalize intensities
        max_int = int_array.max() if len(int_array) > 0 else 1
        int_norm = (int_array / max_int) * 100

        # Create figure
        fig = go.Figure()

        # Add spectrum as bars
        fig.add_trace(go.Bar(
            x=mz_array,
            y=int_norm,
            marker_color='#00ff64',
            width=0.5,
            opacity=0.8,
            hovertemplate='m/z: %{x:.4f}<br>Intensity: %{y:.1f}%<extra></extra>'
        ))

        fig.update_layout(
            title=dict(
                text=f"MS1 Spectrum at RT={actual_rt:.2f}s ({len(mz_array):,} peaks)",
                font=dict(size=14)
            ),
            xaxis_title="m/z",
            yaxis_title="Relative Intensity (%)",
            template="plotly_dark",
            height=300,
            margin=dict(l=60, r=20, t=50, b=50),
            showlegend=False
        )

        fig.update_xaxes(range=[self.view_mz_min, self.view_mz_max] if self.view_mz_min else [0, 2000])
        fig.update_yaxes(range=[0, 105])

        # Update plot
        if self.ms1_spectrum_plot is not None:
            self.ms1_spectrum_plot.update_figure(fig)

        if self.ms1_spectrum_info_label is not None:
            tic_val = float(np.sum(int_array))
            self.ms1_spectrum_info_label.set_text(
                f"RT: {actual_rt:.2f}s | Peaks: {len(mz_array):,} | TIC: {tic_val:.2e}"
            )

    def zoom_to_feature(self, feature_idx: int, padding: float = 0.2):
        """Zoom to a specific feature."""
        if self.feature_map is None or feature_idx >= self.feature_map.size():
            return

        self.selected_feature_idx = feature_idx
        self.selected_id_idx = None
        feature = self.feature_map[feature_idx]

        rt = feature.getRT()
        mz = feature.getMZ()

        hulls = feature.getConvexHulls()
        if hulls:
            all_points = []
            for hull in hulls:
                points = hull.getHullPoints()
                all_points.extend([(p[0], p[1]) for p in points])

            if all_points:
                rt_coords = [p[0] for p in all_points]
                mz_coords = [p[1] for p in all_points]
                feat_rt_min, feat_rt_max = min(rt_coords), max(rt_coords)
                feat_mz_min, feat_mz_max = min(mz_coords), max(mz_coords)
            else:
                feat_rt_min, feat_rt_max = rt - 10, rt + 10
                feat_mz_min, feat_mz_max = mz - 2, mz + 2
        else:
            feat_rt_min, feat_rt_max = rt - 10, rt + 10
            feat_mz_min, feat_mz_max = mz - 2, mz + 2

        rt_range = max(feat_rt_max - feat_rt_min, 20)
        mz_range = max(feat_mz_max - feat_mz_min, 4)

        rt_pad = rt_range * padding
        mz_pad = mz_range * padding

        self.view_rt_min = max(self.rt_min, feat_rt_min - rt_pad)
        self.view_rt_max = min(self.rt_max, feat_rt_max + rt_pad)
        self.view_mz_min = max(self.mz_min, feat_mz_min - mz_pad)
        self.view_mz_max = min(self.mz_max, feat_mz_max + mz_pad)

        self.update_plot()
        ui.notify(f"Zoomed to feature {feature_idx + 1}", type="info")

    def zoom_to_id(self, id_idx: int, padding: float = 0.3):
        """Zoom to a specific peptide identification and show annotated spectrum."""
        if not self.peptide_ids or id_idx >= len(self.peptide_ids):
            return

        self.selected_id_idx = id_idx
        self.selected_feature_idx = None
        pep_id = self.peptide_ids[id_idx]

        rt = pep_id.getRT()
        mz = pep_id.getMZ()

        rt_window = 30
        mz_window = 5

        self.view_rt_min = max(self.rt_min, rt - rt_window)
        self.view_rt_max = min(self.rt_max, rt + rt_window)
        self.view_mz_min = max(self.mz_min, mz - mz_window)
        self.view_mz_max = min(self.mz_max, mz + mz_window)

        self.update_plot()
        self.show_annotated_spectrum(id_idx)
        ui.notify(f"Zoomed to ID {id_idx + 1}", type="info")

    def _data_to_plot_pixel(self, rt: float, mz: float) -> Tuple[int, int]:
        """Convert RT/m/z to pixel coordinates."""
        rt_range = self.view_rt_max - self.view_rt_min
        mz_range = self.view_mz_max - self.view_mz_min

        if rt_range == 0 or mz_range == 0:
            return (0, 0)

        x = int((rt - self.view_rt_min) / rt_range * self.plot_width)
        y = int((1 - (mz - self.view_mz_min) / mz_range) * self.plot_height)

        return (x, y)

    def _is_in_view(self, rt: float, mz: float) -> bool:
        """Check if point is in current view."""
        return (self.view_rt_min <= rt <= self.view_rt_max and
                self.view_mz_min <= mz <= self.view_mz_max)

    def _feature_intersects_view(self, rt_min: float, rt_max: float,
                                  mz_min: float, mz_max: float) -> bool:
        """Check if feature intersects current view."""
        return not (rt_max < self.view_rt_min or rt_min > self.view_rt_max or
                    mz_max < self.view_mz_min or mz_min > self.view_mz_max)

    def _draw_features_on_plot(self, img: Image.Image) -> Image.Image:
        """Draw feature overlays."""
        if self.feature_map is None or self.feature_map.size() == 0:
            return img

        img = img.convert('RGBA')
        overlay = Image.new('RGBA', img.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)

        features_drawn = 0
        max_features = 10000

        for idx, feature in enumerate(self.feature_map):
            if features_drawn >= max_features:
                break

            is_selected = (idx == self.selected_feature_idx)
            rt = feature.getRT()
            mz = feature.getMZ()

            hulls = feature.getConvexHulls()
            if hulls:
                all_points = []
                for hull in hulls:
                    points = hull.getHullPoints()
                    all_points.extend([(p[0], p[1]) for p in points])

                if all_points:
                    rt_coords = [p[0] for p in all_points]
                    mz_coords = [p[1] for p in all_points]
                    feat_rt_min, feat_rt_max = min(rt_coords), max(rt_coords)
                    feat_mz_min, feat_mz_max = min(mz_coords), max(mz_coords)
                else:
                    feat_rt_min, feat_rt_max = rt - 1, rt + 1
                    feat_mz_min, feat_mz_max = mz - 0.5, mz + 0.5
            else:
                feat_rt_min, feat_rt_max = rt - 1, rt + 1
                feat_mz_min, feat_mz_max = mz - 0.5, mz + 0.5

            if not self._feature_intersects_view(feat_rt_min, feat_rt_max,
                                                  feat_mz_min, feat_mz_max):
                continue

            features_drawn += 1

            hull_color = self.selected_color if is_selected else self.hull_color
            bbox_color = self.selected_color if is_selected else self.bbox_color
            centroid_color = self.selected_color if is_selected else self.centroid_color
            line_width = 3 if is_selected else 1

            if self.show_convex_hulls and hulls:
                for hull in hulls:
                    points = hull.getHullPoints()
                    if len(points) >= 3:
                        pixel_points = [self._data_to_plot_pixel(p[0], p[1]) for p in points]
                        pixel_points.append(pixel_points[0])
                        fill_alpha = 100 if is_selected else 50
                        draw.polygon(pixel_points, outline=hull_color,
                                    fill=(*hull_color[:3], fill_alpha))

            if self.show_bounding_boxes:
                top_left = self._data_to_plot_pixel(feat_rt_min, feat_mz_max)
                bottom_right = self._data_to_plot_pixel(feat_rt_max, feat_mz_min)
                draw.rectangle([top_left, bottom_right], outline=bbox_color, width=line_width)

            if self.show_centroids:
                cx, cy = self._data_to_plot_pixel(rt, mz)
                r = 5 if is_selected else 3
                draw.ellipse([cx-r, cy-r, cx+r, cy+r], fill=centroid_color,
                            outline=(255, 255, 255, 255))

        img = Image.alpha_composite(img, overlay)
        return img

    def _draw_spectrum_marker_on_plot(self, img: Image.Image) -> Image.Image:
        """Draw a horizontal line at the selected spectrum's RT."""
        if self.selected_spectrum_idx is None or self.exp is None:
            return img

        spec = self.exp[self.selected_spectrum_idx]
        rt = spec.getRT()

        # Check if RT is in view
        if rt < self.view_rt_min or rt > self.view_rt_max:
            return img

        img = img.convert('RGBA')
        overlay = Image.new('RGBA', img.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)

        # Calculate x position for the RT
        x, _ = self._data_to_plot_pixel(rt, self.view_mz_min)

        # Draw vertical line across the full height
        ms_level = spec.getMSLevel()
        line_color = (0, 212, 255, 200) if ms_level == 1 else (255, 107, 107, 200)  # cyan for MS1, red for MS2

        draw.line([(x, 0), (x, self.plot_height)], fill=line_color, width=2)

        # Draw small label at top
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 10)
        except:
            font = ImageFont.load_default()

        label = f"MS{ms_level} #{self.selected_spectrum_idx}"
        draw.text((x + 4, 4), label, fill=line_color, font=font)

        img = Image.alpha_composite(img, overlay)
        return img

    def _draw_ids_on_plot(self, img: Image.Image) -> Image.Image:
        """Draw peptide ID precursor positions."""
        if not self.peptide_ids or not self.show_ids:
            return img

        img = img.convert('RGBA')
        overlay = Image.new('RGBA', img.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)

        for idx, pep_id in enumerate(self.peptide_ids):
            rt = pep_id.getRT()
            mz = pep_id.getMZ()

            if not self._is_in_view(rt, mz):
                continue

            is_selected = (idx == self.selected_id_idx)
            color = self.id_selected_color if is_selected else self.id_color

            cx, cy = self._data_to_plot_pixel(rt, mz)

            r = 6 if is_selected else 4
            diamond = [(cx, cy - r), (cx + r, cy), (cx, cy + r), (cx - r, cy)]
            draw.polygon(diamond, fill=color, outline=(255, 255, 255, 255))

            if is_selected:
                draw.line([(cx - r - 3, cy), (cx + r + 3, cy)], fill=color, width=2)
                draw.line([(cx, cy - r - 3), (cx, cy + r + 3)], fill=color, width=2)

        img = Image.alpha_composite(img, overlay)
        return img

    def _draw_axes(self, canvas: Image.Image) -> Image.Image:
        """Draw axes on canvas."""
        draw = ImageDraw.Draw(canvas)

        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12)
            title_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 14)
        except:
            try:
                font = ImageFont.truetype("/usr/share/fonts/TTF/DejaVuSans.ttf", 12)
                title_font = ImageFont.truetype("/usr/share/fonts/TTF/DejaVuSans.ttf", 14)
            except:
                font = ImageFont.load_default()
                title_font = font

        plot_left = self.margin_left
        plot_right = self.margin_left + self.plot_width
        plot_top = self.margin_top
        plot_bottom = self.margin_top + self.plot_height

        draw.rectangle([plot_left, plot_top, plot_right, plot_bottom],
                       outline=self.axis_color, width=1)

        # X-axis
        rt_ticks = calculate_nice_ticks(self.view_rt_min, self.view_rt_max, num_ticks=8)
        rt_range = self.view_rt_max - self.view_rt_min

        for tick_val in rt_ticks:
            if self.view_rt_min <= tick_val <= self.view_rt_max:
                x_frac = (tick_val - self.view_rt_min) / rt_range
                x = plot_left + int(x_frac * self.plot_width)

                draw.line([(x, plot_bottom), (x, plot_bottom + 5)], fill=self.tick_color, width=1)
                draw.line([(x, plot_top), (x, plot_bottom)], fill=self.grid_color, width=1)

                label = format_tick_label(tick_val, rt_range)
                bbox = draw.textbbox((0, 0), label, font=font)
                label_width = bbox[2] - bbox[0]
                draw.text((x - label_width // 2, plot_bottom + 8), label, fill=self.label_color, font=font)

        x_title = "RT (s)"
        bbox = draw.textbbox((0, 0), x_title, font=title_font)
        title_width = bbox[2] - bbox[0]
        draw.text((plot_left + self.plot_width // 2 - title_width // 2, plot_bottom + 28),
                  x_title, fill=self.label_color, font=title_font)

        # Y-axis
        mz_ticks = calculate_nice_ticks(self.view_mz_min, self.view_mz_max, num_ticks=8)
        mz_range = self.view_mz_max - self.view_mz_min

        for tick_val in mz_ticks:
            if self.view_mz_min <= tick_val <= self.view_mz_max:
                y_frac = 1 - (tick_val - self.view_mz_min) / mz_range
                y = plot_top + int(y_frac * self.plot_height)

                draw.line([(plot_left - 5, y), (plot_left, y)], fill=self.tick_color, width=1)
                draw.line([(plot_left, y), (plot_right, y)], fill=self.grid_color, width=1)

                label = format_tick_label(tick_val, mz_range)
                bbox = draw.textbbox((0, 0), label, font=font)
                label_width = bbox[2] - bbox[0]
                label_height = bbox[3] - bbox[1]
                draw.text((plot_left - label_width - 10, y - label_height // 2), label, fill=self.label_color, font=font)

        # Y-axis title
        y_title = "m/z"
        txt_img = Image.new('RGBA', (100, 30), (0, 0, 0, 0))
        txt_draw = ImageDraw.Draw(txt_img)
        txt_draw.text((0, 0), y_title, fill=self.label_color, font=title_font)
        txt_img = txt_img.rotate(90, expand=True)

        y_title_x = 5
        y_title_y = plot_top + self.plot_height // 2 - txt_img.height // 2
        canvas.paste(txt_img, (y_title_x, y_title_y), txt_img)

        return canvas

    def render_image(self) -> str:
        """Render current view using datashader."""
        if self.df is None or len(self.df) == 0:
            return ""

        mask = (
            (self.df['rt'] >= self.view_rt_min) &
            (self.df['rt'] <= self.view_rt_max) &
            (self.df['mz'] >= self.view_mz_min) &
            (self.df['mz'] <= self.view_mz_max)
        )
        view_df = self.df[mask]

        if len(view_df) == 0:
            return ""

        ds_canvas = ds.Canvas(
            plot_width=self.plot_width,
            plot_height=self.plot_height,
            x_range=(self.view_rt_min, self.view_rt_max),
            y_range=(self.view_mz_min, self.view_mz_max)
        )

        agg = ds_canvas.points(view_df, 'rt', 'mz', ds.mean('log_intensity'))
        img = tf.shade(agg, cmap=fire, how='linear')
        img = tf.set_background(img, 'black')

        plot_img = img.to_pil()

        if self.feature_map is not None:
            plot_img = self._draw_features_on_plot(plot_img)

        if self.peptide_ids:
            plot_img = self._draw_ids_on_plot(plot_img)

        if self.show_spectrum_marker:
            plot_img = self._draw_spectrum_marker_on_plot(plot_img)

        canvas = Image.new('RGBA', (self.canvas_width, self.canvas_height), (20, 20, 25, 255))
        plot_img_rgba = plot_img.convert('RGBA')
        canvas.paste(plot_img_rgba, (self.margin_left, self.margin_top))

        canvas = self._draw_axes(canvas)

        buffer = io.BytesIO()
        canvas.save(buffer, format='PNG')
        buffer.seek(0)

        return base64.b64encode(buffer.getvalue()).decode('utf-8')

    def render_faims_image(self, cv: float) -> str:
        """Render a single FAIMS CV peak map using datashader."""
        if cv not in self.faims_data or len(self.faims_data[cv]) == 0:
            return ""

        cv_df = self.faims_data[cv]

        mask = (
            (cv_df['rt'] >= self.view_rt_min) &
            (cv_df['rt'] <= self.view_rt_max) &
            (cv_df['mz'] >= self.view_mz_min) &
            (cv_df['mz'] <= self.view_mz_max)
        )
        view_df = cv_df[mask]

        if len(view_df) == 0:
            return ""

        # Smaller plot size for FAIMS panels
        faims_plot_width = self.plot_width // max(1, min(len(self.faims_cvs), 4))
        faims_plot_height = self.plot_height

        ds_canvas = ds.Canvas(
            plot_width=faims_plot_width,
            plot_height=faims_plot_height,
            x_range=(self.view_rt_min, self.view_rt_max),
            y_range=(self.view_mz_min, self.view_mz_max)
        )

        agg = ds_canvas.points(view_df, 'rt', 'mz', ds.mean('log_intensity'))
        img = tf.shade(agg, cmap=fire, how='linear')
        img = tf.set_background(img, 'black')

        plot_img = img.to_pil()

        # Add CV label at top
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 14)
        except:
            font = ImageFont.load_default()

        draw = ImageDraw.Draw(plot_img)
        label = f"CV: {cv:.1f}V"
        bbox = draw.textbbox((0, 0), label, font=font)
        label_width = bbox[2] - bbox[0]
        draw.rectangle([(5, 5), (label_width + 15, 25)], fill=(0, 0, 0, 180))
        draw.text((10, 7), label, fill=(255, 255, 255, 255), font=font)

        # Add border
        draw.rectangle([(0, 0), (faims_plot_width - 1, faims_plot_height - 1)],
                       outline=(100, 100, 100, 255), width=1)

        buffer = io.BytesIO()
        plot_img.save(buffer, format='PNG')
        buffer.seek(0)

        return base64.b64encode(buffer.getvalue()).decode('utf-8')

    def update_faims_plots(self):
        """Update all FAIMS CV peak map panels."""
        if not self.has_faims or not self.show_faims_view:
            return

        for cv in self.faims_cvs:
            if cv in self.faims_images and self.faims_images[cv] is not None:
                img_data = self.render_faims_image(cv)
                if img_data:
                    self.faims_images[cv].set_source(f"data:image/png;base64,{img_data}")

    def update_plot(self):
        """Update displayed plot."""
        if self.df is None:
            return

        if self.status_label:
            self.status_label.set_text("Rendering...")

        img_data = self.render_image()
        if img_data and self.image_element:
            self.image_element.set_source(f"data:image/png;base64,{img_data}")

        # Update FAIMS plots if enabled
        if self.has_faims and self.show_faims_view:
            self.update_faims_plots()

        if self.rt_range_label:
            self.rt_range_label.set_text(f"RT: {self.view_rt_min:.2f} - {self.view_rt_max:.2f} s")
        if self.mz_range_label:
            self.mz_range_label.set_text(f"m/z: {self.view_mz_min:.2f} - {self.view_mz_max:.2f}")

        # Update TIC plot (shows current view range)
        self.update_tic_plot()

        if self.status_label:
            self.status_label.set_text("Ready")

    def reset_view(self):
        """Reset to full view."""
        if self.df is None:
            return
        self.view_rt_min = self.rt_min
        self.view_rt_max = self.rt_max
        self.view_mz_min = self.mz_min
        self.view_mz_max = self.mz_max
        self.selected_feature_idx = None
        self.selected_id_idx = None
        self.update_plot()

    def zoom_in(self, factor=0.5):
        """Zoom in."""
        if self.df is None:
            return
        rt_center = (self.view_rt_min + self.view_rt_max) / 2
        mz_center = (self.view_mz_min + self.view_mz_max) / 2
        rt_range = (self.view_rt_max - self.view_rt_min) * factor / 2
        mz_range = (self.view_mz_max - self.view_mz_min) * factor / 2

        self.view_rt_min = rt_center - rt_range
        self.view_rt_max = rt_center + rt_range
        self.view_mz_min = mz_center - mz_range
        self.view_mz_max = mz_center + mz_range
        self.update_plot()

    def zoom_out(self, factor=2.0):
        """Zoom out."""
        if self.df is None:
            return
        rt_center = (self.view_rt_min + self.view_rt_max) / 2
        mz_center = (self.view_mz_min + self.view_mz_max) / 2
        rt_range = (self.view_rt_max - self.view_rt_min) * factor / 2
        mz_range = (self.view_mz_max - self.view_mz_min) * factor / 2

        self.view_rt_min = max(self.rt_min, rt_center - rt_range)
        self.view_rt_max = min(self.rt_max, rt_center + rt_range)
        self.view_mz_min = max(self.mz_min, mz_center - mz_range)
        self.view_mz_max = min(self.mz_max, mz_center + mz_range)
        self.update_plot()

    def pan(self, rt_frac=0, mz_frac=0):
        """Pan view."""
        if self.df is None:
            return
        rt_shift = (self.view_rt_max - self.view_rt_min) * rt_frac
        mz_shift = (self.view_mz_max - self.view_mz_min) * mz_frac

        if self.view_rt_min + rt_shift < self.rt_min:
            rt_shift = self.rt_min - self.view_rt_min
        if self.view_rt_max + rt_shift > self.rt_max:
            rt_shift = self.rt_max - self.view_rt_max
        if self.view_mz_min + mz_shift < self.mz_min:
            mz_shift = self.mz_min - self.view_mz_min
        if self.view_mz_max + mz_shift > self.mz_max:
            mz_shift = self.mz_max - self.view_mz_max

        self.view_rt_min += rt_shift
        self.view_rt_max += rt_shift
        self.view_mz_min += mz_shift
        self.view_mz_max += mz_shift
        self.update_plot()

    def apply_custom_range(self, rt_min, rt_max, mz_min, mz_max):
        """Apply custom range."""
        if self.df is None:
            return
        try:
            self.view_rt_min = max(self.rt_min, float(rt_min))
            self.view_rt_max = min(self.rt_max, float(rt_max))
            self.view_mz_min = max(self.mz_min, float(mz_min))
            self.view_mz_max = min(self.mz_max, float(mz_max))
            self.update_plot()
        except ValueError:
            ui.notify("Invalid range values", type="warning")


def create_ui():
    """Create NiceGUI interface."""
    global _viewer_instance

    viewer = MzMLViewer()
    _viewer_instance = viewer

    ui.dark_mode().enable()

    with ui.column().classes('w-full items-center p-4'):
        ui.label('mzML Peak Map Viewer').classes('text-3xl font-bold mb-2')
        ui.label('High-performance visualization with Datashader + pyOpenMS').classes('text-gray-400 mb-4')

        # File upload section
        with ui.card().classes('w-full max-w-6xl mb-4'):
            ui.label('Load Data').classes('text-xl font-semibold mb-2')

            with ui.row().classes('w-full items-end gap-4 flex-wrap'):
                # mzML
                with ui.column().classes('flex-1 min-w-64'):
                    ui.label('mzML File (Peak Data)').classes('text-sm text-gray-400')
                    with ui.row().classes('w-full items-end gap-2'):
                        async def handle_mzml_upload(e):
                            content = e.content.read()
                            temp_path = Path('/tmp') / e.name
                            temp_path.write_bytes(content)
                            if viewer.load_mzml(str(temp_path)):
                                viewer.update_plot()

                        ui.upload(label='Upload mzML', on_upload=handle_mzml_upload,
                                  auto_upload=True).props('accept=.mzML,.mzml').classes('w-40')

                        mzml_input = ui.input(placeholder='/path/to/file.mzML').classes('flex-1')

                        async def load_mzml_path():
                            path = mzml_input.value
                            if path and Path(path).exists():
                                if viewer.load_mzml(path):
                                    viewer.update_plot()
                            else:
                                ui.notify("File not found", type="warning")

                        ui.button('Load', on_click=load_mzml_path).props('color=primary dense')

                # FeatureXML
                with ui.column().classes('flex-1 min-w-64'):
                    ui.label('FeatureXML (Features)').classes('text-sm text-gray-400')
                    with ui.row().classes('w-full items-end gap-2'):
                        async def handle_feature_upload(e):
                            content = e.content.read()
                            temp_path = Path('/tmp') / e.name
                            temp_path.write_bytes(content)
                            if viewer.load_featuremap(str(temp_path)):
                                viewer.update_plot()

                        ui.upload(label='Upload featureXML', on_upload=handle_feature_upload,
                                  auto_upload=True).props('accept=.featureXML,.xml').classes('w-40')

                        feature_input = ui.input(placeholder='/path/to/features.featureXML').classes('flex-1')

                        async def load_feature_path():
                            path = feature_input.value
                            if path and Path(path).exists():
                                if viewer.load_featuremap(path):
                                    viewer.update_plot()
                            else:
                                ui.notify("File not found", type="warning")

                        ui.button('Load', on_click=load_feature_path).props('color=primary dense')

                        def clear_features():
                            viewer.clear_features()
                            viewer.update_plot()

                        ui.button('Clear', on_click=clear_features).props('color=negative dense')

                # idXML
                with ui.column().classes('flex-1 min-w-64'):
                    ui.label('idXML (Identifications)').classes('text-sm text-gray-400')
                    with ui.row().classes('w-full items-end gap-2'):
                        async def handle_id_upload(e):
                            content = e.content.read()
                            temp_path = Path('/tmp') / e.name
                            temp_path.write_bytes(content)
                            if viewer.load_idxml(str(temp_path)):
                                viewer.update_plot()

                        ui.upload(label='Upload idXML', on_upload=handle_id_upload,
                                  auto_upload=True).props('accept=.idXML,.xml').classes('w-40')

                        id_input = ui.input(placeholder='/path/to/ids.idXML').classes('flex-1')

                        async def load_id_path():
                            path = id_input.value
                            if path and Path(path).exists():
                                if viewer.load_idxml(path):
                                    viewer.update_plot()
                            else:
                                ui.notify("File not found", type="warning")

                        ui.button('Load', on_click=load_id_path).props('color=primary dense')

                        def clear_ids():
                            viewer.clear_ids()
                            viewer.update_plot()

                        ui.button('Clear', on_click=clear_ids).props('color=negative dense')

        # Info bar
        with ui.row().classes('w-full justify-center gap-6 mb-2 flex-wrap'):
            viewer.info_label = ui.label('No file loaded').classes('text-gray-400')
            viewer.feature_info_label = ui.label('Features: None').classes('text-cyan-400')
            viewer.id_info_label = ui.label('IDs: None').classes('text-orange-400')
            viewer.faims_info_label = ui.label('').classes('text-purple-400')
            viewer.faims_info_label.set_visibility(False)
            viewer.status_label = ui.label('Ready').classes('text-green-400')

        # Range display
        with ui.row().classes('w-full justify-center gap-8 mb-2'):
            viewer.rt_range_label = ui.label('RT: -- - -- s').classes('text-blue-300')
            viewer.mz_range_label = ui.label('m/z: -- - --').classes('text-blue-300')

        # Display options
        with ui.row().classes('w-full justify-center gap-4 mb-2 flex-wrap'):
            ui.label('Show:').classes('text-gray-400')

            def toggle_centroids():
                viewer.show_centroids = centroid_cb.value
                if viewer.df is not None:
                    viewer.update_plot()

            centroid_cb = ui.checkbox('Centroids', value=True, on_change=toggle_centroids).classes('text-green-400')

            def toggle_bboxes():
                viewer.show_bounding_boxes = bbox_cb.value
                if viewer.df is not None:
                    viewer.update_plot()

            bbox_cb = ui.checkbox('Bounding Boxes', value=True, on_change=toggle_bboxes).classes('text-yellow-400')

            def toggle_hulls():
                viewer.show_convex_hulls = hull_cb.value
                if viewer.df is not None:
                    viewer.update_plot()

            hull_cb = ui.checkbox('Convex Hulls', value=True, on_change=toggle_hulls).classes('text-cyan-400')

            def toggle_ids():
                viewer.show_ids = ids_cb.value
                if viewer.df is not None:
                    viewer.update_plot()

            ids_cb = ui.checkbox('Identifications', value=True, on_change=toggle_ids).classes('text-orange-400')

            def toggle_spectrum_marker():
                viewer.show_spectrum_marker = spectrum_marker_cb.value
                if viewer.df is not None:
                    viewer.update_plot()

            spectrum_marker_cb = ui.checkbox('Spectrum Marker', value=True, on_change=toggle_spectrum_marker).classes('text-pink-400')

            # FAIMS toggle (hidden by default, shown when FAIMS data is detected)
            def toggle_faims_view():
                viewer.show_faims_view = faims_toggle.value
                if viewer.faims_container:
                    viewer.faims_container.set_visibility(viewer.show_faims_view)
                if viewer.df is not None and viewer.show_faims_view:
                    viewer.update_faims_plots()

            faims_toggle = ui.checkbox('FAIMS Multi-CV View', value=False, on_change=toggle_faims_view).classes('text-purple-400')
            faims_toggle.set_visibility(False)
            viewer.faims_toggle = faims_toggle

        # TIC Plot (clickable to show MS1 spectrum, zoomable to update peak map)
        with ui.card().classes('w-full max-w-6xl'):
            ui.label('TIC - Click to view spectrum, drag to zoom RT range').classes('text-xs text-gray-500 mb-1')
            viewer.tic_plot = ui.plotly(viewer.create_tic_plot()).classes('w-full')

            # Handle click on TIC plot - show spectrum and center peak map
            def on_tic_click(e):
                try:
                    if e.args and 'points' in e.args and e.args['points']:
                        point = e.args['points'][0]
                        if 'x' in point:
                            rt = point['x']
                            viewer.show_ms1_spectrum(rt)
                            # Also center the peak map on this RT
                            rt_range = viewer.view_rt_max - viewer.view_rt_min
                            viewer.view_rt_min = max(viewer.rt_min, rt - rt_range / 2)
                            viewer.view_rt_max = min(viewer.rt_max, rt + rt_range / 2)
                            viewer.update_plot()
                except Exception:
                    pass

            viewer.tic_plot.on('plotly_click', on_tic_click)

            # Handle zoom/pan on TIC plot - sync RT range to peak map
            def on_tic_relayout(e):
                try:
                    if e.args:
                        args = e.args
                        # Check for x-axis range changes (zoom or pan)
                        if 'xaxis.range[0]' in args and 'xaxis.range[1]' in args:
                            new_rt_min = float(args['xaxis.range[0]'])
                            new_rt_max = float(args['xaxis.range[1]'])
                            # Clamp to data bounds
                            viewer.view_rt_min = max(viewer.rt_min, new_rt_min)
                            viewer.view_rt_max = min(viewer.rt_max, new_rt_max)
                            # Set flag to prevent TIC reset during update
                            viewer._updating_from_tic = True
                            viewer.update_plot()
                            viewer._updating_from_tic = False
                        elif 'xaxis.autorange' in args and args['xaxis.autorange']:
                            # Reset to full range
                            viewer.view_rt_min = viewer.rt_min
                            viewer.view_rt_max = viewer.rt_max
                            viewer._updating_from_tic = True
                            viewer.update_plot()
                            viewer._updating_from_tic = False
                except Exception:
                    viewer._updating_from_tic = False

            viewer.tic_plot.on('plotly_relayout', on_tic_relayout)

        # Main visualization area - peak map with spectrum browser overlay
        with ui.card().classes('w-full max-w-6xl p-2'):
            # Peak map
            with ui.row().classes('w-full items-start gap-0'):
                # Peak map image
                with ui.column().classes('flex-none'):
                    viewer.image_element = ui.image().classes('w-full').style(
                        f'width: {viewer.canvas_width}px; height: {viewer.canvas_height}px; background: #141419;'
                    )

            # 1D Spectrum Browser Plot (directly below peak map, same width)
            with ui.column().classes('w-full mt-2'):
                # Navigation and info row
                with ui.row().classes('w-full items-center gap-2 mb-1'):
                    ui.label('1D Spectrum:').classes('text-sm font-semibold text-gray-300')
                    ui.button('|<', on_click=lambda: viewer.show_spectrum_in_browser(0)).props('dense size=sm').tooltip('First')
                    ui.button('< MS1', on_click=lambda: viewer.navigate_spectrum_by_ms_level(-1, 1)).props('dense size=sm color=cyan').tooltip('Prev MS1')
                    ui.button('<', on_click=lambda: viewer.navigate_spectrum(-1)).props('dense size=sm').tooltip('Prev')

                    viewer.spectrum_nav_label = ui.label('No spectrum').classes('mx-2 text-gray-400 text-sm')

                    ui.button('>', on_click=lambda: viewer.navigate_spectrum(1)).props('dense size=sm').tooltip('Next')
                    ui.button('MS1 >', on_click=lambda: viewer.navigate_spectrum_by_ms_level(1, 1)).props('dense size=sm color=cyan').tooltip('Next MS1')
                    ui.button('>|', on_click=lambda: viewer.show_spectrum_in_browser(viewer.exp.size() - 1 if viewer.exp else 0)).props('dense size=sm').tooltip('Last')

                    ui.label('|').classes('mx-1 text-gray-600')
                    ui.button('< MS2', on_click=lambda: viewer.navigate_spectrum_by_ms_level(-1, 2)).props('dense size=sm color=orange').tooltip('Prev MS2')
                    ui.button('MS2 >', on_click=lambda: viewer.navigate_spectrum_by_ms_level(1, 2)).props('dense size=sm color=orange').tooltip('Next MS2')

                    ui.element('div').classes('flex-grow')  # Spacer

                    viewer.spectrum_browser_info = ui.label('Click TIC or use spectrum table to select').classes('text-xs text-gray-500')

                # Spectrum plot
                viewer.spectrum_browser_plot = ui.plotly(go.Figure()).classes('w-full').style(f'max-width: {viewer.canvas_width}px;')

        # FAIMS Multi-CV Peak Maps (hidden by default)
        faims_container = ui.card().classes('w-full max-w-6xl mt-2 p-2')
        faims_container.set_visibility(False)
        viewer.faims_container = faims_container

        with faims_container:
            ui.label('FAIMS Compensation Voltage Peak Maps').classes('text-lg font-semibold mb-2 text-purple-300')
            ui.label('Separate peak maps for each CV value - zoom/pan is synchronized').classes('text-xs text-gray-500 mb-2')

            # Container for dynamic FAIMS images
            faims_row = ui.row().classes('w-full gap-1 flex-wrap justify-center')

            # Note: Actual images will be created dynamically when FAIMS data is loaded
            # We need to create a method to dynamically populate this container

            def create_faims_images():
                """Create FAIMS image elements dynamically based on detected CVs."""
                faims_row.clear()
                viewer.faims_images = {}

                if not viewer.has_faims:
                    return

                n_cvs = len(viewer.faims_cvs)
                # Calculate width for each panel (max 4 per row)
                panel_width = viewer.plot_width // max(1, min(n_cvs, 4))
                panel_height = viewer.plot_height

                with faims_row:
                    for cv in viewer.faims_cvs:
                        with ui.column().classes('flex-none'):
                            img = ui.image().style(
                                f'width: {panel_width}px; height: {panel_height}px; background: #141419;'
                            )
                            viewer.faims_images[cv] = img

            # Store the function reference for later use
            viewer._create_faims_images = create_faims_images

        # Navigation controls
        with ui.row().classes('justify-center gap-2 mt-2'):
            ui.button('Reset View', on_click=viewer.reset_view).props('color=secondary')
            ui.button('Zoom In', on_click=lambda: viewer.zoom_in(0.5)).props('color=primary')
            ui.button('Zoom Out', on_click=lambda: viewer.zoom_out(2.0)).props('color=primary')
            ui.button('← Pan Left', on_click=lambda: viewer.pan(rt_frac=-0.25)).props('color=accent')
            ui.button('→ Pan Right', on_click=lambda: viewer.pan(rt_frac=0.25)).props('color=accent')
            ui.button('↑ Pan Up', on_click=lambda: viewer.pan(mz_frac=0.25)).props('color=accent')
            ui.button('↓ Pan Down', on_click=lambda: viewer.pan(mz_frac=-0.25)).props('color=accent')

        # MS1 Spectrum Viewer (from TIC click) - now in expansion
        with ui.expansion('TIC Spectrum Viewer', icon='show_chart').classes('w-full max-w-6xl mt-2'):
            viewer.ms1_spectrum_info_label = ui.label(
                'Click on the TIC plot above to display an MS1 spectrum'
            ).classes('text-sm text-gray-400 mb-2')
            viewer.ms1_spectrum_plot = ui.plotly(go.Figure()).classes('w-full')

        # Annotated MS2 Spectrum Viewer
        with ui.card().classes('w-full max-w-6xl mt-4'):
            ui.label('Annotated MS2 Spectrum').classes('text-xl font-semibold mb-2')
            viewer.spectrum_info_label = ui.label(
                'Click an identification to view its annotated MS2 spectrum'
            ).classes('text-sm text-gray-400 mb-2')
            viewer.spectrum_plot = ui.plotly(go.Figure()).classes('w-full')

        # Tables section
        with ui.row().classes('w-full max-w-6xl mt-4 gap-4 flex-wrap'):
            # Feature Table
            with ui.card().classes('flex-1 min-w-96'):
                ui.label('Features').classes('text-xl font-semibold mb-2')
                ui.label('Click a row to zoom to that feature').classes('text-sm text-gray-400 mb-2')

                feature_columns = [
                    {'name': 'idx', 'label': '#', 'field': 'idx', 'sortable': True, 'align': 'left'},
                    {'name': 'rt', 'label': 'RT (s)', 'field': 'rt', 'sortable': True, 'align': 'right'},
                    {'name': 'mz', 'label': 'm/z', 'field': 'mz', 'sortable': True, 'align': 'right'},
                    {'name': 'intensity', 'label': 'Intensity', 'field': 'intensity', 'sortable': True, 'align': 'right'},
                    {'name': 'charge', 'label': 'Z', 'field': 'charge', 'sortable': True, 'align': 'center'},
                    {'name': 'quality', 'label': 'Quality', 'field': 'quality', 'sortable': True, 'align': 'right'},
                ]

                def on_feature_click(e):
                    row = e.args[1]
                    if row and 'idx' in row:
                        viewer.zoom_to_feature(row['idx'])

                viewer.feature_table = ui.table(
                    columns=feature_columns, rows=[], row_key='idx',
                    pagination={'rowsPerPage': 8, 'sortBy': 'intensity', 'descending': True}
                ).classes('w-full').on('rowClick', on_feature_click)
                viewer.feature_table.props('dark flat bordered dense')

            # ID Table
            with ui.card().classes('flex-1 min-w-96'):
                ui.label('Identifications').classes('text-xl font-semibold mb-2')
                ui.label('Click a row to zoom and view annotated spectrum').classes('text-sm text-gray-400 mb-2')

                id_columns = [
                    {'name': 'idx', 'label': '#', 'field': 'idx', 'sortable': True, 'align': 'left'},
                    {'name': 'rt', 'label': 'RT (s)', 'field': 'rt', 'sortable': True, 'align': 'right'},
                    {'name': 'mz', 'label': 'm/z', 'field': 'mz', 'sortable': True, 'align': 'right'},
                    {'name': 'sequence', 'label': 'Sequence', 'field': 'sequence', 'sortable': True, 'align': 'left'},
                    {'name': 'charge', 'label': 'Z', 'field': 'charge', 'sortable': True, 'align': 'center'},
                    {'name': 'score', 'label': 'Score', 'field': 'score', 'sortable': True, 'align': 'right'},
                ]

                def on_id_click(e):
                    row = e.args[1]
                    if row and 'idx' in row:
                        viewer.zoom_to_id(row['idx'])

                viewer.id_table = ui.table(
                    columns=id_columns, rows=[], row_key='idx',
                    pagination={'rowsPerPage': 8, 'sortBy': 'score', 'descending': True}
                ).classes('w-full').on('rowClick', on_id_click)
                viewer.id_table.props('dark flat bordered dense')

        # Spectrum Table (for browsing all spectra)
        with ui.expansion('Spectrum Table', icon='list').classes('w-full max-w-6xl mt-2'):
            ui.label('Click a row to view the spectrum in the 1D viewer above').classes('text-sm text-gray-400 mb-2')

            spectrum_columns = [
                {'name': 'idx', 'label': '#', 'field': 'idx', 'sortable': True, 'align': 'left'},
                {'name': 'rt', 'label': 'RT (s)', 'field': 'rt', 'sortable': True, 'align': 'right'},
                {'name': 'ms_level', 'label': 'MS', 'field': 'ms_level', 'sortable': True, 'align': 'center'},
                {'name': 'n_peaks', 'label': 'Peaks', 'field': 'n_peaks', 'sortable': True, 'align': 'right'},
                {'name': 'tic', 'label': 'TIC', 'field': 'tic', 'sortable': True, 'align': 'right'},
                {'name': 'mz_range', 'label': 'm/z Range', 'field': 'mz_range', 'sortable': False, 'align': 'center'},
                {'name': 'precursor_mz', 'label': 'Prec m/z', 'field': 'precursor_mz', 'sortable': True, 'align': 'right'},
                {'name': 'precursor_z', 'label': 'Prec Z', 'field': 'precursor_z', 'sortable': True, 'align': 'center'},
            ]

            def on_spectrum_click(e):
                row = e.args[1]
                if row and 'idx' in row:
                    viewer.show_spectrum_in_browser(row['idx'])

            viewer.spectrum_table = ui.table(
                columns=spectrum_columns, rows=viewer.spectrum_data, row_key='idx',
                pagination={'rowsPerPage': 10, 'sortBy': 'idx', 'descending': False}
            ).classes('w-full').on('rowClick', on_spectrum_click)
            viewer.spectrum_table.props('dark flat bordered dense')

        # Custom range
        with ui.expansion('Custom Range', icon='tune').classes('w-full max-w-4xl mt-4'):
            with ui.row().classes('w-full gap-4 items-end'):
                rt_min_input = ui.number(label='RT Min (s)', value=0, format='%.2f')
                rt_max_input = ui.number(label='RT Max (s)', value=1000, format='%.2f')
                mz_min_input = ui.number(label='m/z Min', value=0, format='%.2f')
                mz_max_input = ui.number(label='m/z Max', value=2000, format='%.2f')

                def apply_range():
                    viewer.apply_custom_range(
                        rt_min_input.value, rt_max_input.value,
                        mz_min_input.value, mz_max_input.value
                    )

                ui.button('Apply Range', on_click=apply_range).props('color=primary')

        # Legend
        with ui.expansion('Legend & Help', icon='help').classes('w-full max-w-4xl mt-2'):
            with ui.row().classes('gap-8 flex-wrap'):
                with ui.column():
                    ui.label('Overlay Colors:').classes('font-semibold')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:#00ff64;border-radius:50%;border:1px solid white;"></div>', sanitize=False)
                        ui.label('Feature Centroid')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;border:2px solid #ffff00;"></div>', sanitize=False)
                        ui.label('Feature Bounding Box')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:rgba(0,200,255,0.5);border:1px solid #00c8ff;"></div>', sanitize=False)
                        ui.label('Feature Convex Hull')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:#ff9632;transform:rotate(45deg);"></div>', sanitize=False)
                        ui.label('ID Precursor Position')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:#ff64ff;border-radius:50%;"></div>', sanitize=False)
                        ui.label('Selected Item')

                with ui.column():
                    ui.label('Spectrum Annotation:').classes('font-semibold')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:#1f77b4;"></div>', sanitize=False)
                        ui.label('b-ions (blue)')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:#d62728;"></div>', sanitize=False)
                        ui.label('y-ions (red)')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:gray;"></div>', sanitize=False)
                        ui.label('Unmatched peaks')

                with ui.column():
                    ui.label('TIC & Spectra:').classes('font-semibold')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:4px;background:#00d4ff;"></div>', sanitize=False)
                        ui.label('TIC trace')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:rgba(255,255,0,0.2);border:1px solid rgba(255,255,0,0.5);"></div>', sanitize=False)
                        ui.label('Current view range')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:#00ff64;"></div>', sanitize=False)
                        ui.label('MS1 spectrum peaks')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:2px;height:16px;background:#00d4ff;"></div>', sanitize=False)
                        ui.label('MS1 spectrum marker (cyan)')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:2px;height:16px;background:#ff6b6b;"></div>', sanitize=False)
                        ui.label('MS2 spectrum marker (red)')

                with ui.column():
                    ui.label('Keyboard Shortcuts:').classes('font-semibold')
                    ui.markdown('''
| Key | Action |
|-----|--------|
| `+` / `=` | Zoom In |
| `-` | Zoom Out |
| `Arrow Keys` | Pan |
| `Home` | Reset View |
                    ''')

        # Keyboard handlers
        ui.keyboard(
            on_key=lambda e: (
                viewer.zoom_in() if e.key in ['+', '='] and e.action.keydown else
                viewer.zoom_out() if e.key == '-' and e.action.keydown else
                viewer.pan(rt_frac=-0.1) if e.key.arrow_left and e.action.keydown else
                viewer.pan(rt_frac=0.1) if e.key.arrow_right and e.action.keydown else
                viewer.pan(mz_frac=0.1) if e.key.arrow_up and e.action.keydown else
                viewer.pan(mz_frac=-0.1) if e.key.arrow_down and e.action.keydown else
                viewer.reset_view() if e.key == 'Home' and e.action.keydown else
                None
            )
        )

    # Load CLI files after UI is ready
    if _cli_files['mzml']:
        if viewer.load_mzml(_cli_files['mzml']):
            viewer.update_plot()
            viewer.update_tic_plot()
    if _cli_files['featurexml']:
        if viewer.load_featuremap(_cli_files['featurexml']):
            viewer.update_plot()
    if _cli_files['idxml']:
        if viewer.load_idxml(_cli_files['idxml']):
            viewer.update_plot()


@click.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
@click.option('--port', '-p', default=8080, help='Port to run the server on')
@click.option('--host', '-H', default='0.0.0.0', help='Host to bind to')
def main(files, port, host):
    """
    mzML Peak Map Viewer - Fast visualization of mass spectrometry data.

    Pass one or more files to load them automatically:

    \b
    Examples:
        mzml_viewer.py                              # Start empty
        mzml_viewer.py sample.mzML                  # Load mzML
        mzml_viewer.py sample.mzML features.featureXML
        mzml_viewer.py sample.mzML ids.idXML
        mzml_viewer.py data.mzML features.featureXML ids.idXML

    Supported file types (detected by extension):
        .mzML       Mass spectrometry peak data
        .featureXML Detected features with convex hulls
        .idXML      Peptide identifications
    """
    global _cli_files

    for filepath in files:
        path = Path(filepath)
        ext = path.suffix.lower()

        if ext == '.mzml':
            _cli_files['mzml'] = str(path)
            click.echo(f"Will load mzML: {path.name}")
        elif ext == '.featurexml':
            _cli_files['featurexml'] = str(path)
            click.echo(f"Will load featureXML: {path.name}")
        elif ext == '.idxml':
            _cli_files['idxml'] = str(path)
            click.echo(f"Will load idXML: {path.name}")
        elif ext == '.xml':
            name_lower = path.name.lower()
            if 'feature' in name_lower:
                _cli_files['featurexml'] = str(path)
                click.echo(f"Will load as featureXML: {path.name}")
            elif 'id' in name_lower:
                _cli_files['idxml'] = str(path)
                click.echo(f"Will load as idXML: {path.name}")
            else:
                click.echo(f"Unknown XML file type: {path.name} (skipping)")
        else:
            click.echo(f"Unknown file type: {path.name} (skipping)")

    click.echo(f"\nStarting server at http://{host}:{port}")

    @ui.page('/')
    def index():
        create_ui()

    ui.run(
        title='mzML Peak Map Viewer',
        host=host,
        port=port,
        reload=False,
        show=False
    )


if __name__ in {"__main__", "__mp_main__"}:
    main()
