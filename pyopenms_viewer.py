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
import colorcet as cc

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
from nicegui import ui, app, run


# Global viewer instance for CLI file loading
_viewer_instance = None
_cli_files = {'mzml': None, 'featurexml': None, 'idxml': None}


async def file_picker_dialog(title: str = "Select File",
                             start_path: str = None,
                             extensions: List[str] = None) -> Optional[str]:
    """Show a file picker dialog and return the selected file path."""
    import os

    if start_path is None:
        start_path = os.getcwd()

    current_path = start_path
    selected_file = None

    with ui.dialog() as dialog, ui.card().classes('w-96'):
        ui.label(title).classes('text-lg font-bold mb-2')

        # Current path display
        path_label = ui.label(current_path).classes('text-xs text-gray-400 mb-2 break-all')

        # File list container
        file_list = ui.column().classes('w-full max-h-80 overflow-auto border rounded p-2')

        # Selected file display
        selected_label = ui.label('No file selected').classes('text-sm text-cyan-400 mt-2')

        def refresh_file_list():
            nonlocal current_path
            file_list.clear()
            path_label.set_text(current_path)

            try:
                entries = sorted(os.listdir(current_path))
            except PermissionError:
                with file_list:
                    ui.label('Permission denied').classes('text-red-400')
                return

            with file_list:
                # Parent directory
                if current_path != '/':
                    def go_up():
                        nonlocal current_path
                        current_path = str(Path(current_path).parent)
                        refresh_file_list()

                    ui.button('📁 ..', on_click=go_up).props('flat dense align=left').classes('w-full justify-start text-yellow-400')

                # Directories first
                for entry in entries:
                    if entry.startswith('.'):
                        continue
                    full_path = os.path.join(current_path, entry)

                    if os.path.isdir(full_path):
                        def enter_dir(p=full_path):
                            nonlocal current_path
                            current_path = p
                            refresh_file_list()

                        ui.button(f'📁 {entry}', on_click=enter_dir).props('flat dense align=left').classes('w-full justify-start text-yellow-400')

                # Then files
                for entry in entries:
                    if entry.startswith('.'):
                        continue
                    full_path = os.path.join(current_path, entry)

                    if os.path.isfile(full_path):
                        # Filter by extension if specified
                        if extensions:
                            ext = Path(entry).suffix.lower()
                            if ext not in [e.lower() for e in extensions]:
                                continue

                        def select_file(p=full_path, name=entry):
                            nonlocal selected_file
                            selected_file = p
                            selected_label.set_text(f'Selected: {name}')

                        ui.button(f'📄 {entry}', on_click=select_file).props('flat dense align=left').classes('w-full justify-start text-gray-300 hover:text-white')

        refresh_file_list()

        with ui.row().classes('w-full justify-end gap-2 mt-4'):
            ui.button('Cancel', on_click=lambda: dialog.submit(None)).props('flat')
            ui.button('Open', on_click=lambda: dialog.submit(selected_file)).props('color=primary')

    result = await dialog
    return result

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

    # Add experimental spectrum as vertical lines (stem plot) - doesn't get thicker when zooming
    x_stems = []
    y_stems = []
    for mz, intensity in zip(exp_mz, exp_int_norm):
        x_stems.extend([mz, mz, None])
        y_stems.extend([0, intensity, None])

    fig.add_trace(go.Scatter(
        x=x_stems,
        y=y_stems,
        mode='lines',
        line=dict(color='gray', width=1),
        name='Experimental',
        hoverinfo='skip',
        opacity=0.6
    ))

    # Add hover points for experimental peaks
    fig.add_trace(go.Scatter(
        x=exp_mz,
        y=exp_int_norm,
        mode='markers',
        marker=dict(color='gray', size=2),
        showlegend=False,
        hovertemplate='m/z: %{x:.4f}<br>Intensity: %{y:.1f}%<extra></extra>'
    ))

    # Try to generate theoretical spectrum for annotation
    try:
        seq = AASequence.fromString(sequence_str)
        theo_ions = generate_theoretical_spectrum(seq, charge)

        # Match theoretical to experimental and annotate
        matched_peaks = {'b': [], 'y': []}  # Group by ion type

        for ion_type, ions in [('b', theo_ions['b']), ('y', theo_ions['y'])]:
            for theo_mz, ion_name in ions:
                # Find closest experimental peak
                if len(exp_mz) > 0:
                    diffs = np.abs(exp_mz - theo_mz)
                    min_idx = np.argmin(diffs)
                    if diffs[min_idx] <= tolerance_da:
                        matched_peaks[ion_type].append({
                            'mz': exp_mz[min_idx],
                            'intensity': exp_int_norm[min_idx],
                            'label': ion_name
                        })

        # Add matched peaks as colored lines grouped by ion type
        for ion_type, peaks in matched_peaks.items():
            if not peaks:
                continue
            color = ION_COLORS[ion_type]

            # Create stem plot for this ion type
            x_ions = []
            y_ions = []
            for peak in peaks:
                x_ions.extend([peak['mz'], peak['mz'], None])
                y_ions.extend([0, peak['intensity'], None])

            fig.add_trace(go.Scatter(
                x=x_ions,
                y=y_ions,
                mode='lines',
                line=dict(color=color, width=2),
                name=f'{ion_type}-ions',
                hoverinfo='skip'
            ))

            # Add hover points and annotations for matched peaks
            for peak in peaks:
                fig.add_trace(go.Scatter(
                    x=[peak['mz']],
                    y=[peak['intensity']],
                    mode='markers',
                    marker=dict(color=color, size=4),
                    showlegend=False,
                    hovertemplate=f"{peak['label']}<br>m/z: {peak['mz']:.4f}<br>Intensity: {peak['intensity']:.1f}%<extra></extra>"
                ))

                # Add text annotation
                fig.add_annotation(
                    x=peak['mz'],
                    y=peak['intensity'] + 3,
                    text=peak['label'],
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
        )
    )

    fig.update_xaxes(range=[0, max(exp_mz) * 1.05] if len(exp_mz) > 0 else [0, 2000])
    fig.update_yaxes(range=[0, 110])

    return fig


# Available colormaps for peak map visualization
import matplotlib
COLORMAPS = {
    'jet': matplotlib.colormaps['jet'],
    'hot': matplotlib.colormaps['hot'],
    'fire': cc.fire,
    'viridis': matplotlib.colormaps['viridis'],
    'plasma': matplotlib.colormaps['plasma'],
    'inferno': matplotlib.colormaps['inferno'],
    'magma': matplotlib.colormaps['magma'],
}


def get_colormap_background(colormap_name: str) -> str:
    """Get the lowest color from a colormap as a hex string for background."""
    cmap = COLORMAPS.get(colormap_name)
    if cmap is None:
        return 'black'

    # Check if it's a matplotlib colormap (has __call__ method)
    if hasattr(cmap, '__call__'):
        # Matplotlib colormap - get color at 0
        rgba = cmap(0)
        # Convert to hex
        r, g, b = int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)
        return f'#{r:02x}{g:02x}{b:02x}'
    elif isinstance(cmap, list) and len(cmap) > 0:
        # Colorcet list - first element is the lowest color
        return cmap[0]
    else:
        return 'black'


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
        self.show_bounding_boxes = False  # Disabled by default for faster rendering
        self.show_convex_hulls = False    # Disabled by default for faster rendering
        self.show_ids = True
        self.show_spectrum_marker = True  # Always show RT/m/z marker for selected spectrum
        self.colormap = 'jet'  # Default colormap

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
        self.tic_plot = None

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

        # Navigation UI elements
        self.minimap_image = None  # Overview image
        self.minimap_width = 400
        self.minimap_height = 200
        self.breadcrumb_label = None  # Zoom history display
        self.coord_label = None  # RT/m/z coordinate display

        # Zoom history for breadcrumb trail
        self.zoom_history = []  # List of (rt_min, rt_max, mz_min, mz_max, label) tuples
        self.max_zoom_history = 10  # Max history entries

        # UI update flags
        self._updating_from_tic = False  # Prevent circular TIC updates

        # Hover state for visual feedback
        self.hover_feature_idx = None  # Feature being hovered in table
        self.hover_id_idx = None  # ID being hovered in table
        self._hover_update_pending = False  # Debounce hover updates

        # Loading state
        self.loading_indicator = None
        self.loading_label = None  # Label for progress text
        self.is_loading = False

        # 3D visualization
        self.show_3d_view = False
        self.plot_3d = None  # Plotly 3D plot element
        self.scene_3d_container = None
        self.view_3d_status = None  # Status label for 3D view
        self.max_3d_peaks = 5000  # Limit peaks for 3D performance
        self.rt_threshold_3d = 120.0  # Max RT range for 3D (seconds)
        self.mz_threshold_3d = 50.0  # Max m/z range for 3D

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

    def load_mzml_sync(self, filepath: str) -> bool:
        """Load mzML file synchronously without UI updates (for background thread)."""
        try:
            self.exp = MSExperiment()
            MzMLFile().load(filepath, self.exp)

            n_spectra = self.exp.size()
            total_peaks = sum(spec.size() for spec in self.exp)

            if total_peaks == 0:
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
            return True

        except Exception as e:
            print(f"Error loading mzML: {e}")
            return False

    def load_mzml(self, filepath: str) -> bool:
        """Load mzML file and extract peak data (with UI updates)."""
        try:
            filename = Path(filepath).name
            self.set_loading(True, f"Reading {filename}...")
            if self.status_label:
                self.status_label.set_text(f"Reading {filename}...")
            ui.notify(f"Loading {filepath}...", type="info")

            self.exp = MSExperiment()
            MzMLFile().load(filepath, self.exp)

            n_spectra = self.exp.size()
            self.update_loading_progress(f"Loaded {n_spectra:,} spectra, counting peaks...")

            total_peaks = sum(spec.size() for spec in self.exp)

            if total_peaks == 0:
                ui.notify("No peaks found in file!", type="warning")
                return False

            self.update_loading_progress(f"Found {total_peaks:,} peaks, detecting FAIMS...")

            # First pass: detect FAIMS CVs
            cv_set = set()
            for spec in self.exp:
                if spec.getMSLevel() == 1:
                    cv = self._get_cv_from_spectrum(spec)
                    if cv is not None:
                        cv_set.add(cv)

            self.has_faims = len(cv_set) > 1
            self.faims_cvs = sorted(cv_set) if self.has_faims else []

            self.update_loading_progress(f"Extracting {total_peaks:,} peaks...")

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
            spec_count = 0
            progress_interval = max(1, n_spectra // 20)  # Update progress ~20 times
            for spec in self.exp:
                spec_count += 1
                if spec_count % progress_interval == 0:
                    pct = int(100 * spec_count / n_spectra)
                    self.update_loading_progress(f"Extracting peaks... {pct}% ({idx:,} peaks)")

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

            self.update_loading_progress(f"Building data structures...")

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

            self.update_loading_progress(f"Extracting spectrum metadata...")

            # Extract spectrum metadata for browser
            self.spectrum_data = self._extract_spectrum_data()

            self.update_loading_progress(f"Creating DataFrame ({idx:,} peaks)...")

            # Create main DataFrame
            self.df = pd.DataFrame({
                'rt': rts,
                'mz': mzs,
                'intensity': intensities
            })
            if self.has_faims:
                self.df['cv'] = cvs
            self.df['log_intensity'] = np.log1p(self.df['intensity'])

            self.update_loading_progress(f"Finalizing...")

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

            self.set_loading(False)
            return True

        except Exception as e:
            self.set_loading(False)
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
        """Display a spectrum in the 1D browser view. Shows annotations if matching ID exists."""
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

        # Check if there's a matching peptide ID for annotation
        matching_id_idx = self.find_matching_id_for_spectrum(spectrum_idx)

        if matching_id_idx is not None:
            # Use annotated spectrum display
            pep_id = self.peptide_ids[matching_id_idx]
            hits = pep_id.getHits()
            if hits:
                best_hit = hits[0]
                sequence_str = best_hit.getSequence().toString()
                charge = best_hit.getCharge()
                precursors = spec.getPrecursors()
                prec_mz = precursors[0].getMZ() if precursors else pep_id.getMZ()

                # Create annotated spectrum plot
                fig = create_annotated_spectrum_plot(
                    mz_array, int_array,
                    sequence_str, charge, prec_mz
                )

                # Update title to include spectrum index
                title = f"Spectrum #{spectrum_idx} | {sequence_str} (z={charge}+) | RT={rt:.2f}s"
                fig.update_layout(title=dict(text=title, font=dict(size=14)), height=350)

                # Update info label with ID info
                if self.spectrum_browser_info is not None:
                    self.spectrum_browser_info.set_text(
                        f"RT: {rt:.2f}s | ID: {sequence_str} | Charge: {charge}+ | Precursor: {prec_mz:.4f}"
                    )
            else:
                # No hits, fall back to regular display
                matching_id_idx = None

        if matching_id_idx is None:
            # Regular spectrum display (no annotation)
            # Normalize intensities
            max_int = int_array.max() if len(int_array) > 0 else 1
            int_norm = (int_array / max_int) * 100

            # Create figure
            fig = go.Figure()

            # Color based on MS level
            color = '#00d4ff' if ms_level == 1 else '#ff6b6b'

            # Add spectrum as vertical lines (stem plot) - doesn't get thicker when zooming
            # Create x, y arrays for stem plot: each peak is [mz, mz, None], [0, intensity, None]
            x_stems = []
            y_stems = []
            for mz, intensity in zip(mz_array, int_norm):
                x_stems.extend([mz, mz, None])
                y_stems.extend([0, intensity, None])

            fig.add_trace(go.Scatter(
                x=x_stems,
                y=y_stems,
                mode='lines',
                line=dict(color=color, width=1),
                hoverinfo='skip'
            ))

            # Add hover points at peak tops
            fig.add_trace(go.Scatter(
                x=mz_array,
                y=int_norm,
                mode='markers',
                marker=dict(color=color, size=3),
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

            # Update info label
            if self.spectrum_browser_info is not None:
                tic = float(np.sum(int_array))
                mz_range = f"{mz_array.min():.2f} - {mz_array.max():.2f}" if len(mz_array) > 0 else "-"
                self.spectrum_browser_info.set_text(
                    f"RT: {rt:.2f}s | MS Level: {ms_level} | Peaks: {len(mz_array):,} | TIC: {tic:.2e} | m/z: {mz_range}"
                )

        # Update plot
        if self.spectrum_browser_plot is not None:
            self.spectrum_browser_plot.update_figure(fig)

        # Update navigation label
        if self.spectrum_nav_label is not None:
            self.spectrum_nav_label.set_text(f"Spectrum {spectrum_idx + 1} of {self.exp.size()}")

        # Update spectrum table selection
        if self.spectrum_table is not None:
            # Find the row data for this spectrum
            matching_rows = [row for row in self.spectrum_data if row['idx'] == spectrum_idx]
            if matching_rows:
                self.spectrum_table.selected = matching_rows

        # Update TIC to show spectrum marker
        self.update_tic_plot()

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

    def load_featuremap_sync(self, filepath: str) -> bool:
        """Load featureXML file synchronously without UI updates."""
        try:
            self.feature_map = FeatureMap()
            FeatureXMLFile().load(filepath, self.feature_map)
            self.features_file = filepath
            self.selected_feature_idx = None
            self.feature_data = self._extract_feature_data()
            return True
        except Exception as e:
            print(f"Error loading features: {e}")
            return False

    def load_featuremap(self, filepath: str) -> bool:
        """Load featureXML file (with UI updates)."""
        try:
            self.set_loading(True, f"Loading features...")
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

            self.set_loading(False)
            return True

        except Exception as e:
            self.set_loading(False)
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

    def load_idxml_sync(self, filepath: str) -> bool:
        """Load idXML file synchronously without UI updates."""
        try:
            self.protein_ids = []
            self.peptide_ids = []
            IdXMLFile().load(filepath, self.protein_ids, self.peptide_ids)
            self.id_file = filepath
            self.selected_id_idx = None
            self.id_data = self._extract_id_data()
            return True
        except Exception as e:
            print(f"Error loading IDs: {e}")
            return False

    def load_idxml(self, filepath: str) -> bool:
        """Load idXML file with peptide identifications (with UI updates)."""
        try:
            self.set_loading(True, f"Loading identifications...")
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

            self.set_loading(False)
            return True

        except Exception as e:
            self.set_loading(False)
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

    def find_matching_id_for_spectrum(self, spectrum_idx: int, rt_tolerance: float = 5.0, mz_tolerance: float = 0.5) -> Optional[int]:
        """Find peptide ID matching the given spectrum. Returns ID index or None."""
        if self.exp is None or not self.peptide_ids:
            return None

        if spectrum_idx < 0 or spectrum_idx >= self.exp.size():
            return None

        spec = self.exp[spectrum_idx]

        # Only MS2 spectra can have peptide IDs
        if spec.getMSLevel() != 2:
            return None

        spec_rt = spec.getRT()
        precursors = spec.getPrecursors()
        if not precursors:
            return None

        spec_prec_mz = precursors[0].getMZ()

        # Find best matching ID
        best_id_idx = None
        best_rt_diff = float('inf')

        for i, pep_id in enumerate(self.peptide_ids):
            id_rt = pep_id.getRT()
            id_mz = pep_id.getMZ()

            if abs(id_rt - spec_rt) <= rt_tolerance and abs(id_mz - spec_prec_mz) <= mz_tolerance:
                rt_diff = abs(id_rt - spec_rt)
                if rt_diff < best_rt_diff:
                    best_rt_diff = rt_diff
                    best_id_idx = i

        return best_id_idx

    def find_spectrum_for_id(self, id_idx: int, rt_tolerance: float = 5.0, mz_tolerance: float = 0.5) -> Optional[int]:
        """Find spectrum index matching the given peptide ID. Returns spectrum index or None."""
        if self.exp is None or not self.peptide_ids:
            return None

        if id_idx < 0 or id_idx >= len(self.peptide_ids):
            return None

        pep_id = self.peptide_ids[id_idx]
        id_rt = pep_id.getRT()
        id_mz = pep_id.getMZ()

        best_spec_idx = None
        best_rt_diff = float('inf')

        for i in range(self.exp.size()):
            spec = self.exp[i]
            if spec.getMSLevel() != 2:
                continue

            spec_rt = spec.getRT()
            if abs(spec_rt - id_rt) > rt_tolerance:
                continue

            precursors = spec.getPrecursors()
            if precursors:
                prec_mz = precursors[0].getMZ()
                if abs(prec_mz - id_mz) <= mz_tolerance:
                    rt_diff = abs(spec_rt - id_rt)
                    if rt_diff < best_rt_diff:
                        best_rt_diff = rt_diff
                        best_spec_idx = i

        return best_spec_idx

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

        # Add vertical marker for selected spectrum
        if self.selected_spectrum_idx is not None and self.exp is not None:
            spec = self.exp[self.selected_spectrum_idx]
            selected_rt = spec.getRT()
            ms_level = spec.getMSLevel()
            marker_color = '#00d4ff' if ms_level == 1 else '#ff6b6b'
            fig.add_vline(
                x=selected_rt,
                line_dash="solid",
                line_color=marker_color,
                line_width=2,
                annotation_text=f"#{self.selected_spectrum_idx}",
                annotation_position="top",
                annotation_font_color=marker_color,
                annotation_font_size=10
            )

        fig.update_layout(
            title=dict(text="Total Ion Chromatogram (TIC) - Click to select spectrum", font=dict(size=14)),
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

    def find_spectrum_idx_at_rt(self, target_rt: float, ms_level: Optional[int] = None) -> Optional[int]:
        """Find the spectrum index closest to the given RT, optionally filtered by MS level."""
        if self.exp is None:
            return None

        best_idx = None
        best_rt_diff = float('inf')

        for i in range(self.exp.size()):
            spec = self.exp[i]
            if ms_level is not None and spec.getMSLevel() != ms_level:
                continue

            spec_rt = spec.getRT()
            rt_diff = abs(spec_rt - target_rt)
            if rt_diff < best_rt_diff:
                best_rt_diff = rt_diff
                best_idx = i

        return best_idx

    def show_spectrum_at_rt(self, rt: float, ms_level: Optional[int] = None):
        """Display the closest spectrum at the given retention time using the spectrum browser."""
        if self.exp is None:
            ui.notify("Load mzML file first", type="warning")
            return

        spec_idx = self.find_spectrum_idx_at_rt(rt, ms_level=ms_level)
        if spec_idx is None:
            level_str = f"MS{ms_level} " if ms_level else ""
            ui.notify(f"No {level_str}spectrum found near RT={rt:.1f}s", type="warning")
            return

        self.show_spectrum_in_browser(spec_idx)

    def show_ms1_spectrum(self, rt: float):
        """Display MS1 spectrum at the given retention time using the spectrum browser."""
        self.show_spectrum_at_rt(rt, ms_level=1)

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

        # Find matching spectrum and show it with annotations
        spec_idx = self.find_spectrum_for_id(id_idx)
        if spec_idx is not None:
            self.show_spectrum_in_browser(spec_idx)
        else:
            ui.notify(f"No matching MS2 spectrum found for this ID", type="warning")

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
        """Draw a crosshair at the selected spectrum's RT and precursor m/z (for MS2)."""
        if self.selected_spectrum_idx is None or self.exp is None:
            return img

        spec = self.exp[self.selected_spectrum_idx]
        rt = spec.getRT()
        ms_level = spec.getMSLevel()

        # Get precursor m/z for MS2 spectra
        precursor_mz = None
        if ms_level == 2:
            precursors = spec.getPrecursors()
            if precursors:
                precursor_mz = precursors[0].getMZ()

        # Check if RT is in view (still render if precursor_mz is in view)
        rt_in_view = self.view_rt_min <= rt <= self.view_rt_max
        mz_in_view = precursor_mz is not None and self.view_mz_min <= precursor_mz <= self.view_mz_max

        if not rt_in_view and not mz_in_view:
            return img

        img = img.convert('RGBA')
        overlay = Image.new('RGBA', img.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)

        # Colors: cyan for MS1, red for MS2
        line_color = (0, 212, 255, 200) if ms_level == 1 else (255, 107, 107, 200)
        crosshair_color = (255, 200, 50, 180)  # Yellow-orange for crosshair intersection

        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 10)
        except:
            font = ImageFont.load_default()

        # Calculate x position for the RT
        x, _ = self._data_to_plot_pixel(rt, self.view_mz_min)

        # Draw vertical line at RT (if in view)
        if rt_in_view:
            draw.line([(x, 0), (x, self.plot_height)], fill=line_color, width=2)

            # Draw label at top
            label = f"MS{ms_level} #{self.selected_spectrum_idx}"
            draw.text((x + 4, 4), label, fill=line_color, font=font)

        # For MS2 spectra, draw horizontal line at precursor m/z to create crosshair
        if precursor_mz is not None and mz_in_view:
            _, y = self._data_to_plot_pixel(rt, precursor_mz)
            draw.line([(0, y), (self.plot_width, y)], fill=line_color, width=2)

            # Draw precursor m/z label on the right
            mz_label = f"Prec: {precursor_mz:.4f}"
            bbox = draw.textbbox((0, 0), mz_label, font=font)
            label_width = bbox[2] - bbox[0]
            draw.text((self.plot_width - label_width - 4, y - 14), mz_label, fill=line_color, font=font)

            # Draw crosshair intersection marker (if both RT and m/z are in view)
            if rt_in_view:
                # Draw a small circle/crosshair at intersection
                r = 6
                draw.ellipse([(x - r, y - r), (x + r, y + r)], outline=crosshair_color, width=2)
                # Inner cross for visibility
                draw.line([(x - r - 2, y), (x + r + 2, y)], fill=crosshair_color, width=1)
                draw.line([(x, y - r - 2), (x, y + r + 2)], fill=crosshair_color, width=1)

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

    def _draw_hover_overlay(self, img: Image.Image) -> Image.Image:
        """Draw hover highlights for features and IDs."""
        if self.hover_feature_idx is None and self.hover_id_idx is None:
            return img

        img = img.convert('RGBA')
        overlay = Image.new('RGBA', img.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)

        hover_color = (100, 255, 200, 180)  # Bright cyan-green for hover

        # Draw hover highlight for feature
        if self.hover_feature_idx is not None and self.feature_map is not None:
            try:
                feature = self.feature_map[self.hover_feature_idx]
                rt = feature.getRT()
                mz = feature.getMZ()

                # Get bounding box
                convex_hulls = feature.getConvexHulls()
                if convex_hulls:
                    hull = convex_hulls[0]
                    hull_points = hull.getHullPoints()
                    if hull_points:
                        rt_vals = [p.getX() for p in hull_points]
                        mz_vals = [p.getY() for p in hull_points]
                        rt_min_f, rt_max_f = min(rt_vals), max(rt_vals)
                        mz_min_f, mz_max_f = min(mz_vals), max(mz_vals)
                    else:
                        rt_min_f, rt_max_f = rt - 5, rt + 5
                        mz_min_f, mz_max_f = mz - 0.5, mz + 0.5
                else:
                    rt_min_f, rt_max_f = rt - 5, rt + 5
                    mz_min_f, mz_max_f = mz - 0.5, mz + 0.5

                # Draw bounding box preview
                x1, y1 = self._data_to_plot_pixel(rt_min_f, mz_max_f)
                x2, y2 = self._data_to_plot_pixel(rt_max_f, mz_min_f)

                # Clamp to plot area
                x1 = max(0, min(self.plot_width, x1))
                x2 = max(0, min(self.plot_width, x2))
                y1 = max(0, min(self.plot_height, y1))
                y2 = max(0, min(self.plot_height, y2))

                if x1 != x2 and y1 != y2:
                    # Draw dashed-style bounding box (corners + center cross)
                    draw.rectangle([(x1, y1), (x2, y2)], outline=hover_color, width=3)

                    # Draw centroid marker
                    cx, cy = self._data_to_plot_pixel(rt, mz)
                    r = 8
                    draw.ellipse([(cx - r, cy - r), (cx + r, cy + r)], outline=hover_color, width=3)

            except (IndexError, RuntimeError):
                pass

        # Draw hover highlight for ID
        if self.hover_id_idx is not None and self.peptide_ids:
            try:
                pep_id = self.peptide_ids[self.hover_id_idx]
                rt = pep_id.getRT()
                mz = pep_id.getMZ()

                if self._is_in_view(rt, mz):
                    cx, cy = self._data_to_plot_pixel(rt, mz)

                    # Draw pulsing ring effect
                    for r in [10, 14, 18]:
                        alpha = int(180 * (18 - r) / 8)  # Fade out
                        ring_color = (100, 255, 200, alpha)
                        draw.ellipse([(cx - r, cy - r), (cx + r, cy + r)], outline=ring_color, width=2)

            except (IndexError, RuntimeError):
                pass

        img = Image.alpha_composite(img, overlay)
        return img

    def set_hover_feature(self, idx: Optional[int]):
        """Set the feature being hovered."""
        if idx != self.hover_feature_idx:
            self.hover_feature_idx = idx
            self.hover_id_idx = None  # Clear other hover
            self.update_plot()

    def set_hover_id(self, idx: Optional[int]):
        """Set the ID being hovered."""
        if idx != self.hover_id_idx:
            self.hover_id_idx = idx
            self.hover_feature_idx = None  # Clear other hover
            self.update_plot()

    def clear_hover(self):
        """Clear all hover states."""
        if self.hover_feature_idx is not None or self.hover_id_idx is not None:
            self.hover_feature_idx = None
            self.hover_id_idx = None
            self.update_plot()

    def set_loading(self, loading: bool, message: str = "Loading..."):
        """Set loading state and update indicator."""
        self.is_loading = loading
        if self.loading_indicator:
            if loading:
                self.loading_indicator.style('display: flex;')
                if self.loading_label:
                    self.loading_label.set_text(message)
            else:
                self.loading_indicator.style('display: none;')

    def update_loading_progress(self, message: str):
        """Update the loading progress message."""
        if self.loading_label:
            self.loading_label.set_text(message)
        if self.status_label:
            self.status_label.set_text(message)

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
        img = tf.shade(agg, cmap=COLORMAPS[self.colormap], how='linear')
        # Use dynspread to make points more visible (dynamically adjusts based on density)
        img = tf.dynspread(img, threshold=0.5, max_px=3)
        img = tf.set_background(img, get_colormap_background(self.colormap))

        plot_img = img.to_pil()

        if self.feature_map is not None:
            plot_img = self._draw_features_on_plot(plot_img)

        if self.peptide_ids:
            plot_img = self._draw_ids_on_plot(plot_img)

        if self.show_spectrum_marker:
            plot_img = self._draw_spectrum_marker_on_plot(plot_img)

        # Draw hover highlights (last so they appear on top)
        plot_img = self._draw_hover_overlay(plot_img)

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
        img = tf.shade(agg, cmap=COLORMAPS[self.colormap], how='linear')
        img = tf.dynspread(img, threshold=0.5, max_px=3)
        img = tf.set_background(img, get_colormap_background(self.colormap))

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

        # Update minimap
        self.update_minimap()

        # Update breadcrumb
        self.update_breadcrumb()

        # Update 3D view if enabled
        if self.show_3d_view:
            self.update_3d_view()

    def render_minimap(self) -> Optional[str]:
        """Render the minimap showing full data extent with view rectangle overlay."""
        if self.df is None or len(self.df) == 0:
            return None

        # Create minimap canvas (no margins for simplicity)
        cvs = ds.Canvas(plot_width=self.minimap_width, plot_height=self.minimap_height,
                        x_range=(self.rt_min, self.rt_max),
                        y_range=(self.mz_min, self.mz_max))

        # Aggregate using log of count for better visualization
        agg = cvs.points(self.df, 'rt', 'mz', agg=ds.sum('intensity'))

        # Apply color map
        img = tf.shade(agg, cmap=COLORMAPS[self.colormap], how='log')
        img = tf.dynspread(img, threshold=0.5, max_px=2)
        img = tf.set_background(img, get_colormap_background(self.colormap))

        # Convert to PIL
        plot_img = img.to_pil()

        # Draw view rectangle
        if (self.view_rt_min is not None and self.view_rt_max is not None and
            self.view_mz_min is not None and self.view_mz_max is not None):

            draw = ImageDraw.Draw(plot_img)

            # Convert data coords to pixel coords
            rt_range = self.rt_max - self.rt_min
            mz_range = self.mz_max - self.mz_min

            if rt_range > 0 and mz_range > 0:
                x1 = int((self.view_rt_min - self.rt_min) / rt_range * self.minimap_width)
                x2 = int((self.view_rt_max - self.rt_min) / rt_range * self.minimap_width)
                y1 = int((self.mz_max - self.view_mz_max) / mz_range * self.minimap_height)
                y2 = int((self.mz_max - self.view_mz_min) / mz_range * self.minimap_height)

                # Clamp to minimap bounds
                x1, x2 = max(0, x1), min(self.minimap_width - 1, x2)
                y1, y2 = max(0, y1), min(self.minimap_height - 1, y2)

                # Draw rectangle outline (yellow)
                draw.rectangle([x1, y1, x2, y2], outline=(255, 255, 0, 255), width=2)

        # Convert to base64
        buffer = io.BytesIO()
        plot_img.save(buffer, format='PNG')
        return base64.b64encode(buffer.getvalue()).decode('utf-8')

    def update_minimap(self):
        """Update the minimap display."""
        if self.minimap_image is None:
            return

        img_data = self.render_minimap()
        if img_data:
            self.minimap_image.set_source(f"data:image/png;base64,{img_data}")

    def minimap_click_to_view(self, x_frac: float, y_frac: float):
        """Center the main view on the clicked position in minimap."""
        if self.df is None:
            return

        # Convert minimap fractions to data coordinates
        rt_click = self.rt_min + x_frac * (self.rt_max - self.rt_min)
        mz_click = self.mz_max - y_frac * (self.mz_max - self.mz_min)

        # Center current view on this point
        rt_half_range = (self.view_rt_max - self.view_rt_min) / 2
        mz_half_range = (self.view_mz_max - self.view_mz_min) / 2

        new_rt_min = rt_click - rt_half_range
        new_rt_max = rt_click + rt_half_range
        new_mz_min = mz_click - mz_half_range
        new_mz_max = mz_click + mz_half_range

        # Clamp to data bounds
        if new_rt_min < self.rt_min:
            new_rt_max += (self.rt_min - new_rt_min)
            new_rt_min = self.rt_min
        if new_rt_max > self.rt_max:
            new_rt_min -= (new_rt_max - self.rt_max)
            new_rt_max = self.rt_max

        if new_mz_min < self.mz_min:
            new_mz_max += (self.mz_min - new_mz_min)
            new_mz_min = self.mz_min
        if new_mz_max > self.mz_max:
            new_mz_min -= (new_mz_max - self.mz_max)
            new_mz_max = self.mz_max

        # Final clamp
        self.view_rt_min = max(self.rt_min, new_rt_min)
        self.view_rt_max = min(self.rt_max, new_rt_max)
        self.view_mz_min = max(self.mz_min, new_mz_min)
        self.view_mz_max = min(self.mz_max, new_mz_max)

        self.update_plot()

    def push_zoom_history(self):
        """Save current view state to zoom history."""
        if self.view_rt_min is None:
            return

        # Create label for this view state
        rt_range = self.view_rt_max - self.view_rt_min
        mz_range = self.view_mz_max - self.view_mz_min
        full_rt = self.rt_max - self.rt_min
        full_mz = self.mz_max - self.mz_min

        # Check if this is approximately full view
        if rt_range >= full_rt * 0.95 and mz_range >= full_mz * 0.95:
            label = "Full"
        else:
            label = f"RT {self.view_rt_min:.0f}-{self.view_rt_max:.0f}"

        state = (self.view_rt_min, self.view_rt_max, self.view_mz_min, self.view_mz_max, label)

        # Don't add if same as last entry
        if self.zoom_history and self.zoom_history[-1][:4] == state[:4]:
            return

        self.zoom_history.append(state)

        # Limit history size
        if len(self.zoom_history) > self.max_zoom_history:
            self.zoom_history = self.zoom_history[-self.max_zoom_history:]

    def go_to_zoom_history(self, index: int):
        """Jump to a specific point in zoom history."""
        if index < 0 or index >= len(self.zoom_history):
            return

        state = self.zoom_history[index]
        self.view_rt_min, self.view_rt_max, self.view_mz_min, self.view_mz_max, _ = state

        # Truncate history to this point (forward history is lost)
        self.zoom_history = self.zoom_history[:index + 1]

        self.update_plot()

    def update_breadcrumb(self):
        """Update the breadcrumb trail display."""
        if self.breadcrumb_label is None:
            return

        if not self.zoom_history:
            self.breadcrumb_label.set_text("Full view")
            return

        # Build breadcrumb string
        parts = []
        for i, (_, _, _, _, label) in enumerate(self.zoom_history):
            parts.append(label)

        breadcrumb_text = " → ".join(parts)
        self.breadcrumb_label.set_text(breadcrumb_text)

    def pixel_to_data_coords(self, pixel_x: int, pixel_y: int) -> Tuple[float, float]:
        """Convert pixel coordinates to RT/m/z data coordinates."""
        # Account for margins
        plot_x = pixel_x - self.margin_left
        plot_y = pixel_y - self.margin_top

        # Clamp to plot area
        plot_x = max(0, min(self.plot_width, plot_x))
        plot_y = max(0, min(self.plot_height, plot_y))

        # Convert to data coordinates
        rt = self.view_rt_min + (plot_x / self.plot_width) * (self.view_rt_max - self.view_rt_min)
        mz = self.view_mz_max - (plot_y / self.plot_height) * (self.view_mz_max - self.view_mz_min)

        return rt, mz

    def update_coord_display(self, pixel_x: int, pixel_y: int):
        """Update the coordinate display label."""
        if self.coord_label is None or self.df is None:
            return

        rt, mz = self.pixel_to_data_coords(pixel_x, pixel_y)
        self.coord_label.set_text(f"RT: {rt:.2f}s  m/z: {mz:.4f}")

    # ==================== 3D Visualization Methods ====================

    def is_small_region(self) -> bool:
        """Check if current view is small enough for 3D visualization."""
        if self.view_rt_min is None or self.view_mz_min is None:
            return False
        rt_range = self.view_rt_max - self.view_rt_min
        mz_range = self.view_mz_max - self.view_mz_min
        return rt_range <= self.rt_threshold_3d and mz_range <= self.mz_threshold_3d

    def update_3d_view(self):
        """Update the 3D visualization with current view data using pyopenms-viz."""
        if not self.show_3d_view or self.plot_3d is None or self.df is None:
            return

        # Check if region is small enough
        if not self.is_small_region():
            # Show message that region is too large
            if self.view_3d_status:
                rt_range = self.view_rt_max - self.view_rt_min
                mz_range = self.view_mz_max - self.view_mz_min
                self.view_3d_status.set_text(
                    f'Zoom in more for 3D (current: RT={rt_range:.0f}s, m/z={mz_range:.0f} | need: RT≤{self.rt_threshold_3d:.0f}s, m/z≤{self.mz_threshold_3d:.0f})'
                )
            return

        # Get peaks in current view
        mask = (
            (self.df['rt'] >= self.view_rt_min) &
            (self.df['rt'] <= self.view_rt_max) &
            (self.df['mz'] >= self.view_mz_min) &
            (self.df['mz'] <= self.view_mz_max)
        )
        view_df = self.df[mask].copy()

        if len(view_df) == 0:
            if self.view_3d_status:
                self.view_3d_status.set_text('No peaks in view')
            return

        # Subsample if too many peaks
        num_peaks_total = len(view_df)
        if len(view_df) > self.max_3d_peaks:
            view_df = view_df.nlargest(self.max_3d_peaks, 'intensity')
        num_peaks_shown = len(view_df)

        try:
            # Use pyopenms-viz for 3D plotting
            from pyopenms_viz._plotly.core import PLOTLYPeakMapPlot

            # Rename columns to match pyopenms-viz expectations
            plot_df = view_df.rename(columns={'rt': 'RT', 'mz': 'mz', 'intensity': 'int'})

            # Create 3D plot (no title - header shows info)
            plot = PLOTLYPeakMapPlot(
                plot_df,
                x='RT',
                y='mz',
                z='int',
                plot_3d=True,
                title=''
            )
            plot.plot()

            # Get the plotly figure
            fig = plot.fig

            # Update layout for dark theme - maximize space usage
            fig.update_layout(
                paper_bgcolor='#1a1a1f',
                plot_bgcolor='#1a1a1f',
                font=dict(color='#cccccc'),
                scene=dict(
                    xaxis=dict(title='RT (s)', backgroundcolor='#1a1a1f', gridcolor='#333333'),
                    yaxis=dict(title='m/z', backgroundcolor='#1a1a1f', gridcolor='#333333'),
                    zaxis=dict(title='Intensity', backgroundcolor='#1a1a1f', gridcolor='#333333'),
                    bgcolor='#1a1a1f',
                    aspectmode='manual',
                    aspectratio=dict(x=1.5, y=1, z=0.8)
                ),
                margin=dict(l=0, r=0, t=0, b=0),
                height=500,
                showlegend=True,
                legend=dict(x=0, y=1, bgcolor='rgba(26,26,31,0.8)'),
                modebar=dict(orientation='v', bgcolor='rgba(0,0,0,0)')
            )

            # Add feature markers if available
            if self.feature_map is not None and self.show_centroids:
                self._add_features_to_3d_plot(fig)

            # Update the plotly element
            self.plot_3d.update_figure(fig)

            # Update status
            if self.view_3d_status:
                if num_peaks_shown < num_peaks_total:
                    self.view_3d_status.set_text(f'Showing {num_peaks_shown:,} of {num_peaks_total:,} peaks (top intensity)')
                else:
                    self.view_3d_status.set_text(f'Showing {num_peaks_shown:,} peaks')

        except Exception as e:
            if self.view_3d_status:
                self.view_3d_status.set_text(f'3D plot error: {str(e)[:50]}')

    def _add_features_to_3d_plot(self, fig):
        """Add feature bounding boxes to the 3D plotly figure."""
        import plotly.graph_objects as go

        if self.feature_map is None:
            return

        # Collect all box edges for a single trace (more efficient)
        box_x = []
        box_y = []
        box_z = []

        for i, feature in enumerate(self.feature_map):
            rt = feature.getRT()
            mz = feature.getMZ()

            # Check if feature is in current view
            if not (self.view_rt_min <= rt <= self.view_rt_max and
                    self.view_mz_min <= mz <= self.view_mz_max):
                continue

            # Get RT and m/z bounds from convex hull
            hulls = feature.getConvexHulls()
            if hulls and len(hulls) > 0:
                # Get bounds from convex hull points
                hull_points = hulls[0].getHullPoints()
                if len(hull_points) > 0:
                    rt_vals = [p[0] for p in hull_points]
                    mz_vals = [p[1] for p in hull_points]
                    rt_min, rt_max = min(rt_vals), max(rt_vals)
                    mz_min, mz_max = min(mz_vals), max(mz_vals)
                else:
                    # Fallback to small box around centroid
                    rt_min, rt_max = rt - 5, rt + 5
                    mz_min, mz_max = mz - 0.5, mz + 0.5
            else:
                # Fallback to small box around centroid
                rt_min, rt_max = rt - 5, rt + 5
                mz_min, mz_max = mz - 0.5, mz + 0.5

            # Draw box edges on the baseline (z=0)
            # Bottom rectangle (4 edges)
            z_base = 0
            # Edge 1: rt_min to rt_max at mz_min
            box_x.extend([rt_min, rt_max, None])
            box_y.extend([mz_min, mz_min, None])
            box_z.extend([z_base, z_base, None])
            # Edge 2: rt_max at mz_min to mz_max
            box_x.extend([rt_max, rt_max, None])
            box_y.extend([mz_min, mz_max, None])
            box_z.extend([z_base, z_base, None])
            # Edge 3: rt_max to rt_min at mz_max
            box_x.extend([rt_max, rt_min, None])
            box_y.extend([mz_max, mz_max, None])
            box_z.extend([z_base, z_base, None])
            # Edge 4: rt_min at mz_max to mz_min
            box_x.extend([rt_min, rt_min, None])
            box_y.extend([mz_max, mz_min, None])
            box_z.extend([z_base, z_base, None])

        if box_x:
            # Add all bounding boxes as a single trace
            fig.add_trace(go.Scatter3d(
                x=box_x,
                y=box_y,
                z=box_z,
                mode='lines',
                line=dict(color='#00ff66', width=3),
                name='Features',
                hoverinfo='skip'
            ))

    # ==================== Search and Filter Methods ====================

    def search_global(self, query: str) -> List[Dict[str, Any]]:
        """
        Search across spectra, IDs, and features.
        Returns list of results with type, description, and action.
        """
        results = []
        query = query.strip()
        if not query:
            return results

        query_lower = query.lower()

        # Try to parse as spectrum number (e.g., "#123" or "123")
        spec_match = query.replace('#', '').strip()
        if spec_match.isdigit():
            spec_idx = int(spec_match)
            if self.exp and 0 <= spec_idx < self.exp.size():
                spec = self.exp[spec_idx]
                results.append({
                    'type': 'spectrum',
                    'icon': 'analytics',
                    'label': f"Spectrum #{spec_idx} (MS{spec.getMSLevel()}, RT={spec.getRT():.1f}s)",
                    'action': lambda idx=spec_idx: self.show_spectrum_in_browser(idx)
                })

        # Try to parse as m/z value (e.g., "500.25" or "mz:500.25")
        mz_query = query_lower.replace('mz:', '').replace('m/z:', '').strip()
        try:
            mz_val = float(mz_query)
            if 50 < mz_val < 10000:  # Reasonable m/z range
                # Find spectra with precursor near this m/z
                if self.exp:
                    for i in range(min(self.exp.size(), 1000)):  # Limit search
                        spec = self.exp[i]
                        if spec.getMSLevel() > 1:
                            precs = spec.getPrecursors()
                            if precs and abs(precs[0].getMZ() - mz_val) < 0.5:
                                results.append({
                                    'type': 'spectrum',
                                    'icon': 'analytics',
                                    'label': f"MS2 #{i} precursor m/z={precs[0].getMZ():.4f}",
                                    'action': lambda idx=i: self.show_spectrum_in_browser(idx)
                                })
                                if len(results) >= 10:
                                    break
                # Find IDs near this m/z
                for i, pid in enumerate(self.peptide_ids[:100]):
                    if abs(pid.getMZ() - mz_val) < 0.5:
                        hits = pid.getHits()
                        seq = hits[0].getSequence().toString() if hits else "?"
                        results.append({
                            'type': 'id',
                            'icon': 'biotech',
                            'label': f"ID: {seq} (m/z={pid.getMZ():.4f})",
                            'action': lambda idx=i: self.zoom_to_id(idx)
                        })
                        if len(results) >= 15:
                            break
        except ValueError:
            pass

        # Try to parse as RT value (e.g., "rt:100" or "100s")
        rt_query = query_lower.replace('rt:', '').replace('s', '').strip()
        try:
            rt_val = float(rt_query)
            if self.rt_min <= rt_val <= self.rt_max:
                results.append({
                    'type': 'navigate',
                    'icon': 'place',
                    'label': f"Go to RT={rt_val:.1f}s",
                    'action': lambda rt=rt_val: self.go_to_rt(rt)
                })
        except ValueError:
            pass

        # Search peptide sequences
        if len(query) >= 2:
            for i, pid in enumerate(self.peptide_ids[:500]):
                hits = pid.getHits()
                if hits:
                    seq = hits[0].getSequence().toString()
                    if query_lower in seq.lower():
                        results.append({
                            'type': 'id',
                            'icon': 'biotech',
                            'label': f"{seq} (RT={pid.getRT():.1f}s, m/z={pid.getMZ():.2f})",
                            'action': lambda idx=i: self.zoom_to_id(idx)
                        })
                        if len(results) >= 20:
                            break

        return results[:20]  # Limit total results

    def go_to_rt(self, rt: float):
        """Center view on given RT value."""
        if self.df is None:
            return
        self.push_zoom_history()
        rt_range = self.view_rt_max - self.view_rt_min
        self.view_rt_min = max(self.rt_min, rt - rt_range / 2)
        self.view_rt_max = min(self.rt_max, rt + rt_range / 2)
        self.push_zoom_history()
        self.update_plot()
        # Also show closest spectrum
        self.show_spectrum_at_rt(rt)

    def go_to_mz(self, mz: float):
        """Center view on given m/z value."""
        if self.df is None:
            return
        self.push_zoom_history()
        mz_range = self.view_mz_max - self.view_mz_min
        self.view_mz_min = max(self.mz_min, mz - mz_range / 2)
        self.view_mz_max = min(self.mz_max, mz + mz_range / 2)
        self.push_zoom_history()
        self.update_plot()

    def go_to_location(self, rt: Optional[float] = None, mz: Optional[float] = None, spectrum_idx: Optional[int] = None):
        """Navigate to a specific location."""
        if spectrum_idx is not None and self.exp and 0 <= spectrum_idx < self.exp.size():
            self.show_spectrum_in_browser(spectrum_idx)
            spec = self.exp[spectrum_idx]
            rt = spec.getRT()

        if self.df is None:
            return

        self.push_zoom_history()

        if rt is not None:
            rt_range = self.view_rt_max - self.view_rt_min
            self.view_rt_min = max(self.rt_min, rt - rt_range / 2)
            self.view_rt_max = min(self.rt_max, rt + rt_range / 2)

        if mz is not None:
            mz_range = self.view_mz_max - self.view_mz_min
            self.view_mz_min = max(self.mz_min, mz - mz_range / 2)
            self.view_mz_max = min(self.mz_max, mz + mz_range / 2)

        self.push_zoom_history()
        self.update_plot()

    def filter_spectrum_data(self, ms_level: Optional[int] = None, rt_min: Optional[float] = None,
                             rt_max: Optional[float] = None, min_peaks: Optional[int] = None,
                             min_tic: Optional[float] = None) -> List[Dict]:
        """Filter spectrum data based on criteria."""
        filtered = []
        for row in self.spectrum_data:
            if ms_level is not None and row['ms_level'] != ms_level:
                continue
            if rt_min is not None and row['rt'] < rt_min:
                continue
            if rt_max is not None and row['rt'] > rt_max:
                continue
            if min_peaks is not None and row['n_peaks'] < min_peaks:
                continue
            if min_tic is not None and row['tic'] < min_tic:
                continue
            filtered.append(row)
        return filtered

    def filter_id_data(self, sequence_pattern: Optional[str] = None, min_score: Optional[float] = None,
                       charge: Optional[int] = None) -> List[Dict]:
        """Filter ID data based on criteria."""
        filtered = []
        for row in self.id_data:
            if sequence_pattern and sequence_pattern.lower() not in row['sequence'].lower():
                continue
            if min_score is not None and row['score'] < min_score:
                continue
            if charge is not None and row['charge'] != charge:
                continue
            filtered.append(row)
        return filtered

    def filter_feature_data(self, min_intensity: Optional[float] = None,
                            min_quality: Optional[float] = None, charge: Optional[int] = None) -> List[Dict]:
        """Filter feature data based on criteria."""
        filtered = []
        for row in self.feature_data:
            if min_intensity is not None and row['intensity'] < min_intensity:
                continue
            if min_quality is not None and row['quality'] < min_quality:
                continue
            if charge is not None and row['charge'] != charge:
                continue
            filtered.append(row)
        return filtered

    def reset_view(self):
        """Reset to full view."""
        if self.df is None:
            return
        # Clear zoom history on reset
        self.zoom_history = []
        self.view_rt_min = self.rt_min
        self.view_rt_max = self.rt_max
        self.view_mz_min = self.mz_min
        self.view_mz_max = self.mz_max
        self.selected_feature_idx = None
        self.selected_id_idx = None
        self.push_zoom_history()  # Add full view as first entry
        self.update_plot()

    def zoom_in(self, factor=0.5):
        """Zoom in."""
        if self.df is None:
            return
        self.push_zoom_history()  # Save current state before zoom
        rt_center = (self.view_rt_min + self.view_rt_max) / 2
        mz_center = (self.view_mz_min + self.view_mz_max) / 2
        rt_range = (self.view_rt_max - self.view_rt_min) * factor / 2
        mz_range = (self.view_mz_max - self.view_mz_min) * factor / 2

        self.view_rt_min = rt_center - rt_range
        self.view_rt_max = rt_center + rt_range
        self.view_mz_min = mz_center - mz_range
        self.view_mz_max = mz_center + mz_range
        self.push_zoom_history()  # Save new state
        self.update_plot()

    def zoom_out(self, factor=2.0):
        """Zoom out."""
        if self.df is None:
            return
        self.push_zoom_history()  # Save current state before zoom
        rt_center = (self.view_rt_min + self.view_rt_max) / 2
        mz_center = (self.view_mz_min + self.view_mz_max) / 2
        rt_range = (self.view_rt_max - self.view_rt_min) * factor / 2
        mz_range = (self.view_mz_max - self.view_mz_min) * factor / 2

        self.view_rt_min = max(self.rt_min, rt_center - rt_range)
        self.view_rt_max = min(self.rt_max, rt_center + rt_range)
        self.view_mz_min = max(self.mz_min, mz_center - mz_range)
        self.view_mz_max = min(self.mz_max, mz_center + mz_range)
        self.push_zoom_history()  # Save new state
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

    def zoom_at_point(self, x_frac: float, y_frac: float, zoom_in: bool = True):
        """Zoom centered on a specific point (given as fraction of plot area).

        Args:
            x_frac: Horizontal position (0=left, 1=right) in plot area
            y_frac: Vertical position (0=top, 1=bottom) in plot area
            zoom_in: True to zoom in, False to zoom out
        """
        if self.df is None:
            return

        # Save current state to zoom history
        self.push_zoom_history()

        # Convert fractions to data coordinates
        rt_point = self.view_rt_min + x_frac * (self.view_rt_max - self.view_rt_min)
        mz_point = self.view_mz_max - y_frac * (self.view_mz_max - self.view_mz_min)  # Y is inverted

        # Zoom factor
        factor = 0.7 if zoom_in else 1.4

        # Current ranges
        rt_range = self.view_rt_max - self.view_rt_min
        mz_range = self.view_mz_max - self.view_mz_min

        # New ranges
        new_rt_range = rt_range * factor
        new_mz_range = mz_range * factor

        # Keep the point under cursor at same position
        new_rt_min = rt_point - x_frac * new_rt_range
        new_rt_max = rt_point + (1 - x_frac) * new_rt_range
        new_mz_min = mz_point - (1 - y_frac) * new_mz_range
        new_mz_max = mz_point + y_frac * new_mz_range

        # Clamp to data bounds
        self.view_rt_min = max(self.rt_min, new_rt_min)
        self.view_rt_max = min(self.rt_max, new_rt_max)
        self.view_mz_min = max(self.mz_min, new_mz_min)
        self.view_mz_max = min(self.mz_max, new_mz_max)

        self.update_plot()

    def pan_by_pixels(self, dx: float, dy: float):
        """Pan the view by pixel amounts.

        Args:
            dx: Horizontal pixel delta (positive = pan right/increase RT)
            dy: Vertical pixel delta (positive = pan down/decrease mz)
        """
        if self.df is None:
            return

        # Convert pixels to data units
        rt_per_pixel = (self.view_rt_max - self.view_rt_min) / self.plot_width
        mz_per_pixel = (self.view_mz_max - self.view_mz_min) / self.plot_height

        rt_shift = -dx * rt_per_pixel  # Negative because dragging right should decrease RT view
        mz_shift = dy * mz_per_pixel   # Positive because dragging down should decrease mz view

        # Clamp shifts to stay within bounds
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


def create_ui():
    """Create NiceGUI interface."""
    global _viewer_instance

    viewer = MzMLViewer()
    _viewer_instance = viewer

    ui.dark_mode().enable()

    with ui.column().classes('w-full items-center p-4'):
        ui.label('mzML Peak Map Viewer').classes('text-3xl font-bold mb-2')
        ui.label('High-performance visualization with Datashader + pyOpenMS').classes('text-gray-400 mb-4')

        # File loading section (local filesystem paths)
        with ui.card().classes('w-full max-w-6xl mb-4'):
            ui.label('Load Data').classes('text-xl font-semibold mb-2')

            with ui.row().classes('w-full items-end gap-4 flex-wrap'):
                # mzML
                with ui.column().classes('flex-1 min-w-64'):
                    ui.label('mzML File (Peak Data)').classes('text-sm text-gray-400')
                    with ui.row().classes('w-full items-end gap-2'):
                        mzml_input = ui.input(placeholder='/path/to/file.mzML').classes('flex-1')

                        async def browse_mzml():
                            path = await file_picker_dialog(
                                title="Select mzML File",
                                extensions=['.mzml', '.mzML']
                            )
                            if path:
                                mzml_input.value = path

                        ui.button('Browse', on_click=browse_mzml).props('dense outline')

                        async def load_mzml_path():
                            path = mzml_input.value
                            if not path:
                                ui.notify("Please enter a file path", type="warning")
                                return
                            if not Path(path).exists():
                                ui.notify("File not found", type="warning")
                                return

                            # Show loading state
                            viewer.set_loading(True, f"Loading {Path(path).name}...")
                            if viewer.status_label:
                                viewer.status_label.set_text(f"Loading {Path(path).name}...")

                            try:
                                # Run heavy I/O in background thread to keep UI responsive
                                success = await run.io_bound(viewer.load_mzml_sync, path)
                                if success:
                                    # Update all UI elements after successful load
                                    viewer.update_plot()
                                    viewer.update_tic_plot()

                                    # Update info labels
                                    info_text = f"Loaded: {Path(path).name} | Spectra: {viewer.exp.size():,} | Peaks: {len(viewer.df):,}"
                                    if viewer.has_faims:
                                        info_text += f" | FAIMS: {len(viewer.faims_cvs)} CVs"
                                    if viewer.info_label:
                                        viewer.info_label.set_text(info_text)

                                    # Update spectrum table
                                    if viewer.spectrum_table is not None:
                                        viewer.spectrum_table.update_rows(viewer.spectrum_data)

                                    # Update FAIMS UI
                                    if viewer.has_faims:
                                        if viewer.faims_info_label:
                                            cv_str = ", ".join([f"{cv:.1f}V" for cv in viewer.faims_cvs])
                                            viewer.faims_info_label.set_text(f"FAIMS CVs: {cv_str}")
                                            viewer.faims_info_label.set_visibility(True)
                                        if viewer.faims_toggle:
                                            viewer.faims_toggle.set_visibility(True)

                                    ui.notify(f"Loaded {len(viewer.df):,} peaks", type="positive")
                                else:
                                    ui.notify("No peaks found in file", type="warning")
                            except Exception as ex:
                                ui.notify(f"Error loading file: {ex}", type="negative")
                            finally:
                                viewer.set_loading(False)
                                if viewer.status_label:
                                    viewer.status_label.set_text("Ready")

                        ui.button('Load', on_click=load_mzml_path).props('color=primary dense')

                # FeatureXML
                with ui.column().classes('flex-1 min-w-64'):
                    ui.label('FeatureXML (Features)').classes('text-sm text-gray-400')
                    with ui.row().classes('w-full items-end gap-2'):
                        feature_input = ui.input(placeholder='/path/to/features.featureXML').classes('flex-1')

                        async def browse_feature():
                            path = await file_picker_dialog(
                                title="Select FeatureXML File",
                                extensions=['.featurexml', '.featureXML', '.xml']
                            )
                            if path:
                                feature_input.value = path

                        ui.button('Browse', on_click=browse_feature).props('dense outline')

                        async def load_feature_path():
                            path = feature_input.value
                            if not path or not Path(path).exists():
                                ui.notify("File not found", type="warning")
                                return

                            viewer.set_loading(True, f"Loading features...")
                            try:
                                success = await run.io_bound(viewer.load_featuremap_sync, path)
                                if success:
                                    viewer.update_plot()
                                    if viewer.feature_info_label:
                                        viewer.feature_info_label.set_text(f"Features: {viewer.feature_map.size():,}")
                                    if viewer.feature_table is not None:
                                        viewer.feature_table.update_rows(viewer.feature_data)
                                    ui.notify(f"Loaded {viewer.feature_map.size():,} features", type="positive")
                            except Exception as ex:
                                ui.notify(f"Error loading features: {ex}", type="negative")
                            finally:
                                viewer.set_loading(False)
                                if viewer.status_label:
                                    viewer.status_label.set_text("Ready")

                        ui.button('Load', on_click=load_feature_path).props('color=primary dense')

                        def clear_features():
                            viewer.clear_features()
                            viewer.update_plot()

                        ui.button('Clear', on_click=clear_features).props('color=negative dense')

                # idXML
                with ui.column().classes('flex-1 min-w-64'):
                    ui.label('idXML (Identifications)').classes('text-sm text-gray-400')
                    with ui.row().classes('w-full items-end gap-2'):
                        id_input = ui.input(placeholder='/path/to/ids.idXML').classes('flex-1')

                        async def browse_id():
                            path = await file_picker_dialog(
                                title="Select idXML File",
                                extensions=['.idxml', '.idXML', '.xml']
                            )
                            if path:
                                id_input.value = path

                        ui.button('Browse', on_click=browse_id).props('dense outline')

                        async def load_id_path():
                            path = id_input.value
                            if not path or not Path(path).exists():
                                ui.notify("File not found", type="warning")
                                return

                            viewer.set_loading(True, f"Loading identifications...")
                            try:
                                success = await run.io_bound(viewer.load_idxml_sync, path)
                                if success:
                                    viewer.update_plot()
                                    if viewer.id_info_label:
                                        viewer.id_info_label.set_text(f"IDs: {len(viewer.peptide_ids):,}")
                                    if viewer.id_table is not None:
                                        viewer.id_table.update_rows(viewer.id_data)
                                    ui.notify(f"Loaded {len(viewer.peptide_ids):,} peptide IDs", type="positive")
                            except Exception as ex:
                                ui.notify(f"Error loading IDs: {ex}", type="negative")
                            finally:
                                viewer.set_loading(False)
                                if viewer.status_label:
                                    viewer.status_label.set_text("Ready")

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

        # Search and Go To section
        with ui.card().classes('w-full max-w-4xl mb-2 p-2'):
            with ui.row().classes('w-full items-center gap-4'):
                # Global search
                with ui.column().classes('flex-grow'):
                    ui.label('Search').classes('text-xs text-gray-400')
                    search_input = ui.input(placeholder='Peptide sequence, m/z, RT, or spectrum #...').classes('w-full').props('dense outlined')

                    # Search results container
                    search_results_container = ui.column().classes('w-full')
                    search_results_container.set_visibility(False)

                    def on_search(e):
                        query = search_input.value
                        search_results_container.clear()
                        if not query or len(query) < 1:
                            search_results_container.set_visibility(False)
                            return

                        results = viewer.search_global(query)
                        if results:
                            search_results_container.set_visibility(True)
                            with search_results_container:
                                with ui.card().classes('w-full p-1 mt-1').style('max-height: 200px; overflow-y: auto;'):
                                    for result in results:
                                        def make_click_handler(action):
                                            def handler():
                                                action()
                                                search_results_container.set_visibility(False)
                                                search_input.value = ''
                                            return handler

                                        with ui.row().classes('w-full items-center p-1 hover:bg-gray-700 cursor-pointer rounded').on('click', make_click_handler(result['action'])):
                                            ui.icon(result['icon'], size='xs').classes('text-gray-400 mr-2')
                                            ui.label(result['label']).classes('text-sm')
                        else:
                            search_results_container.set_visibility(False)

                    search_input.on('keyup', on_search)

                # Divider
                ui.label('|').classes('text-gray-600')

                # Go To controls
                with ui.row().classes('items-end gap-2'):
                    goto_rt = ui.number(label='RT (s)', format='%.1f').props('dense outlined').classes('w-24')
                    goto_mz = ui.number(label='m/z', format='%.2f').props('dense outlined').classes('w-24')
                    goto_spec = ui.number(label='Spectrum #', format='%.0f').props('dense outlined').classes('w-24')

                    def do_goto():
                        rt = goto_rt.value if goto_rt.value else None
                        mz = goto_mz.value if goto_mz.value else None
                        spec = int(goto_spec.value) if goto_spec.value else None
                        if rt is not None or mz is not None or spec is not None:
                            viewer.go_to_location(rt=rt, mz=mz, spectrum_idx=spec)
                            ui.notify(f"Navigated to location", type="info")

                    ui.button('Go', on_click=do_goto).props('dense color=primary')

        # FAIMS toggle (hidden by default, shown when FAIMS data is detected) and 3D View toggle
        with ui.row().classes('w-full justify-center gap-4 mb-2'):
            def toggle_faims_view():
                viewer.show_faims_view = faims_toggle.value
                if viewer.faims_container:
                    viewer.faims_container.set_visibility(viewer.show_faims_view)
                if viewer.df is not None and viewer.show_faims_view:
                    viewer.update_faims_plots()

            faims_toggle = ui.checkbox('FAIMS Multi-CV View', value=False, on_change=toggle_faims_view).classes('text-purple-400')
            faims_toggle.set_visibility(False)
            viewer.faims_toggle = faims_toggle

            # 3D View button
            def toggle_3d_view():
                viewer.show_3d_view = not viewer.show_3d_view
                if viewer.scene_3d_container:
                    viewer.scene_3d_container.set_visibility(viewer.show_3d_view)
                if viewer.show_3d_view and viewer.df is not None:
                    viewer.update_3d_view()
                # Update button appearance
                if viewer.show_3d_view:
                    view_3d_btn.props('color=purple')
                else:
                    view_3d_btn.props('color=grey')

            view_3d_btn = ui.button('3D Peak View', on_click=toggle_3d_view).props('dense color=grey').classes('text-sm')
            ui.label('(RT<120s, m/z<50)').classes('text-xs text-gray-500')

        # TIC Plot (clickable to show MS1 spectrum, zoomable to update peak map)
        with ui.card().classes('w-full max-w-6xl'):
            ui.label('TIC - Click to view spectrum, drag to zoom RT range').classes('text-xs text-gray-500 mb-1')
            viewer.tic_plot = ui.plotly(viewer.create_tic_plot()).classes('w-full')

            # Handle click on TIC plot - show closest spectrum and center peak map
            def on_tic_click(e):
                try:
                    if e.args and 'points' in e.args and e.args['points']:
                        point = e.args['points'][0]
                        if 'x' in point:
                            rt = point['x']
                            # Show closest spectrum (any MS level) - this also highlights in table
                            viewer.show_spectrum_at_rt(rt)
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

        # Main visualization area - peak map with spectrum browser overlay (collapsible)
        with ui.expansion('2D Peak Map', icon='grid_on', value=False).classes('w-full max-w-6xl'):
            # Display options row
            with ui.row().classes('w-full items-center gap-4 mb-2 flex-wrap'):
                ui.label('Overlay:').classes('text-xs text-gray-400')

                def toggle_centroids():
                    viewer.show_centroids = centroid_cb.value
                    if viewer.df is not None:
                        viewer.update_plot()

                centroid_cb = ui.checkbox('Centroids', value=True, on_change=toggle_centroids).props('dense').classes('text-green-400')

                def toggle_bboxes():
                    viewer.show_bounding_boxes = bbox_cb.value
                    if viewer.df is not None:
                        viewer.update_plot()

                bbox_cb = ui.checkbox('Bounding Boxes', value=False, on_change=toggle_bboxes).props('dense').classes('text-yellow-400')

                def toggle_hulls():
                    viewer.show_convex_hulls = hull_cb.value
                    if viewer.df is not None:
                        viewer.update_plot()

                hull_cb = ui.checkbox('Convex Hulls', value=False, on_change=toggle_hulls).props('dense').classes('text-cyan-400')

                def toggle_ids():
                    viewer.show_ids = ids_cb.value
                    if viewer.df is not None:
                        viewer.update_plot()

                ids_cb = ui.checkbox('Identifications', value=True, on_change=toggle_ids).props('dense').classes('text-orange-400')

                ui.label('|').classes('text-gray-600 mx-2')
                ui.label('Colormap:').classes('text-xs text-gray-400')

                def change_colormap(e):
                    viewer.colormap = e.value
                    if viewer.df is not None:
                        viewer.update_plot()
                        viewer.update_minimap()

                colormap_options = list(COLORMAPS.keys())
                ui.select(colormap_options, value='jet', on_change=change_colormap).props('dense outlined').classes('w-28')

            # Breadcrumb trail and coordinate display row
            with ui.row().classes('w-full items-center justify-between mb-1'):
                with ui.row().classes('items-center gap-2'):
                    ui.icon('navigation', size='xs').classes('text-gray-400')
                    viewer.breadcrumb_label = ui.label('Full view').classes('text-xs text-gray-400')
                viewer.coord_label = ui.label('RT: --  m/z: --').classes('text-xs text-cyan-400 font-mono')

            ui.label('Scroll to zoom, drag to select region, double-click to reset').classes('text-xs text-gray-500 mb-1')

            # Loading indicator overlay
            with ui.element('div').classes('relative'):
                viewer.loading_indicator = ui.element('div').classes(
                    'absolute inset-0 flex items-center justify-center bg-black/60 z-50'
                ).style('display: none;')
                with viewer.loading_indicator:
                    with ui.column().classes('items-center gap-2'):
                        ui.spinner('dots', size='xl', color='cyan')
                        viewer.loading_label = ui.label('Loading...').classes('text-cyan-400 text-sm font-mono')

            # Peak map with mouse interaction and minimap
            with ui.row().classes('w-full items-start gap-2'):
                # Peak map image with mouse handlers
                with ui.column().classes('flex-none'):
                    viewer.image_element = ui.image().classes('w-full').style(
                        f'width: {viewer.canvas_width}px; height: {viewer.canvas_height}px; background: #141419; cursor: crosshair;'
                    )

                    # Mouse wheel zoom handler
                    def on_wheel(e):
                        try:
                            # Get mouse position relative to image
                            offset_x = e.args.get('offsetX', 0)
                            offset_y = e.args.get('offsetY', 0)
                            delta_y = e.args.get('deltaY', 0)

                            # Convert to plot area coordinates (account for margins)
                            plot_x = offset_x - viewer.margin_left
                            plot_y = offset_y - viewer.margin_top

                            # Check if within plot area
                            if 0 <= plot_x <= viewer.plot_width and 0 <= plot_y <= viewer.plot_height:
                                x_frac = plot_x / viewer.plot_width
                                y_frac = plot_y / viewer.plot_height
                                zoom_in = delta_y < 0  # Scroll up = zoom in
                                viewer.zoom_at_point(x_frac, y_frac, zoom_in)
                        except Exception:
                            pass

                    viewer.image_element.on('wheel.prevent', on_wheel)

                    # Drag to select region for zoom (no continuous updates - just on release)
                    drag_state = {'dragging': False, 'start_x': 0, 'start_y': 0}

                    def on_mousedown(e):
                        offset_x = e.args.get('offsetX', 0)
                        offset_y = e.args.get('offsetY', 0)
                        # Only start drag if within plot area
                        plot_x = offset_x - viewer.margin_left
                        plot_y = offset_y - viewer.margin_top
                        if 0 <= plot_x <= viewer.plot_width and 0 <= plot_y <= viewer.plot_height:
                            drag_state['dragging'] = True
                            drag_state['start_x'] = offset_x
                            drag_state['start_y'] = offset_y

                    def on_mouseup(e):
                        if drag_state['dragging']:
                            drag_state['dragging'] = False
                            end_x = e.args.get('offsetX', 0)
                            end_y = e.args.get('offsetY', 0)

                            # Calculate selection in plot coordinates
                            start_plot_x = drag_state['start_x'] - viewer.margin_left
                            start_plot_y = drag_state['start_y'] - viewer.margin_top
                            end_plot_x = end_x - viewer.margin_left
                            end_plot_y = end_y - viewer.margin_top

                            # Ensure within bounds
                            start_plot_x = max(0, min(viewer.plot_width, start_plot_x))
                            start_plot_y = max(0, min(viewer.plot_height, start_plot_y))
                            end_plot_x = max(0, min(viewer.plot_width, end_plot_x))
                            end_plot_y = max(0, min(viewer.plot_height, end_plot_y))

                            # Only zoom if dragged a meaningful distance (>10 pixels)
                            dx = abs(end_plot_x - start_plot_x)
                            dy = abs(end_plot_y - start_plot_y)

                            if dx > 10 and dy > 10:
                                # Save current state to zoom history before changing
                                viewer.push_zoom_history()

                                # Convert to data coordinates
                                rt_range = viewer.view_rt_max - viewer.view_rt_min
                                mz_range = viewer.view_mz_max - viewer.view_mz_min

                                x1_frac = min(start_plot_x, end_plot_x) / viewer.plot_width
                                x2_frac = max(start_plot_x, end_plot_x) / viewer.plot_width
                                y1_frac = min(start_plot_y, end_plot_y) / viewer.plot_height
                                y2_frac = max(start_plot_y, end_plot_y) / viewer.plot_height

                                new_rt_min = viewer.view_rt_min + x1_frac * rt_range
                                new_rt_max = viewer.view_rt_min + x2_frac * rt_range
                                # Y is inverted (top = high m/z)
                                new_mz_max = viewer.view_mz_max - y1_frac * mz_range
                                new_mz_min = viewer.view_mz_max - y2_frac * mz_range

                                viewer.view_rt_min = new_rt_min
                                viewer.view_rt_max = new_rt_max
                                viewer.view_mz_min = new_mz_min
                                viewer.view_mz_max = new_mz_max

                                # Save new state to zoom history
                                viewer.push_zoom_history()
                                viewer.update_plot()

                    def on_dblclick(e):
                        viewer.reset_view()

                    # Mousemove handler for coordinate display
                    def on_mousemove(e):
                        try:
                            offset_x = e.args.get('offsetX', 0)
                            offset_y = e.args.get('offsetY', 0)
                            viewer.update_coord_display(offset_x, offset_y)
                        except Exception:
                            pass

                    def on_mouseleave_coord(e):
                        drag_state['dragging'] = False
                        if viewer.coord_label:
                            viewer.coord_label.set_text('RT: --  m/z: --')

                    viewer.image_element.on('mousedown', on_mousedown)
                    viewer.image_element.on('mouseup', on_mouseup)
                    viewer.image_element.on('mouseleave', on_mouseleave_coord)
                    viewer.image_element.on('mousemove', on_mousemove)
                    viewer.image_element.on('dblclick', on_dblclick)

                # Minimap panel (to the right of peak map)
                with ui.column().classes('flex-none'):
                    ui.label('Overview').classes('text-xs text-gray-400 mb-1')
                    viewer.minimap_image = ui.image().style(
                        f'width: {viewer.minimap_width}px; height: {viewer.minimap_height}px; '
                        f'background: #141419; cursor: pointer; border: 1px solid #333;'
                    )

                    # Minimap click handler
                    def on_minimap_click(e):
                        try:
                            offset_x = e.args.get('offsetX', 0)
                            offset_y = e.args.get('offsetY', 0)
                            x_frac = offset_x / viewer.minimap_width
                            y_frac = offset_y / viewer.minimap_height
                            viewer.minimap_click_to_view(x_frac, y_frac)
                        except Exception:
                            pass

                    viewer.minimap_image.on('click', on_minimap_click)

                    # Back button for zoom history
                    def go_back():
                        if len(viewer.zoom_history) > 1:
                            viewer.go_to_zoom_history(len(viewer.zoom_history) - 2)

                    ui.button('← Back', on_click=go_back).props('dense size=sm color=grey').classes('mt-1').tooltip('Go to previous view')

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

        # 3D View Container (hidden by default, shown when toggle enabled) - centered
        viewer.scene_3d_container = ui.column().classes('w-full max-w-6xl mt-2 mx-auto')
        viewer.scene_3d_container.set_visibility(False)
        with viewer.scene_3d_container:
            with ui.card().classes('w-full').style('background: #1a1a1f; padding: 0.5rem;'):
                with ui.row().classes('w-full items-center justify-between'):
                    with ui.row().classes('items-center gap-2'):
                        ui.label('3D Peak View').classes('text-sm font-semibold text-purple-400')
                        ui.label('(pyopenms-viz)').classes('text-xs text-gray-600')
                    viewer.view_3d_status = ui.label('').classes('text-xs text-yellow-400')

                # Create empty plotly figure for 3D view
                empty_fig = go.Figure()
                empty_fig.update_layout(
                    paper_bgcolor='#1a1a1f',
                    plot_bgcolor='#1a1a1f',
                    height=500,
                    margin=dict(l=0, r=0, t=0, b=0)
                )
                viewer.plot_3d = ui.plotly(empty_fig).classes('w-full').style('margin-top: -10px;')

        # Navigation controls
        with ui.row().classes('justify-center gap-2 mt-2'):
            ui.button('Reset View', on_click=viewer.reset_view).props('color=secondary')
            ui.button('Zoom In', on_click=lambda: viewer.zoom_in(0.5)).props('color=primary')
            ui.button('Zoom Out', on_click=lambda: viewer.zoom_out(2.0)).props('color=primary')
            ui.button('← Pan Left', on_click=lambda: viewer.pan(rt_frac=-0.25)).props('color=accent')
            ui.button('→ Pan Right', on_click=lambda: viewer.pan(rt_frac=0.25)).props('color=accent')
            ui.button('↑ Pan Up', on_click=lambda: viewer.pan(mz_frac=0.25)).props('color=accent')
            ui.button('↓ Pan Down', on_click=lambda: viewer.pan(mz_frac=-0.25)).props('color=accent')

        # Spectrum Table (moved up for better visibility)
        with ui.expansion('Spectrum Table', icon='list', value=True).classes('w-full max-w-6xl mt-2'):
            ui.label('Click a row to view the spectrum in the 1D viewer above').classes('text-sm text-gray-400 mb-2')

            # Filters row
            with ui.row().classes('w-full items-end gap-2 mb-2 flex-wrap'):
                ui.label('Filter:').classes('text-xs text-gray-400')
                spec_ms_filter = ui.select(['All', 'MS1', 'MS2'], value='All', label='MS Level').props('dense outlined').classes('w-20')
                spec_rt_min = ui.number(label='RT Min', format='%.0f').props('dense outlined').classes('w-20')
                spec_rt_max = ui.number(label='RT Max', format='%.0f').props('dense outlined').classes('w-20')
                spec_min_peaks = ui.number(label='Min Peaks', format='%.0f').props('dense outlined').classes('w-24')

                def apply_spectrum_filter():
                    ms_level = None
                    if spec_ms_filter.value == 'MS1':
                        ms_level = 1
                    elif spec_ms_filter.value == 'MS2':
                        ms_level = 2

                    filtered = viewer.filter_spectrum_data(
                        ms_level=ms_level,
                        rt_min=spec_rt_min.value if spec_rt_min.value else None,
                        rt_max=spec_rt_max.value if spec_rt_max.value else None,
                        min_peaks=int(spec_min_peaks.value) if spec_min_peaks.value else None
                    )
                    viewer.spectrum_table.update_rows(filtered)
                    ui.notify(f"Showing {len(filtered)} spectra", type="info")

                def reset_spectrum_filter():
                    spec_ms_filter.value = 'All'
                    spec_rt_min.value = None
                    spec_rt_max.value = None
                    spec_min_peaks.value = None
                    viewer.spectrum_table.update_rows(viewer.spectrum_data)

                ui.button('Apply', on_click=apply_spectrum_filter).props('dense size=sm color=primary')
                ui.button('Reset', on_click=reset_spectrum_filter).props('dense size=sm color=grey')

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
                pagination={'rowsPerPage': 10, 'sortBy': 'idx', 'descending': False},
                selection='single',
                on_select=lambda e: viewer.show_spectrum_in_browser(e.selection[0]['idx']) if e.selection else None
            ).classes('w-full').on('rowClick', on_spectrum_click)
            viewer.spectrum_table.props('dark flat bordered dense')

        # Feature Table
        with ui.expansion('Features', icon='scatter_plot').classes('w-full max-w-6xl mt-2'):
            ui.label('Click a row to zoom to that feature').classes('text-sm text-gray-400 mb-2')

            # Feature filters row
            with ui.row().classes('w-full items-end gap-2 mb-2 flex-wrap'):
                ui.label('Filter:').classes('text-xs text-gray-400')
                feat_min_intensity = ui.number(label='Min Intensity', format='%.0f').props('dense outlined').classes('w-28')
                feat_min_quality = ui.number(label='Min Quality', format='%.2f').props('dense outlined').classes('w-24')
                feat_charge = ui.select(['All', '1', '2', '3', '4', '5+'], value='All', label='Charge').props('dense outlined').classes('w-20')

                def apply_feature_filter():
                    charge_val = None
                    if feat_charge.value and feat_charge.value != 'All':
                        if feat_charge.value == '5+':
                            charge_val = 5  # Will match 5 or greater
                        else:
                            charge_val = int(feat_charge.value)

                    filtered = viewer.filter_feature_data(
                        min_intensity=feat_min_intensity.value if feat_min_intensity.value else None,
                        min_quality=feat_min_quality.value if feat_min_quality.value else None,
                        charge=charge_val
                    )
                    viewer.feature_table.update_rows(filtered)
                    ui.notify(f"Showing {len(filtered)} features", type="info")

                def reset_feature_filter():
                    feat_min_intensity.value = None
                    feat_min_quality.value = None
                    feat_charge.value = 'All'
                    viewer.feature_table.update_rows(viewer.feature_data)

                ui.button('Apply', on_click=apply_feature_filter).props('dense size=sm color=primary')
                ui.button('Reset', on_click=reset_feature_filter).props('dense size=sm color=grey')

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

            def on_feature_hover(e):
                """Handle feature row hover for visual feedback."""
                try:
                    row = e.args[1] if len(e.args) > 1 else None
                    if row and 'idx' in row:
                        viewer.set_hover_feature(row['idx'])
                except Exception:
                    pass

            def on_feature_leave(e):
                """Clear feature hover state."""
                viewer.clear_hover()

            viewer.feature_table = ui.table(
                columns=feature_columns, rows=[], row_key='idx',
                pagination={'rowsPerPage': 8, 'sortBy': 'intensity', 'descending': True}
            ).classes('w-full hover-highlight').on('rowClick', on_feature_click)
            viewer.feature_table.on('row-dblclick', on_feature_hover)  # Use dblclick as hover proxy
            viewer.feature_table.props('dark flat bordered dense')

        # ID Table
        with ui.expansion('Identifications', icon='biotech').classes('w-full max-w-6xl mt-2'):
            ui.label('Click a row to zoom and view annotated spectrum').classes('text-sm text-gray-400 mb-2')

            # ID filters row
            with ui.row().classes('w-full items-end gap-2 mb-2 flex-wrap'):
                ui.label('Filter:').classes('text-xs text-gray-400')
                id_seq_pattern = ui.input(label='Sequence', placeholder='e.g. PEPTIDE').props('dense outlined').classes('w-32')
                id_min_score = ui.number(label='Min Score', format='%.2f').props('dense outlined').classes('w-24')
                id_charge = ui.select(['All', '1', '2', '3', '4', '5+'], value='All', label='Charge').props('dense outlined').classes('w-20')

                def apply_id_filter():
                    charge_val = None
                    if id_charge.value and id_charge.value != 'All':
                        if id_charge.value == '5+':
                            charge_val = 5  # Will match 5 or greater
                        else:
                            charge_val = int(id_charge.value)

                    filtered = viewer.filter_id_data(
                        sequence_pattern=id_seq_pattern.value if id_seq_pattern.value else None,
                        min_score=id_min_score.value if id_min_score.value else None,
                        charge=charge_val
                    )
                    viewer.id_table.update_rows(filtered)
                    ui.notify(f"Showing {len(filtered)} identifications", type="info")

                def reset_id_filter():
                    id_seq_pattern.value = ''
                    id_min_score.value = None
                    id_charge.value = 'All'
                    viewer.id_table.update_rows(viewer.id_data)

                ui.button('Apply', on_click=apply_id_filter).props('dense size=sm color=primary')
                ui.button('Reset', on_click=reset_id_filter).props('dense size=sm color=grey')

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

            def on_id_hover(e):
                """Handle ID row hover for visual feedback."""
                try:
                    row = e.args[1] if len(e.args) > 1 else None
                    if row and 'idx' in row:
                        viewer.set_hover_id(row['idx'])
                except Exception:
                    pass

            viewer.id_table = ui.table(
                columns=id_columns, rows=[], row_key='idx',
                pagination={'rowsPerPage': 8, 'sortBy': 'score', 'descending': True}
            ).classes('w-full hover-highlight').on('rowClick', on_id_click)
            viewer.id_table.on('row-dblclick', on_id_hover)  # Use dblclick as preview
            viewer.id_table.props('dark flat bordered dense')

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
@click.option('--open/--no-open', '-o/-n', default=True, help='Open browser automatically (default: open)')
@click.option('--native', is_flag=True, default=False, help='Run as native desktop app (requires: pip install pywebview)')
def main(files, port, host, open, native):
    """
    mzML Peak Map Viewer - Fast visualization of mass spectrometry data.

    Pass one or more files to load them automatically:

    \b
    Examples:
        mzml_viewer.py sample.mzML                  # Open with browser
        mzml_viewer.py sample.mzML --native         # Open as desktop app
        mzml_viewer.py sample.mzML --no-open        # Start server only
        mzml_viewer.py sample.mzML features.featureXML
        mzml_viewer.py sample.mzML ids.idXML
        mzml_viewer.py data.mzML features.featureXML ids.idXML

    Supported file types (detected by extension):
        .mzML       Mass spectrometry peak data
        .featureXML Detected features with convex hulls
        .idXML      Peptide identifications

    Native mode (--native):
        Runs as a standalone desktop application using pywebview.
        Install with: pip install pywebview
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

    if native:
        click.echo(f"\nStarting native desktop app...")
    else:
        click.echo(f"\nStarting server at http://{host}:{port}")
        if open:
            click.echo("Opening browser...")

    @ui.page('/')
    def index():
        create_ui()

    ui.run(
        title='mzML Peak Map Viewer',
        host=host,
        port=port,
        reload=False,
        show=open and not native,
        native=native,
        window_size=(1400, 900) if native else None,
    )


if __name__ in {"__main__", "__mp_main__"}:
    main()
